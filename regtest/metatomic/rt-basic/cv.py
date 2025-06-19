from typing import Dict, List, Optional

import torch
import metatensor.torch as mts
from metatensor.torch import Labels, TensorBlock, TensorMap
from metatomic.torch import (
    AtomisticModel,
    ModelCapabilities,
    ModelMetadata,
    ModelOutput,
    NeighborListOptions,
    System,
)


class TestCollectiveVariable(torch.nn.Module):
    r"""
    This class computes a simple CV which is then used to test metatomic integration
    with PLUMED.

    The per-atom CV is defined as a sum over all pairs for an atom:

        CV^1_i = \sum_j 1/r_ij
        CV^2_i = \sum_j 1/r_ij^2

    The global CV is a sum over all atoms of the per-atom CV:

        CV^1 = \sum_i CV^1_i
        CV^2 = \sum_i CV^2_i

    If ``multiple_properties=True``, a only CV^1 is returned, otherwise both CV^1 and
    CV^2 are returned.
    """

    def __init__(self, cutoff, multiple_properties):
        super().__init__()

        self._nl_request = NeighborListOptions(
            cutoff=cutoff, full_list=True, strict=True
        )
        self._multiple_properties = multiple_properties

    def forward(
        self,
        systems: List[System],
        outputs: Dict[str, ModelOutput],
        selected_atoms: Optional[Labels],
    ) -> Dict[str, TensorMap]:
        if "features" not in outputs:
            return {}

        device = torch.device("cpu")
        if len(systems) > 0:
            device = systems[0].positions.device

        output = outputs["features"]

        if output.per_atom:
            samples_list: List[List[int]] = []
            for s, system in enumerate(systems):
                for i in range(len(system)):
                    samples_list.append([s, i])

            sample_values = torch.tensor(samples_list, device=device, dtype=torch.int32)
            samples = Labels(
                ["system", "atom"],
                sample_values.reshape(-1, 2),
            )
        else:
            samples = Labels(
                "system", torch.arange(len(systems), device=device).reshape(-1, 1)
            )

        if self._multiple_properties:
            properties = Labels("cv", torch.tensor([[0], [1]], device=device))
        else:
            properties = Labels("cv", torch.tensor([[0]], device=device))

        values = torch.zeros(
            (len(samples), len(properties)), dtype=torch.float32, device=device
        )
        system_start = 0
        for system_i, system in enumerate(systems):
            system_stop = system_start + len(system)

            neighbors = system.get_neighbor_list(self._nl_request)

            atom_index = neighbors.samples.column("first_atom")
            distances = torch.linalg.vector_norm(neighbors.values.reshape(-1, 3), dim=1)
            inv_dist = 1.0 / distances

            if output.per_atom:
                sliced = values[system_start:system_stop, 0]
                sliced += sliced.index_add(0, atom_index, inv_dist)
            else:
                values[system_i, 0] += inv_dist.sum()

            if self._multiple_properties:
                if output.per_atom:
                    sliced = values[system_start:system_stop, 1]
                    sliced += sliced.index_add(0, atom_index, inv_dist**2)
                else:
                    values[system_i, 1] += inv_dist.sum() ** 2

            system_start = system_stop

        block = TensorBlock(
            values=values,
            samples=samples,
            components=[],
            properties=properties,
        )
        cv = TensorMap(
            keys=Labels("_", torch.tensor([[0]], device=device)),
            blocks=[block],
        )

        if selected_atoms is not None:
            if output.per_atom:
                cv = mts.slice(cv, axis="samples", selection=selected_atoms)
            else:
                raise ValueError(
                    "selected atoms is only supported with per-atom output"
                )

        return {"features": cv}

    def requested_neighbor_lists(self) -> List[NeighborListOptions]:
        return [self._nl_request]


CUTOFF = 3.5

capabilities = ModelCapabilities(
    outputs={"features": ModelOutput(per_atom=True)},
    interaction_range=CUTOFF,
    supported_devices=["cpu", "mps", "cuda"],
    length_unit="A",
    atomic_types=[6, 8],
    dtype="float32",
)

# export all variations of the model
cv = TestCollectiveVariable(cutoff=CUTOFF, multiple_properties=False)
cv.eval()
model = AtomisticModel(cv, ModelMetadata(), capabilities)
model.save("scalar-per-atom.pt")

cv = TestCollectiveVariable(cutoff=CUTOFF, multiple_properties=True)
cv.eval()
model = AtomisticModel(cv, ModelMetadata(), capabilities)
model.save("vector-per-atom.pt")

capabilities = ModelCapabilities(
    outputs={"features": ModelOutput(per_atom=False)},
    interaction_range=CUTOFF,
    supported_devices=["cpu", "mps", "cuda"],
    length_unit="A",
    atomic_types=[6, 8],
    dtype="float32",
)

cv = TestCollectiveVariable(cutoff=CUTOFF, multiple_properties=False)
cv.eval()
model = AtomisticModel(cv, ModelMetadata(), capabilities)
model.save("scalar-global.pt")

cv = TestCollectiveVariable(cutoff=CUTOFF, multiple_properties=True)
cv.eval()
model = AtomisticModel(cv, ModelMetadata(), capabilities)
model.save("vector-global.pt")
