from typing import Dict, List, Optional

import metatensor.torch as mts
import torch
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

    def __init__(self, cutoff, multiple_properties, feature_key="feature"):
        super().__init__()

        self._nl_request = NeighborListOptions(
            cutoff=cutoff, full_list=True, strict=True
        )
        self._multiple_properties = multiple_properties
        self._feature_key = feature_key

    def forward(
        self,
        systems: List[System],
        outputs: Dict[str, ModelOutput],
        selected_atoms: Optional[Labels],
    ) -> Dict[str, TensorMap]:
        if self._feature_key not in outputs:
            return {}

        device = torch.device("cpu")
        if len(systems) > 0:
            device = systems[0].positions.device

        output = outputs[self._feature_key]

        if output.sample_kind == "atom":
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
            assert output.sample_kind == "system"
            samples = Labels(
                "system", torch.arange(len(systems), device=device).reshape(-1, 1)
            )

        if self._multiple_properties:
            properties = Labels("cv", torch.tensor([[0], [1]], device=device))
        else:
            properties = Labels("cv", torch.tensor([[0]], device=device))

        values = torch.zeros(
            (len(samples), len(properties)), dtype=torch.float64, device=device
        )
        position_values: List[torch.Tensor] = []
        position_samples: List[List[int]] = []
        strain_values = torch.zeros(
            (len(samples), 3, 3, len(properties)), dtype=torch.float64, device=device
        )
        system_start = 0
        for system_i, system in enumerate(systems):
            system_stop = system_start + len(system)

            neighbors = system.get_neighbor_list(self._nl_request)

            atom_index = neighbors.samples.column("first_atom")
            distances = torch.linalg.vector_norm(neighbors.values.reshape(-1, 3), dim=1)
            inv_dist = 1.0 / distances

            if output.sample_kind == "atom":
                sliced = values[system_start:system_stop, 0]
                sliced += sliced.index_add(0, atom_index, inv_dist)
            else:
                values[system_i, 0] += inv_dist.sum()

            if self._multiple_properties:
                if output.sample_kind == "atom":
                    sliced = values[system_start:system_stop, 1]
                    sliced += sliced.index_add(0, atom_index, inv_dist**2)
                else:
                    values[system_i, 1] += inv_dist.sum() ** 2

            if len(output.explicit_gradients) != 0:
                # Derivatives of the same CVs as above, without relying on the
                # engine to differentiate a standard feature output.
                vectors = neighbors.values.reshape(-1, 3)
                second = neighbors.samples.column("second_atom")
                rows = len(system) if output.sample_kind == "atom" else 1
                derivatives = torch.zeros(
                    (rows, len(system), 3, len(properties)),
                    dtype=torch.float64,
                    device=device,
                )
                for pair in range(len(distances)):
                    a = int(atom_index[pair])
                    b = int(second[pair])
                    row = a if output.sample_kind == "atom" else 0
                    sample = (
                        system_start + row if output.sample_kind == "atom" else system_i
                    )
                    for prop in range(len(properties)):
                        factor = -(inv_dist[pair] ** 3)
                        if prop == 1:
                            if output.sample_kind == "atom":
                                factor = -2 * inv_dist[pair] ** 4
                            else:
                                factor = -2 * inv_dist.sum() * inv_dist[pair] ** 3
                        gradient = factor * vectors[pair]
                        derivatives[row, a, :, prop] -= gradient
                        derivatives[row, b, :, prop] += gradient
                        strain_values[sample, :, :, prop] += torch.outer(
                            vectors[pair], gradient
                        )
                position_values.append(derivatives.reshape(-1, 3, len(properties)))
                for row in range(rows):
                    sample = (
                        system_start + row if output.sample_kind == "atom" else system_i
                    )
                    for atom in range(len(system)):
                        position_samples.append([sample, system_i, atom])

            system_start = system_stop

        block = TensorBlock(
            values=values,
            samples=samples,
            components=[],
            properties=properties,
        )
        if "positions" in output.explicit_gradients:
            block.add_gradient(
                "positions",
                TensorBlock(
                    values=torch.cat(position_values)
                    if len(position_values) != 0
                    else torch.zeros(
                        (0, 3, len(properties)), dtype=torch.float64, device=device
                    ),
                    samples=Labels(
                        ["sample", "system", "atom"],
                        torch.tensor(
                            position_samples, dtype=torch.int32, device=device
                        ).reshape(-1, 3),
                    ),
                    components=[
                        Labels("xyz", torch.arange(3, dtype=torch.int32).reshape(-1, 1))
                    ],
                    properties=properties,
                ),
            )
        if "strain" in output.explicit_gradients:
            block.add_gradient(
                "strain",
                TensorBlock(
                    values=strain_values,
                    samples=Labels(
                        "sample",
                        torch.arange(len(samples), dtype=torch.int32).reshape(-1, 1),
                    ),
                    components=[
                        Labels(
                            "xyz_1", torch.arange(3, dtype=torch.int32).reshape(-1, 1)
                        ),
                        Labels(
                            "xyz_2", torch.arange(3, dtype=torch.int32).reshape(-1, 1)
                        ),
                    ],
                    properties=properties,
                ),
            )
        cv = TensorMap(
            keys=Labels("_", torch.tensor([[0]], device=device)),
            blocks=[block],
        )

        if selected_atoms is not None:
            if output.sample_kind == "atom":
                cv = mts.slice(cv, axis="samples", selection=selected_atoms)
            else:
                raise ValueError(
                    "selected atoms is only supported with per-atom output"
                )

        return {self._feature_key: cv}

    def requested_neighbor_lists(self) -> List[NeighborListOptions]:
        return [self._nl_request]


CUTOFF = 3.5

capabilities = ModelCapabilities(
    outputs={
        "feature": ModelOutput(
            sample_kind="atom", explicit_gradients=["positions", "strain"]
        )
    },
    interaction_range=CUTOFF,
    supported_devices=["cpu", "mps", "cuda"],
    length_unit="A",
    atomic_types=[6, 8],
    dtype="float64",
)

# export all variations of the model
cv = TestCollectiveVariable(cutoff=CUTOFF, multiple_properties=False)
cv.eval()
model = AtomisticModel(cv, ModelMetadata(), capabilities)
model.save("scalar-per-atom.pt")
# AtomisticModel adds deprecated aliases to this mutable capabilities object.
capabilities.outputs = {
    "feature": ModelOutput(
        sample_kind="atom", explicit_gradients=["positions", "strain"]
    )
}

cv = TestCollectiveVariable(cutoff=CUTOFF, multiple_properties=True)
cv.eval()
model = AtomisticModel(cv, ModelMetadata(), capabilities)
model.save("vector-per-atom.pt")

capabilities = ModelCapabilities(
    outputs={
        "feature": ModelOutput(
            sample_kind="system", explicit_gradients=["positions", "strain"]
        )
    },
    interaction_range=CUTOFF,
    supported_devices=["cpu", "mps", "cuda"],
    length_unit="A",
    atomic_types=[6, 8],
    dtype="float64",
)

cv = TestCollectiveVariable(cutoff=CUTOFF, multiple_properties=False)
cv.eval()
model = AtomisticModel(cv, ModelMetadata(), capabilities)
model.save("scalar-global.pt")
# AtomisticModel adds deprecated aliases to this mutable capabilities object.
capabilities.outputs = {
    "feature": ModelOutput(
        sample_kind="system", explicit_gradients=["positions", "strain"]
    )
}

cv = TestCollectiveVariable(cutoff=CUTOFF, multiple_properties=True)
cv.eval()
model = AtomisticModel(cv, ModelMetadata(), capabilities)
model.save("vector-global.pt")


cv = TestCollectiveVariable(
    cutoff=CUTOFF, multiple_properties=False, feature_key="feature/variant"
)
capabilities.outputs = {
    "feature/variant": ModelOutput(
        sample_kind="system", explicit_gradients=["positions", "strain"]
    )
}
cv.eval()
model = AtomisticModel(cv, ModelMetadata(), capabilities)
model.save("variant-global.pt")
