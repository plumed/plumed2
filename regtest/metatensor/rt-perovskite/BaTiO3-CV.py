from typing import Dict, List, Optional

import torch
from metatensor.torch import Labels, TensorBlock, TensorMap
from metatensor.torch.atomistic import (
    MetatensorAtomisticModel,
    ModelCapabilities,
    ModelMetadata,
    ModelOutput,
    System,
)
from rascaline.torch import SphericalExpansion


class BaTiO3_CV(torch.nn.Module):
    def __init__(self):
        super().__init__()

        # only compute the representation for λ=1, Ti centers and O neighbors
        self.selected_keys = Labels(
            ["o3_lambda", "o3_sigma", "center_type", "neighbor_type"],
            values=torch.tensor([[1, 1, 56, 8]]),
        )
        self.neighbors = Labels(
            ["neighbor_type"],
            values=torch.tensor([[8], [22], [56]]),
        )

        self.calculator = SphericalExpansion(
            cutoff=3.0,
            max_angular=1,
            max_radial=1,
            radial_basis={"Gto": {}},
            cutoff_function={"ShiftedCosine": {"width": 0.5}},
            center_atom_weight=1.0,
            atomic_gaussian_width=0.5,
        )

    def forward(
        self,
        systems: List[System],
        outputs: Dict[str, ModelOutput],
        selected_atoms: Optional[Labels],
    ) -> Dict[str, TensorMap]:

        if "plumed::cv" not in outputs:
            return {}

        spherical_expansion = self.calculator(
            systems,
            selected_keys=self.selected_keys,
            selected_samples=selected_atoms,
        )

        if len(systems[0]) == 0:
            # PLUMED is trying to determine the size of the output
            CV = torch.zeros((1, 1), dtype=torch.float64)
        else:
            CV = 100 * spherical_expansion.block().values.mean(dim=0).norm()

        block = TensorBlock(
            values=CV.reshape(1, 1),
            samples=Labels("system", torch.tensor([[0]])),
            components=[],
            properties=Labels("cv", torch.tensor([[0]])),
        )
        cv = TensorMap(
            keys=Labels("_", torch.tensor([[0]])),
            blocks=[block],
        )

        return {"plumed::cv": cv}


cv = BaTiO3_CV()
cv.eval()


capabilities = ModelCapabilities(
    outputs={"plumed::cv": ModelOutput(per_atom=False)},
    interaction_range=3.0,
    supported_devices=["cpu"],
    length_unit="A",
    atomic_types=[56, 22, 8],
    dtype="float64",
)

metadata = ModelMetadata(
    name="Polarisation-inspired CV in BaTiO3",
    description="""
This is an alternative implementation of the CV used in https://arxiv.org/abs/2310.12579
to study phase transition in BaTiO3.
""",
    references={
        "model": ["https://arxiv.org/abs/2310.12579"],
    },
)


model = MetatensorAtomisticModel(cv, metadata, capabilities)
model.export("BaTiO3_CV.pt", collect_extensions="extensions")
