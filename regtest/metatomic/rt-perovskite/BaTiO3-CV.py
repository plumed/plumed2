from typing import Dict, List, Optional

import torch
from featomic.torch import SphericalExpansion
from metatensor.torch import Labels, TensorBlock, TensorMap
from metatomic.torch import (
    AtomisticModel,
    ModelCapabilities,
    ModelMetadata,
    ModelOutput,
    System,
)


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
            cutoff={
                "radius": 3.0,
                "smoothing": {"type": "ShiftedCosine", "width": 0.5},
            },
            density={"type": "Gaussian", "width": 0.5},
            basis={
                "type": "TensorProduct",
                "max_angular": 1,
                "radial": {"type": "Gto", "max_radial": 0},
            },
        )

    def forward(
        self,
        systems: List[System],
        outputs: Dict[str, ModelOutput],
        selected_atoms: Optional[Labels],
    ) -> Dict[str, TensorMap]:
        if "feature" not in outputs:
            return {}

        gradients: List[str] = []
        if outputs["feature"].explicit_gradients:
            gradients = ["positions", "cell"]
        spherical_expansion = self.calculator(
            systems,
            gradients=gradients,
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
        if gradients:
            positions = torch.zeros((len(systems[0]), 3), dtype=torch.float64)
            strain = torch.zeros((3, 3), dtype=torch.float64)
            if len(systems[0]) != 0:
                expansion = spherical_expansion.block()
                mean = expansion.values.mean(dim=0)
                weight = torch.zeros_like(mean)
                if mean.norm() > 0:
                    weight = 100 * mean / (mean.norm() * len(expansion.samples))
                gradient = expansion.gradient("positions")
                values = (gradient.values * weight).sum(dim=(2, 3))
                positions.index_add_(0, gradient.samples.column("atom"), values)
                cell_gradient = expansion.gradient("cell")
                cell_values = (cell_gradient.values * weight).sum(dim=(3, 4)).sum(dim=0)
                strain = (
                    systems[0].positions.T @ positions + systems[0].cell.T @ cell_values
                )
            if "positions" in outputs["feature"].explicit_gradients:
                samples = torch.zeros((len(systems[0]), 3), dtype=torch.int32)
                samples[:, 2] = torch.arange(len(systems[0]), dtype=torch.int32)
                block.add_gradient(
                    "positions",
                    TensorBlock(
                        positions[:, :, None],
                        Labels(["sample", "system", "atom"], samples),
                        [
                            Labels(
                                "xyz", torch.arange(3, dtype=torch.int32).reshape(-1, 1)
                            )
                        ],
                        block.properties,
                    ),
                )
            if "strain" in outputs["feature"].explicit_gradients:
                block.add_gradient(
                    "strain",
                    TensorBlock(
                        strain[None, :, :, None],
                        Labels("sample", torch.tensor([[0]])),
                        [
                            Labels(
                                "xyz_1",
                                torch.arange(3, dtype=torch.int32).reshape(-1, 1),
                            ),
                            Labels(
                                "xyz_2",
                                torch.arange(3, dtype=torch.int32).reshape(-1, 1),
                            ),
                        ],
                        block.properties,
                    ),
                )
        cv = TensorMap(
            keys=Labels("_", torch.tensor([[0]])),
            blocks=[block],
        )

        return {"feature": cv}


cv = BaTiO3_CV()
cv.eval()

capabilities = ModelCapabilities(
    outputs={
        "feature": ModelOutput(
            sample_kind="system", explicit_gradients=["positions", "strain"]
        )
    },
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


model = AtomisticModel(cv, metadata, capabilities)
model.save("BaTiO3_CV.pt", collect_extensions="extensions")
