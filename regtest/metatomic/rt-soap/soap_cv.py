from typing import Dict, List, Optional

import torch
from featomic.torch import SoapPowerSpectrum
from metatensor.torch import Labels, TensorBlock, TensorMap
from metatomic.torch import (
    AtomisticModel,
    ModelCapabilities,
    ModelMetadata,
    ModelOutput,
    System,
)


class SOAP_CV(torch.nn.Module):
    def __init__(self, species):
        super().__init__()

        self.neighbor_type_pairs = Labels(
            names=["neighbor_1_type", "neighbor_2_type"],
            values=torch.tensor(
                [[t1, t2] for t1 in species for t2 in species if t1 <= t2]
            ),
        )
        self.calculator = SoapPowerSpectrum(
            cutoff={
                "radius": 4.0,
                "smoothing": {"type": "ShiftedCosine", "width": 0.5},
            },
            density={"type": "Gaussian", "width": 0.3},
            basis={
                "type": "TensorProduct",
                "max_angular": 6,
                "radial": {"type": "Gto", "max_radial": 5},
            },
        )

        torch.manual_seed(-230623)
        self.register_buffer("pca_projection", torch.rand(2520, 3, dtype=torch.float64))

    def forward(
        self,
        systems: List[System],
        outputs: Dict[str, ModelOutput],
        selected_atoms: Optional[Labels],
    ) -> Dict[str, TensorMap]:
        if "feature" not in outputs:
            return {}

        if outputs["feature"].sample_kind != "atom":
            raise ValueError("only per-atom output is supported")

        soap_block: Optional[TensorBlock] = None
        if len(systems[0]) == 0:
            # PLUMED is trying to determine the size of the output
            projected = torch.zeros((0, 3), dtype=torch.float64)
            samples = Labels(["system", "atom"], torch.zeros((0, 2), dtype=torch.int32))
        else:
            gradients: List[str] = []
            if outputs["feature"].explicit_gradients:
                gradients = ["positions", "cell"]
            soap = self.calculator(
                systems, selected_samples=selected_atoms, gradients=gradients
            )
            soap = soap.keys_to_samples("center_type")
            soap = soap.keys_to_properties(self.neighbor_type_pairs)

            soap_block = soap.block()
            projected = soap_block.values @ self.pca_projection

            samples = soap_block.samples.remove("center_type")

        block = TensorBlock(
            values=projected,
            samples=samples,
            components=[],
            properties=Labels("soap_pca", torch.tensor([[0], [1], [2]])),
        )
        for parameter in outputs["feature"].explicit_gradients:
            if soap_block is None:
                components = [
                    Labels("xyz", torch.arange(3, dtype=torch.int32).reshape(-1, 1))
                ]
                gradient_samples = Labels(
                    ["sample", "system", "atom"], torch.zeros((0, 3), dtype=torch.int32)
                )
                gradient_values = torch.zeros((0, 3, 3), dtype=torch.float64)
                if parameter == "strain":
                    components = [
                        Labels(
                            "xyz_1", torch.arange(3, dtype=torch.int32).reshape(-1, 1)
                        ),
                        Labels(
                            "xyz_2", torch.arange(3, dtype=torch.int32).reshape(-1, 1)
                        ),
                    ]
                    gradient_samples = Labels(
                        "sample", torch.arange(0, dtype=torch.int32).reshape(-1, 1)
                    )
                    gradient_values = torch.zeros((0, 3, 3, 3), dtype=torch.float64)
            elif parameter == "positions":
                gradient = soap_block.gradient("positions")
                components = gradient.components
                gradient_samples = gradient.samples
                gradient_values = gradient.values @ self.pca_projection
            else:
                # featomic exposes cell derivatives at fixed positions. Strain
                # deforms both: R -> R (I + strain), H -> H (I + strain).
                cell_gradient = soap_block.gradient("cell")
                cell_values = cell_gradient.values @ self.pca_projection
                gradient_values = torch.einsum(
                    "ka,skbp->sabp", systems[0].cell, cell_values
                )
                positions_gradient = soap_block.gradient("positions")
                positions_values = positions_gradient.values @ self.pca_projection
                for row in range(len(positions_gradient.samples)):
                    sample = int(positions_gradient.samples.values[row, 0])
                    atom = int(positions_gradient.samples.values[row, 2])
                    gradient_values[sample] += torch.einsum(
                        "a,bp->abp", systems[0].positions[atom], positions_values[row]
                    )
                components = [
                    Labels("xyz_1", torch.arange(3, dtype=torch.int32).reshape(-1, 1)),
                    Labels("xyz_2", torch.arange(3, dtype=torch.int32).reshape(-1, 1)),
                ]
                gradient_samples = cell_gradient.samples
            block.add_gradient(
                parameter,
                TensorBlock(
                    values=gradient_values,
                    samples=gradient_samples,
                    components=components,
                    properties=block.properties,
                ),
            )
        cv = TensorMap(
            keys=Labels("_", torch.tensor([[0]])),
            blocks=[block],
        )

        return {"feature": cv}


cv = SOAP_CV(species=[1, 6, 7, 8])
cv.eval()


capabilities = ModelCapabilities(
    outputs={
        "feature": ModelOutput(
            sample_kind="atom", explicit_gradients=["positions", "strain"]
        )
    },
    interaction_range=4.0,
    supported_devices=["cpu"],
    length_unit="nm",
    atomic_types=[6, 1, 7, 8],
    dtype="float64",
)

metadata = ModelMetadata(name="Collective Variable test")
model = AtomisticModel(cv, metadata, capabilities)
model.save("soap_cv.pt", collect_extensions="extensions")
