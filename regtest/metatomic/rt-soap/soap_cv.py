from typing import Dict, List, Optional

import torch
from metatensor.torch import Labels, TensorBlock, TensorMap
from metatomic.torch import (
    AtomisticModel,
    ModelCapabilities,
    ModelMetadata,
    ModelOutput,
    System,
)
from featomic.torch import SoapPowerSpectrum


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
        if "features/missing" in outputs:
            # only there for testing purposes
            raise ValueError("This output is not supported")

        if "features/soap" not in outputs:
            return {}

        if not outputs["features/soap"].per_atom:
            raise ValueError("per_atom=False is not supported")

        if len(systems[0]) == 0:
            # PLUMED is trying to determine the size of the output
            projected = torch.zeros((0, 3), dtype=torch.float64)
            samples = Labels(["system", "atom"], torch.zeros((0, 2), dtype=torch.int32))
        else:
            soap = self.calculator(systems, selected_samples=selected_atoms)
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
        cv = TensorMap(
            keys=Labels("_", torch.tensor([[0]])),
            blocks=[block],
        )

        return {"features/soap": cv}


cv = SOAP_CV(species=[1, 6, 7, 8])
cv.eval()


capabilities = ModelCapabilities(
    outputs={
        "features/soap": ModelOutput(per_atom=True, description="SOAP PCA features"),
        "features/missing": ModelOutput(per_atom=True, description="Test features"),
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
