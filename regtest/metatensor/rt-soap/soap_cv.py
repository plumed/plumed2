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
from rascaline.torch import SoapPowerSpectrum


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
            cutoff=4.0,
            max_angular=6,
            max_radial=6,
            radial_basis={"Gto": {}},
            cutoff_function={"ShiftedCosine": {"width": 0.5}},
            center_atom_weight=1.0,
            atomic_gaussian_width=0.3,
        )

        self.pca_projection = torch.rand(2520, 3, dtype=torch.float64)

    def forward(
        self,
        systems: List[System],
        outputs: Dict[str, ModelOutput],
        selected_atoms: Optional[Labels],
    ) -> Dict[str, TensorMap]:

        if "collective_variable" not in outputs:
            return {}

        output = outputs["collective_variable"]

        soap = self.calculator(systems, selected_samples=selected_atoms)
        soap = soap.keys_to_samples("center_type")
        soap = soap.keys_to_properties(self.neighbor_type_pairs)

        if not output.per_atom:
            raise ValueError("per_atom=False is not supported")

        soap_block = soap.block()
        projected = soap_block.values @ self.pca_projection

        block = TensorBlock(
            values=projected,
            samples=soap_block.samples,
            components=[],
            properties=Labels("soap_pca", torch.tensor([[0], [1], [2]])),
        )
        cv = TensorMap(keys=Labels("_", torch.tensor([[0]])), blocks=[block])

        return {"collective_variable": cv}


cv = SOAP_CV(species=[1, 6, 7, 8])
cv.eval()


capabilites = ModelCapabilities(
    outputs={
        "collective_variable": ModelOutput(
            quantity="",
            unit="",
            per_atom=True,
            explicit_gradients=["postions"],
        )
    },
    interaction_range=4.0,
    supported_devices=["cpu", "cuda"],
    length_unit="nm",
    atomic_types=[6, 1, 7, 8],
    # dtype=TODO
)

metadata = ModelMetadata(
    name="Collective Variable test",
    description="""
A simple collective variable for testing purposes
""",
    authors=["..."],
    references={
        "implementation": ["ref to SOAP code"],
        "architecture": ["ref to SOAP"],
        "model": ["ref to paper"],
    },
)


model = MetatensorAtomisticModel(cv, metadata, capabilites)
model.export("soap_cv.pt", collect_extensions="extensions")
