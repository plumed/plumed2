"""
Per-atom CV: a random projection of the SOAP power spectrum, computed with featomic,
with explicit gradients with respect to positions and strain.
"""

import featomic
import numpy as np
from metatensor import Labels, TensorBlock, TensorMap

from metatomic import ModelCapabilities, ModelMetadata, NumpyModel, Quantity


SPECIES = [1, 6, 7, 8]


class FeatomicSystem(featomic.systems.SystemBase):
    """A metatomic system given to featomic, which computes the neighbors itself"""

    def __init__(self, system):
        super().__init__()
        if np.any(system.pbc) and not np.all(system.pbc):
            raise ValueError("featomic needs all directions to be periodic, or none")
        self._types = np.asarray(system.types, dtype=np.int32)
        self._positions = np.asarray(system.positions, dtype=np.float64)
        self._cell = np.asarray(system.cell, dtype=np.float64)

    def size(self):
        return len(self._types)

    def types(self):
        return self._types

    def positions(self):
        return self._positions

    def cell(self):
        return self._cell

    def compute_neighbors(self, cutoff):
        raise NotImplementedError("featomic computes the neighbors of this system")

    def pairs(self):
        raise NotImplementedError("featomic computes the neighbors of this system")

    def pairs_containing(self, atom):
        raise NotImplementedError("featomic computes the neighbors of this system")


class SoapProjection(NumpyModel):
    def __init__(self):
        super().__init__()
        self.calculator = featomic.SoapPowerSpectrum(
            cutoff={"radius": 4.0, "smoothing": {"type": "ShiftedCosine", "width": 0.5}},
            density={"type": "Gaussian", "width": 0.3},
            basis={
                "type": "TensorProduct",
                "max_angular": 6,
                "radial": {"type": "Gto", "max_radial": 5},
            },
        )
        self.neighbor_pairs = Labels(
            ["neighbor_1_type", "neighbor_2_type"],
            np.array([[a, b] for a in SPECIES for b in SPECIES if a <= b]),
        )
        self.projection = np.random.default_rng(230623).random((2520, 3))

    def capabilities(self):
        return ModelCapabilities(
            atomic_types=[6, 1, 7, 8],
            interaction_range=4.0,
            length_unit="nm",
            supported_devices=["cpu"],
            dtype="float64",
            outputs=[
                Quantity(
                    name="feature",
                    unit="",
                    sample_kind="atom",
                    gradients=["positions", "strain"],
                )
            ],
        )

    def metadata(self):
        return ModelMetadata(name="SOAP projection")

    def execute_inner(self, systems, selected_atoms, requested_outputs):
        (output,) = requested_outputs
        properties = Labels.range("soap_pca", 3)
        if all(len(system) == 0 for system in systems):
            # PLUMED is asking for the shape of the output
            empty = Labels(["system", "atom"], np.zeros((0, 2), dtype=np.int32))
            block = TensorBlock(np.zeros((0, 3)), empty, [], properties)
            return [TensorMap(Labels.range("_", 1), [block])]

        soap = self.calculator.compute(
            [FeatomicSystem(system) for system in systems],
            gradients=list(output.gradients),
            selected_samples=selected_atoms,
        )
        soap = soap.keys_to_samples("center_type").keys_to_properties(
            self.neighbor_pairs
        )
        soap = soap.block()

        block = TensorBlock(
            soap.values @ self.projection,
            soap.samples.remove("center_type"),
            [],
            properties,
        )
        for parameter, gradient in soap.gradients():
            block.add_gradient(
                parameter,
                TensorBlock(
                    gradient.values @ self.projection,
                    gradient.samples,
                    gradient.components,
                    properties,
                ),
            )
        return [TensorMap(Labels.range("_", 1), [block])]


model = SoapProjection()
