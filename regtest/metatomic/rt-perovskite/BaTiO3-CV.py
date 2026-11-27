"""
Alternative implementation of the CV used in https://arxiv.org/abs/2310.12579 to study
the phase transitions of BaTiO3: the norm of the mean λ=1 spherical expansion of the
oxygen density around the barium atoms, computed with featomic, with explicit gradients
with respect to positions and strain.
"""

import featomic
import numpy as np
from metatensor import Labels, TensorBlock, TensorMap

from metatomic import ModelCapabilities, ModelMetadata, NumpyModel, Quantity


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


class BaTiO3Polarization(NumpyModel):
    def __init__(self):
        super().__init__()
        self.calculator = featomic.SphericalExpansion(
            cutoff={"radius": 3.0, "smoothing": {"type": "ShiftedCosine", "width": 0.5}},
            density={"type": "Gaussian", "width": 0.5},
            basis={
                "type": "TensorProduct",
                "max_angular": 1,
                "radial": {"type": "Gto", "max_radial": 0},
            },
        )
        # λ=1, barium centers and oxygen neighbors
        self.selected_keys = Labels(
            ["o3_lambda", "o3_sigma", "center_type", "neighbor_type"],
            np.array([[1, 1, 56, 8]]),
        )

    def capabilities(self):
        return ModelCapabilities(
            atomic_types=[56, 22, 8],
            interaction_range=3.0,
            length_unit="A",
            supported_devices=["cpu"],
            dtype="float64",
            outputs=[
                Quantity(
                    name="feature",
                    unit="",
                    sample_kind="system",
                    gradients=["positions", "strain"],
                )
            ],
        )

    def metadata(self):
        return ModelMetadata(
            name="Polarisation-inspired CV in BaTiO3",
            references={"model": ["https://arxiv.org/abs/2310.12579"]},
        )

    def execute_inner(self, systems, selected_atoms, requested_outputs):
        (system,) = systems
        (output,) = requested_outputs
        properties = Labels.range("cv", 1)

        # CV = 100 |m|, with m the mean of the expansion v_i over the centers, and
        # dCV/dx = sum_i w . dv_i/dx with w = 100 m / (|m| n_centers)
        cv = 0.0
        weight = np.zeros((3, 1))
        expansion = None
        if len(system) != 0:
            expansion = self.calculator.compute(
                FeatomicSystem(system),
                gradients=list(output.gradients),
                selected_keys=self.selected_keys,
                selected_samples=selected_atoms,
            ).block()
            mean = expansion.values.mean(axis=0)
            cv = 100 * np.linalg.norm(mean)
            if cv > 0:
                weight = 100 * mean / (np.linalg.norm(mean) * len(expansion.samples))

        block = TensorBlock(np.array([[cv]]), Labels.range("system", 1), [], properties)

        if "positions" in output.gradients:
            positions = np.zeros((len(system), 3))
            if expansion is not None:
                gradient = expansion.gradient("positions")
                np.add.at(
                    positions,
                    gradient.samples.column("atom"),
                    np.einsum("gxmp,mp->gx", gradient.values, weight),
                )
            atoms = np.arange(len(system), dtype=np.int32)
            block.add_gradient(
                "positions",
                TensorBlock(
                    positions[:, :, None],
                    Labels(
                        ["sample", "system", "atom"],
                        np.column_stack([0 * atoms, 0 * atoms, atoms]),
                    ),
                    [Labels.range("xyz", 3)],
                    properties,
                ),
            )

        if "strain" in output.gradients:
            strain = np.zeros((3, 3))
            if expansion is not None:
                gradient = expansion.gradient("strain").values
                strain = np.einsum("sabmp,mp->ab", gradient, weight)
            block.add_gradient(
                "strain",
                TensorBlock(
                    strain[None, :, :, None],
                    Labels.range("sample", 1),
                    [Labels.range("xyz_1", 3), Labels.range("xyz_2", 3)],
                    properties,
                ),
            )

        return [TensorMap(Labels.range("_", 1), [block])]


model = BaTiO3Polarization()
