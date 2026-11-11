"""
Per-atom weighted sum of squared pair distances, with explicit gradients with
respect to positions and strain. The values are returned in the reverse order of
the requested atoms, to check that PLUMED follows the samples of the output.
"""

import metatensor
import numpy as np
from metatensor import Labels, TensorBlock, TensorMap

from metatomic import (
    ModelCapabilities,
    ModelMetadata,
    NumpyModel,
    PairListOptions,
    Quantity,
)


CUTOFF = 3.2

# PLUMED gives pair lists as metatensor's C++ SimpleDataArray
metatensor.register_external_data_wrapper(
    "metatensor::SimpleDataArray", metatensor.ExternalCpuArray
)


class PerAtomDistances(NumpyModel):
    def capabilities(self):
        return ModelCapabilities(
            atomic_types=[1],
            interaction_range=CUTOFF,
            length_unit="angstrom",
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
        return ModelMetadata(name="per-atom squared distances")

    def requested_pair_lists(self, requested_outputs):
        return [PairListOptions(cutoff=CUTOFF, full_list=False)]

    def execute_inner(self, systems, selected_atoms, requested_outputs):
        (system,) = systems
        pairs = system.pairs(self.requested_pair_lists(requested_outputs)[0])
        vectors = np.asarray(pairs.values)[:, :, 0]
        first, second = np.asarray(pairs.samples.values)[:, :2].T

        atoms = np.arange(len(system))
        if selected_atoms is not None:
            atoms = np.asarray(selected_atoms.values)[:, 1]
        atoms = atoms[::-1]

        values = np.zeros((len(atoms), 1))
        positions = np.zeros((len(atoms), len(system), 3, 1))
        strain = np.zeros((len(atoms), 3, 3, 1))
        for row, center in enumerate(atoms):
            weight = float(center + 1)
            for a, b, vector in zip(first, second, vectors):
                if center not in (a, b):
                    continue
                values[row, 0] += weight * vector @ vector
                positions[row, a, :, 0] -= 2 * weight * vector
                positions[row, b, :, 0] += 2 * weight * vector
                strain[row, :, :, 0] += 2 * weight * np.outer(vector, vector)

        properties = Labels.range("feature", 1)
        block = TensorBlock(
            values,
            Labels(
                ["system", "atom"],
                np.column_stack([np.zeros(len(atoms), dtype=np.int32), atoms]),
            ),
            [],
            properties,
        )

        gradients = requested_outputs[0].gradients
        if "positions" in gradients:
            samples = [
                [row, 0, atom] for row in range(len(atoms)) for atom in range(len(system))
            ]
            block.add_gradient(
                "positions",
                TensorBlock(
                    positions.reshape(-1, 3, 1),
                    Labels(["sample", "system", "atom"], np.array(samples, dtype=np.int32)),
                    [Labels.range("xyz", 3)],
                    properties,
                ),
            )
        if "strain" in gradients:
            block.add_gradient(
                "strain",
                TensorBlock(
                    strain,
                    Labels.range("sample", len(atoms)),
                    [Labels.range("xyz_1", 3), Labels.range("xyz_2", 3)],
                    properties,
                ),
            )
        return [TensorMap(Labels.range("_", 1), [block])]


model = PerAtomDistances()
