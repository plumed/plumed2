"""Per-atom squared pair distances with explicit position and strain gradients."""

import builtins
import sys

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

metatensor.register_external_data_wrapper(
    "metatensor::SimpleDataArray", metatensor.ExternalCpuArray
)


class PerAtomDistances(NumpyModel):
    def capabilities(self):
        return ModelCapabilities(
            atomic_types=[1],
            interaction_range=3.2,
            length_unit="angstrom",
            supported_devices=["cpu"],
            dtype="float64",
            outputs=[
                Quantity(
                    name="feature",
                    unit="",
                    sample_kind="atom",
                    gradients=["positions", "strain"],
                ),
                Quantity(
                    name="feature/no_positions",
                    unit="",
                    sample_kind="atom",
                    gradients=["strain"],
                ),
                Quantity(
                    name="feature/no_strain",
                    unit="",
                    sample_kind="atom",
                    gradients=["positions"],
                ),
            ],
        )

    def metadata(self):
        return ModelMetadata(name="per-atom explicit-gradient regression")

    def requested_pair_lists(self, requested_outputs):
        return [PairListOptions(cutoff=3.2, full_list=False)]

    def execute_inner(self, systems, selected_atoms, requested_outputs):
        # The optional pycv smoke test checks this marker from its own callback.
        builtins._metatomic_test_interpreter = id(sys.modules)
        (system,) = systems
        pairs = system.pairs(self.requested_pair_lists(requested_outputs)[0])
        vectors = np.asarray(pairs.values)[:, :, 0]
        first, second = np.asarray(pairs.samples.values)[:, :2].T
        atoms = (
            np.arange(len(system))
            if selected_atoms is None
            else np.asarray(selected_atoms.values)[:, 1]
        )
        # Deliberately return a different order from the requested selection.
        atoms = atoms[::-1]
        values = np.zeros((len(atoms), 1))
        positions = np.zeros((len(atoms), len(system), 3, 1))
        strain = np.zeros((len(atoms), 3, 3, 1))
        for row, center in enumerate(atoms):
            weight = float(center + 1)
            for a, b, vector in zip(first, second, vectors):
                if center != a and center != b:
                    continue
                values[row, 0] += weight * np.dot(vector, vector)
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
        if "positions" in requested_outputs[0].gradients:
            samples = [
                [row, 0, atom]
                for row in range(len(atoms))
                for atom in range(len(system))
            ]
            block.add_gradient(
                "positions",
                TensorBlock(
                    positions.reshape(-1, 3, 1),
                    Labels(
                        ["sample", "system", "atom"], np.array(samples, dtype=np.int32)
                    ),
                    [Labels.range("xyz", 3)],
                    properties,
                ),
            )
        if "strain" in requested_outputs[0].gradients:
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
