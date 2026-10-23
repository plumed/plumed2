# Python model for the metatomic-core plugin path: a scalar `feature` equal to
# the sum of the distances in a half pair list (i.e. the distance between the
# two atoms, including periodic images within the cutoff), with explicit
# gradients with respect to positions and strain.
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


CUTOFF = 4.0

# PLUMED stores pair lists in metatensor's C++ SimpleDataArray, which the Python
# metatensor package can only read once this origin is registered. This can be
# removed once metatensor registers it by default (registering twice is fine).
metatensor.register_external_data_wrapper(
    "metatensor::SimpleDataArray", metatensor.ExternalCpuArray
)


class TwoAtomDistance(NumpyModel):
    def capabilities(self):
        return ModelCapabilities(
            atomic_types=[1, 8],
            interaction_range=CUTOFF,
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
        return ModelMetadata(name="two-atom distance")

    def requested_pair_lists(self, requested_outputs):
        return [PairListOptions(cutoff=CUTOFF, full_list=False)]

    def execute_inner(self, systems, selected_atoms, requested_outputs):
        assert selected_atoms is None, "selected_atoms is not supported by this model"
        assert len(systems) == 1, "this model only supports a single system"

        system = systems[0]
        pairs = system.pairs(self.requested_pair_lists(requested_outputs)[0])
        vectors = np.asarray(pairs.values)[:, :, 0]
        first, second = np.asarray(pairs.samples.values)[:, :2].T

        distances = np.linalg.norm(vectors, axis=1)
        directions = vectors / distances[:, None]

        properties = Labels.range("distance", 1)
        block = TensorBlock(
            values=np.array([[distances.sum()]]),
            samples=Labels.range("system", 1),
            components=[],
            properties=properties,
        )

        gradients = requested_outputs[0].gradients
        if "positions" in gradients:
            positions = np.zeros((len(system), 3))
            np.add.at(positions, second, directions)
            np.add.at(positions, first, -directions)

            atoms = np.arange(len(system), dtype=np.int32)
            samples = np.stack([np.zeros_like(atoms), np.zeros_like(atoms), atoms], 1)
            block.add_gradient(
                "positions",
                TensorBlock(
                    values=positions[:, :, None],
                    samples=Labels(["sample", "system", "atom"], samples),
                    components=[Labels.range("xyz", 3)],
                    properties=properties,
                ),
            )

        if "strain" in gradients:
            strain = np.einsum("pa,pb->ab", vectors, directions)
            block.add_gradient(
                "strain",
                TensorBlock(
                    values=strain[None, :, :, None],
                    samples=Labels.range("sample", 1),
                    components=[Labels.range("xyz_1", 3), Labels.range("xyz_2", 3)],
                    properties=properties,
                ),
            )

        return [TensorMap(Labels.range("_", 1), [block])]


model = TwoAtomDistance()
