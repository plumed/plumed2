# Models whose output indices are out of range, to check that PLUMED rejects
# them rather than writing outside its own buffers.
#
# Nothing upstream bounds these: metatomic-core validates the names and
# components of a gradient block, not the values of its samples, so the
# `sample` and `atom` columns below reach PLUMED unchecked even with
# CHECK_CONSISTENCY enabled.
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

metatensor.register_external_data_wrapper(
    "metatensor::SimpleDataArray", metatensor.ExternalCpuArray
)


class BrokenDistance(NumpyModel):
    """
    The rt-core distance model, with one index deliberately out of range. The
    `broken` option picks which one.
    """

    def set_options(self, options):
        broken = options.pop("broken")
        if broken not in ("atom", "names") or len(options) != 0:
            raise ValueError(f"unexpected options: broken={broken}, {options}")
        self._broken = broken

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
        return ModelMetadata(name="broken two-atom distance")

    def requested_pair_lists(self, requested_outputs):
        return [PairListOptions(cutoff=CUTOFF, full_list=False)]

    def execute_inner(self, systems, selected_atoms, requested_outputs):
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

        positions = np.zeros((len(system), 3))
        np.add.at(positions, second, directions)
        np.add.at(positions, first, -directions)

        atoms = np.arange(len(system), dtype=np.int32)
        rows = np.zeros_like(atoms)
        names = ["sample", "system", "atom"]

        if self._broken == "names":
            # the third column is no longer the atom
            names = ["sample", "system", "not_an_atom"]
        elif len(system) > 0:
            # PLUMED runs the model once on an empty system to discover how many
            # properties it returns, so that call is left alone.
            assert self._broken == "atom"
            # one past the last atom of the system
            atoms = atoms.copy()
            atoms[-1] = len(system)

        samples = np.stack([rows, np.zeros_like(atoms), atoms], 1)
        block.add_gradient(
            "positions",
            TensorBlock(
                values=positions[:, :, None],
                samples=Labels(names, samples),
                components=[Labels.range("xyz", 3)],
                properties=properties,
            ),
        )

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


model = BrokenDistance()
