# Frozen complete-path functions

This optional module provides stateless functions for fixed complete-path
biases. It does not perform learning, define a PIMD integrator, or supply a
reference distribution. PATH_LOGMEANEXP provides stable replica-root
aggregation of probability ratios, and PROBABILITY_MIX combines frozen
path potentials with a global partition-function ratio. CONDITIONAL_PATH
uses a distinct coordinate-dependent normalizer. Consult each action for
its derivative, domain, synchronization and stationary-sampling contract.

The module is distributed under the GNU Lesser General Public License,
version 3 or later, like PLUMED. Enable it with `--enable-modules=+pathbias`.

## Terminology and composition

The order of averaging determines the path potential:

- **Coordinate-centroid bias:** evaluate the CV on the mean bead coordinates,
  then evaluate the bias.
- **Bead-averaged CV bias:** average the bead CVs with ENSEMBLE, then evaluate
  the bias.
- **Bead-averaged bias energy:** evaluate one common field on each bead, then
  average its energies. This is the LAMMPS `bead_density` convention.
- **Mean probability-ratio bias:** average `exp(-v_b/kBT)`, then take `-kBT`
  times the logarithm. PATH_LOGMEANEXP supplies the stable logarithm and its
  derivatives; it does not average the bias energies or an estimated density.

PROBABILITY_MIX combines two frozen path potentials using a global
normalizer. CONDITIONAL_PATH uses a coordinate-dependent normalizer with
its own derivative. Their `COUPLING` parameters do not contract coordinates.
LAMMPS `path_contraction` instead changes the geometry supplied to the CV;
its intermediate values are generally not a linear mixture of endpoint CVs.
These descriptions leave existing action and keyword names unchanged.
