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
