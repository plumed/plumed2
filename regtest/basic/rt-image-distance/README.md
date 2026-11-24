IMAGE_DISTANCE validation
=========================

`main.cpp` uses the public MD interface and BIASVALUE so that values, atomic
forces and virials are tested through real force propagation. Central finite
differences use 1e-6 perturbations, 2e-7 absolute tolerance. Box tests deform
positions and every cell vector together, checking all nine components,
including shear. Tests avoid active-branch switches when differentiating.
The separate driver fixture exercises WHOLEMOLECULES and NUMERICAL_DERIVATIVES.

The oracle directly enumerates all ordered pairs and original-basis integer
shifts in [-12,12]^3 for analytic fixtures and [-34,34]^3 for random cases. This differs from the production reduced-basis adaptive
search. Explicit analytic fixtures include zero-shift exclusion, self images,
coherent translations, an extended displacement requiring shift -3, and a skew
basis whose shortest vector requires (-3,1,0). The production guarantee is the
reciprocal-bound proof in the action documentation, not a claim that this fixed
oracle works for arbitrary cells. All random cases use [-34,34]^3.

The random generator is deterministic. Its boxes have diagonal entries in [3,5]
and lower off-diagonal entries in [-3,3]; coordinates lie in [-8,8]. A conservative
bound for inverse column 0 is (1/3,1/3,8/9), column 1 is (0,1/3,1/3), column 2
is (0,0,1/3). Thus for any pair |s_k| <= 16 sum_l |inverse(l,k)| <= 224/9,
and R <= 5, so every improvement has |n_k| < 31. The wider cube of radius 34
therefore certifies all random comparisons. The arbitrary-cell proof remains
the adaptive reciprocal bound in the action documentation.

Repeated tie tests check determinism without finite differences at a cusp.
Near contact checks finite forces; zero contact deliberately throws. Input
checks include empty selection, missing/partially singular cell, nonfinite
coordinates and box, and integer-bound overflow. A 20-step gradient descent
under LOWER_WALLS checks that the propagated forces reduce the wall energy.
WHOLEMOLECULES connectivity order matters; once reconstructed, ATOMS order does
not affect values. Individual wrapping without reconstruction is outside the
documented coordinate convention.

To reproduce in a configured, built tree:

```
PLUMED_PROGRAM_NAME="$PWD/src/lib/plumed" make -C regtest/basic/rt-image-distance
PLUMED_PROGRAM_NAME="$PWD/src/lib/plumed" make -C regtest/basic/rt-image-distance-driver
cd regtest/basic/rt-image-distance/tmp
./exe --benchmark
```

The optional benchmark reports time/frame for 10 frames each at N=10,100,500,
1000. It uses compact random coordinates in a 4x5x6 cell. No timing thresholds
are used in regression tests. Setup is excluded; force-interface overhead is
included. Production cost depends on geometry and N squared. No MPI is needed.
