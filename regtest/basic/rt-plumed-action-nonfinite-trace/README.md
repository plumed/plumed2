# Opt-in non-finite action trace regression

`PLUMED_NONFINITE_ACTION_TRACE` accepts comma- or space-separated action labels,
or `*`. It is sampled once per process. With no selection the diagnostic is
disabled and does not change the force calculation.

After a selected action calculates, the diagnostic checks its current values
and stored scalar derivatives. After every active action applies, it checks
selected scalar forces and interface force buffers, including atomic position
and box channels even when their internal labels are not selected. Messages
identify the step, phase, triggering action, target value and element.

The regression covers a singular CUSTOM derivative, a singular atomic COLVAR
derivative, finite scalar inputs whose product overflows an atomic force,
overflow of an intermediate scalar force, overflow of a box force with finite
atomic forces, and a finite `select` expression at its branch boundary.
Synthetic singularities are expected failures, not supported physical inputs.

This is not an exhaustive sanitizer: grid-derivative layouts, hidden action
state and MD-engine-side unit conversion or force accumulation are outside
the scalar-derivative checks. A finite trace does not establish scientific
validity or guarantee that a later calculation will remain finite.

All spatial MPI ranks must use the same selection. The first failing rank
is reported to all participants before raising an exception.
