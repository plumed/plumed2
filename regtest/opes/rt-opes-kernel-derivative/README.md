# Kernel derivatives

The truncated kernels are force-shifted: with `c = exp(-rc2 / 2)` and `rc2` the
squared cutoff, the kernel value is `h * (exp(-r2 / 2) - c * (1 + (rc2 - r2) / 2))`,
and its derivative along component `i` is `-h * (exp(-r2 / 2) - c) * delta_i / sigma_i^2`.
Value and derivative both go to zero at the cutoff, and the forces are minus the
gradient of the bias.

This test deposits several kernels and compares Cartesian forces against
central finite differences of the bias energy, with updates disabled during
the comparison. Both OPES_METAD and OPES_METAD_EXPLORE are tested. A cutoff
of 3 makes the error of an inconsistent derivative large enough to detect
reliably: the test fails with the kernels used before plumed/plumed2#1463,
whose value was energy-shifted, `h * (exp(-r2 / 2) - c)`, while the derivative
was that of the force-shifted kernel.
