# Shifted Gaussian kernel derivatives

The kernel value is `h * (exp(-r2 / 2) - epsilon)`. The cutoff shift is
constant, so its derivative is zero. For component `i`, the derivative is
`-h * exp(-r2 / 2) * delta_i / sigma_i^2`; multiplying by the shifted kernel
value instead introduces a spurious cutoff-dependent term.

This test deposits several kernels and compares Cartesian forces against
central finite differences of the bias energy, with updates disabled during
the comparison. Both OPES_METAD and OPES_METAD_EXPLORE are tested. A cutoff
of 3 makes the old derivative error large enough to detect reliably.

Force references for existing fixed-trajectory OPES tests are updated for
the corrected derivative. Their energies, kernel files, inputs and numerical
tolerances are unchanged. The regression fails with the unmodified upstream
kernel and passes after this correction.
