# OPES kernel merge stability

For positive kernel weights, merge the centered second moment:

```
w1 = h1 / (h1 + h2)
w2 = h2 / (h1 + h2)
delta = c2 - c1
c = c1 + w2 * delta
variance = w1 * sigma1^2 + w2 * sigma2^2 + w1 * w2 * delta^2
```

This is algebraically equivalent to the raw-second-moment formula but avoids
subtracting two large, nearly equal squared centers. The existing periodic
image selection and final wrapping are retained. No width floor or clamp is
introduced.

The regression uses identical and moving deposition centers, in both
OPES_METAD and OPES_METAD_EXPLORE. Adding a constant 100000000 to a nonperiodic
CV must preserve finite bias, forces and virial, and their translation
invariance. Narrow kernels with this offset lost their variance in the old
implementation; unshifted controls remain finite. This adversarial numerical
test is not evidence that any particular production failure has this cause.
