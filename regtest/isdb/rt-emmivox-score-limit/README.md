This reproducer uses one nitrogen atom and one map voxel. It checks the EMMIVOX
score and its coordinate derivative for positive, negative and zero residuals,
including points near zero and either side of the small-residual branch.

Run it with a PLUMED build that enables LibTorch:

```sh
plumed driver --plumed calibration.dat --ixyz traj.xyz
python3 check.py prepare
plumed driver --plumed plumed.dat --ixyz traj.xyz
python3 check.py check
```

The first run calibrates the forward-model density. The checker accounts for
float32 rounding and compares the actual action outputs with independent double
precision Simpson quadrature of a positive integrand. It normalizes coordinate
derivatives by a well-conditioned residual so the atom model cancels. It needs
only Python's standard library, not Python PyTorch or SciPy.

For z=deviation/sigma, the normalized score density is proportional to
S(z)=integral_0^1 exp(-z*z*t*t/2)dt and its derivative is z times the weighted
mean of t*t. The unfixed float32 erf/subtraction formula can produce zero,
wrong-sign or very large derivatives near z=0; the exact derivative is z/3 to
leading order. The score itself must also approach its finite zero-residual limit.
