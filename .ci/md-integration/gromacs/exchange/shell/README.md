# Converged shell state control

This deliberately small polarization fixture checks evaluation order, not a
physical polarization parametrization. The tolerance 0.1 allows convergence
with GROMACS single precision. The zero-valued trace action observes distance,
ABMD minimum and auxiliary context without adding energy or forces.

`validate_shells.py` runs flat, same-kernel nested and explicitly loaded nested
PLUMED inputs. Shell coordinates must change, but the ratchet minimum must stay
fixed throughout each step's extra force evaluations. The first force evaluation
of each step remains ordinary. Separate basic regressions verify nonzero nested
forces and off-stride bias activation in both auxiliary modes.
