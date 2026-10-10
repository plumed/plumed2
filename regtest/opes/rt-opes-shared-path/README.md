# Shared ring-polymer OPES deposition weights

A configurable number of synchronized beads have distinct coordinates and optional bead-local walls.
Compare every KERNELS logweight with the independently measured pre-update
energy: local energy/kBT by default, full-path mean energy/kBT with
WALKERS_SHARED_BIAS. The four-rank native test uses two ranks per bead;
the executable accepts the bead count as its first argument (default 2).
The MPI size must be divisible by that count; any count of at least 2 is
supported, including odd counts. Coordinates span the same interval for
each layout so the test is not tied to the number of beads.

Frozen-field coordinate finite differences check the applied force, with and
without a wall. Matching-mode STATE restart must work; both directions of a
local/shared-path history change must fail. These are estimator and interface
tests, not evidence of equilibrium or long-time numerical stability.

An unsupported adaptive width or excluded region on only one walker must be
rejected collectively, without stranding other walkers inside an MPI call.
