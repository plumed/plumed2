# Multi-simulation auxiliary evaluation controls

The MCMC input is a two-atom four-state expanded ensemble. With replica states
[0,1], PLUMED ENSEMBLE of the lambda indices and a restraint of strength2 give
local bias0.25. On the unguarded host, replica0's state1 candidate is evaluated
simultaneously with replica1's state2 candidate, giving mean1.5 and bias2.25.
The conditional configuration [1,1] instead gives mean1 and bias1.0. Independent
per-replica acceptance cannot use that simultaneous-candidate value.

The shell control reuses the common polarization fixture. A tiny initial shell
separation converges after the first force call, while the other simulation
needs more. ENSEMBLE remains active even with a zero-strength restraint, so an
extra force call makes a collective that the converged simulation will not make.
The aligned-iteration control completes on the old host; unequal counts hang.

The validator checks clear rejection of both unsupported multi-simulation paths
and successful ordinary replica-coupled MD. Single-simulation MCMC and shell
support, including multiple ranks and nested forces, have separate controls.
The host conservatively rejects independent multi-simulation use too because
custom bias couplings cannot be inferred reliably; separate mdrun invocations
are the supported alternative.
