# Explicit auxiliary evaluation context

`setEvaluationMode` (command API 14) provides context before dependency
preparation. Send an integer after initialization and before `prepareCalc` or
`prepareDependencies`:

- 0: ordinary evaluation, with existing calc/update semantics.
- 1: auxiliary force evaluation. Histories are frozen; normal force STRIDE
  activation and scaling remain in effect. Use `performCalcNoUpdate`.
- 2: full auxiliary potential. Histories are frozen and all bias pilots are
  activated, including off-stride biases. Use `performCalcNoForces`.

Restore mode 0 after the evaluation, including error paths. Nonzero modes reject
`performCalc` and `update`; mode 1 rejects `performCalcNoForces`, and mode 2
rejects `performCalcNoUpdate`. `clear` resets the mode. The mode remains in force
until explicitly changed, so preparation and nested kernels see the same context.
It does not imply foreign coordinates: `setExchangeStep` remains separate.
Modes 1 and 2 are local auxiliary evaluations: mode 1 allows shell coordinates
to move with the ordinary neighbor-list approximation; mode 2 reuses the current
ordinary coordinates while varying lambda/ENERGY. A caller supplying foreign
coordinates must additionally set the exchange flag before preparation.
`Action::isGeometryTrial()` distinguishes such foreign trials (and unclassified
mode-0 no-force calls) from local auxiliary calls. Geometry-list refresh and
reduced-list restrictions use this predicate; persistent-state protection still
uses `isCalculationOnly()`. Lambda calls can reuse initialized current-coordinate
lists even on an ordinary refresh step. Local repeated path calculations retain
the frame membership already selected for this step.
`isCalculationOnly()` means state changes are forbidden; it does not mean that
forces are unnecessary. Nested evaluations propagate these two requirements
independently, and force-bearing loaded guests require the new command.

The 2026 GROMACS patch uses mode 1 for repeated force calls at the same step
(shell relaxation), and mode 2 for lambda MCMC candidates. Older kernels reject
lambda MCMC with a clear message. Shell forces retain the legacy call sequence
with one main-rank GROMACS warning per provider: the old kernel cannot protect
adaptive state or check unsafe inter-replica communication during these calls.
Ordinary MD and common bias
exchange/HREX do not require the new command. Command API 14 identifies an actual
new capability, not the presence of the action fixes: maintenance backports of
those fixes retain their API numbers. Existing normal-mode calls, including
repeated full calculations at the same numerical step, keep their semantics.

The host supplies context explicitly rather than inferring it from the numerical
step: an existing caller can intentionally repeat full calculations at one step.

External integrations that reject API versions above a fixed maximum need their
accepted range updated to include 14. The maintained LAMMPS integration patch
includes this change; its existing MD commands and MPI injection are unchanged.

The GROMACS implementation restricts lambda MCMC and repeated shell forces to
one simulation per invocation (which can contain multiple ranks). Simultaneously
changing each replica's candidate does not evaluate the one-replica conditional
potential of a coupled model. Unequal shell iteration counts can mismatch
inter-replica collective calls. Ordinary multi-simulation MD and common exchange
are unaffected. This conservative host restriction also rejects independent
multi-simulation MCMC/shell runs: arbitrary custom bias couplings cannot be
reliably inferred. Such independent runs can use separate mdrun invocations.
A general coupled sampler or synchronized shell solver needs a distinct protocol.
