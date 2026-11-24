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
`isCalculationOnly()` means state changes are forbidden; it does not mean that
forces are unnecessary. Nested evaluations propagate these two requirements
independently, and force-bearing loaded guests require the new command.

The 2026 GROMACS patch uses mode 1 for repeated force calls at the same step
(shell relaxation), and mode 2 for lambda MCMC candidates. Older kernels reject
these specific features with a clear message. Ordinary MD and common bias
exchange/HREX do not require the new command. Command API 14 identifies an actual
new capability, not the presence of the action fixes: maintenance backports of
those fixes retain their API numbers. Existing normal-mode calls, including
repeated full calculations at the same numerical step, keep their semantics.

The host supplies context explicitly rather than inferring it from the numerical
step: an existing caller can intentionally repeat full calculations at one step.

External integrations that reject API versions above a fixed maximum need their
accepted range updated to include 14. The maintained LAMMPS integration patch
includes this change; its existing MD commands and MPI injection are unchanged.
