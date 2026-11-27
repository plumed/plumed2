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

Local auxiliary modes 1/2 permit independent simulations and models whose
inter-replica communication occurs only in update(), such as multiple-walker
metadynamics and MAXENT. While an action executes prepare(), calculate() or
apply(), the kernel rejects communication through its multi_sim_comm (including
PLMD::Communicator copies/splits) in these modes and names the action. This
prevents wrong conditional lambda energies and unmatched shell collectives.
The diagnostic is written directly to stderr, then the local MPI job is aborted:
throwing on replica leaders alone could strand other ranks before host error
handling. Metadata rank/size queries remain allowed. GREX's synchronized mode-0
evaluations, ordinary updates and checkpoint communication remain supported.

Actions/plugins calling raw MPI through Get_comm() bypass this check and must
not communicate across replicas in local auxiliary modes. The check does not
make a genuinely coupled sampler or asynchronous shell solver supported; those
need a coordinated evaluation algorithm. Older kernels use the documented
legacy shell fallback; lambda MCMC requires the API14 protocol. Simultaneous
expanded-ensemble sampling and replica exchange remain unsupported by the host.

For PINES with `PL_REFRESH=1`, ordinary evaluations still request all candidate
atoms and select pairs from the current coordinates. Refreshing only after
foreign coordinates would restore the reproduced ordinary pair-crossing error.
The solvated regtest/pines/rt-pines input (NaCl and 508 explicit waters, 1526 atoms)
was benchmarked with five stationary serial runs of 500 steps after 20 warmup steps.
Compared with pre-fix `832f9c9a1`, the mean atom request rose from 780.5 to 1526 (1.955x)
and the median PINES evaluation rose from 0.277 to 0.410 ms (1.478x). These are action
evaluation/request measurements, not total MD runtime or network timings.
Avoiding full requests safely needs a displacement-bound/candidate-validity
scheme; reusing an arbitrary previous-frame list is insufficient. Local lambda
trials can reuse initialized membership without this foreign-coordinate refresh.
