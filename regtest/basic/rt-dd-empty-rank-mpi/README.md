# Empty-rank domain exchange

This native two-rank test places both atoms on rank 0, then both on rank 1,
then one on each rank. Each layout is tested with synchronous and asynchronous
exchange for the first step (positions, masses and charges) and the next step
(positions only): 12 analytic cases in total. A distance restraint has energy
0.5 and atomic x forces +1 and -1; every force component must be finite and
agree within 1e-12. This checks both exchange paths against the same analytic
result without rounding the values before comparison.

In the asynchronous case with rank 0 owning all atoms, the next empty message
arrives after the receive buffer is full. A debug standard library must not
reject construction of that zero-count receive pointer. Keep both empty
messages in the protocol, including their tags and order.

The MD arrays use explicit pointer shapes supported by both the v2.9 and
newer wrappers and retain dummy storage on empty ranks. `Atoms::setAtomsNlocal`
also retains at least one element in its send buffers. The receive capacity
is the global atom count; a globally empty system cannot enter this exchange
because `setAtomsGatindex` leaves `shuffledAtoms` zero on every rank. Thus this
regression exercises zero remaining receive capacity, not an unreachable
zero-global-atoms exchange. The fix nevertheless avoids arithmetic on a null
`data()` pointer if the receive storage is empty.

Run through the native regression harness with MPI enabled, in both release
and `--enable-debug --enable-debug-glibcxx` builds. No OPES action or LAMMPS
interface is used.
