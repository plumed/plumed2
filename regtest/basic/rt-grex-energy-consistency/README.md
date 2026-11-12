Two replicas exercise GREX without an MD-engine acceptance calculation.

1. Each physical Hamiltonian is U_i(x)=(i+1)*x^2, with ENERGY biased as U^2.
   Supply own and cross energies, then check the local total energy differences
   against analytic values. Using the partner's serialized own ENERGY instead
   of the local Hamiltonian's cross energy fails this assertion.
2. Deposit a METAD hill on the exchange step and compare the bias difference
   against the same post-update stretched Gaussian evaluated at both states.
   Subtracting the cached pre-update local bias fails this assertion.

The test calls no update during the foreign calculation. It is compatible with
native MPI and the existing internal MPI test harness. A serial-only kernel
with its MPI adapter can run the unchanged rank body through a native MPI host.
