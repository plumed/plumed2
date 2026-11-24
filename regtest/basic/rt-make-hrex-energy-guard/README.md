Exercise the legacy GREX foreign-energy command without MPI. Ordinary ENERGY
biasing and coordinate-only HREX remain supported. Defining ENERGY in HREX,
even without a bias, is conservatively rejected at cacheLocalUSwap, before
foreign-state communication. No new MD-engine command is required.

The same test source is mergeable to master: its reference expects support
there, where GREX substitutes the correct foreign physical ENERGY. The actual
cross-Hamiltonian calculation is covered by rt-grex-energy-consistency.
