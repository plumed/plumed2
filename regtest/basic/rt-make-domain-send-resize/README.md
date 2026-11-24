# Outstanding domain send buffers

Two embedded thread-MPI ranks exchange one atom each. Rank 1 deliberately
waits before receiving rank 0's buffers. Rank 0 destroys its PLUMED instance
as soon as its own receives complete. Destruction must wait for its outstanding
sends before releasing storage. The resize companion instead grows the local
atom buffers before destroying the instance. Both verify the received positions.

AddressSanitizer catches the old use after free in both cases; ordinary builds
also check numerical values. The five-second coordination deadline prevents a
missing send from hanging the test indefinitely.

With a native-MPI kernel, config selects the corresponding MPI profiling-interface
fixture and its reference, launching two ranks. This checks numerical bias/forces
and completion of the same domain sends before destruction/resizing. Every
fixture therefore executes on either supported backend, including strict
PLUMED_ALL_TESTS coverage jobs.
