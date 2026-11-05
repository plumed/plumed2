# MPI wrapper integration checks

These checks exercise MPI wrapper and host-injection behavior that needs a
configured kernel, a provider, an MPI implementation, or synthetic ABI inputs.
They are run explicitly; this directory is not an automatically discovered
regtest suite. Shared MPI helper sources used by hosts remain in
`regtest/scripts/`.

`check-mpi-host-injection.sh MPI_CXX MPIEXEC SERIAL_KERNEL [PROVIDER_LIBRARY]`
is the focused injection check used by CI. It supports Open MPI and MPICH,
testing MPI host injection through both PLUMED creation paths and
an optional provider DSO. It also checks that `--no-mpi` permits injection
and configuration queries without initialization, finalization, or abort calls,
including failure paths and leading-option precedence. The remaining scripts
are focused manual checks:

```sh
bash regtest/mpiwrap/check-mpi-host-injection.sh MPI_CXX MPIEXEC SERIAL_KERNEL [PROVIDER_LIBRARY]
bash regtest/mpiwrap/check-mkmpi.sh
bash regtest/mpiwrap/check-mpi-provider-capacity.sh MPI_CXX
PLUMED_MPI_WRAPPER=/absolute/provider.so PLUMED_KERNEL=/absolute/kernel.so \
  bash regtest/mpiwrap/check-mpi-provider-chdir.sh MPI_CLI_HELPER
bash regtest/mpiwrap/check-mpi-provider-force.sh KERNEL PROVIDER
bash regtest/mpiwrap/check-mpi-provider-selection.sh PROGRAM PROVIDER MPI_COMPILER
bash regtest/mpiwrap/check-mpi-wide-native-handle.sh
bash regtest/mpiwrap/check-mpi-nested-kernel.sh RUNTIME_PROGRAM KERNEL
bash regtest/mpiwrap/check-thread-mpi-benchmark.sh PLUMED_PROGRAM
bash regtest/mpiwrap/check-wrapper-header-compatibility.sh
```

`check-mkmpi.sh` verifies provider build argument handling and publication;
capacity and wide-handle checks exercise native handle storage; provider
selection and force checks cover provider loading policy; the CHDIR and nested
kernel checks cover nested kernel behavior; the benchmark check compares
thread-MPI results with serial results; and the header check verifies standalone
wrapper language compatibility. These checks are independent of the standard
regtest fixtures and require only the inputs listed in their usage messages.
