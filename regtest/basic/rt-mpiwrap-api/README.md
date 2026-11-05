# Common mapped-MPI API regression

`main.cpp` uses the shared `src/mpiwrap` launcher and four ranks. The ordinary
regtest runner supplies `-n 4` from `mpiprocs=4`; it uses `mpiexec` for a native
build and the internal thread-MPI launcher for a serial build. No optional
PLUMED module is needed. Each callback owns its buffers, handles, and objects.

The same body exercises all 118 mapped functions. The tests verify values,
rank/group/topology relationships, derived datatype layout, request completion
indices/statuses, buffer attachment, packing, persistent requests, and inquiries.
Completion scenarios run twice, with explicit statuses and ignored statuses.
The wrapper uses canonical `nullptr` for ignored status arrays; native builds
use their actual `MPI_STATUSES_IGNORE` value. Data, completion indices, and
request-handle updates are checked in both modes.
Ready sends have matching receives posted before a barrier. A full native MPI
run also exercises standard reduction and error-handler callback registration.
Providers intentionally omit these callback constructors; the test verifies
their clean unavailable exceptions instead. Error-handler callbacks are never
installed on a communicator or triggered.

`thread-mpi-contract.json` is an explicit expected capability list. The test
checks **both** its available and unavailable entries, so losing a binding cannot
turn a passing semantic test into a silently skipped test. Every missing entry
is called with structurally valid dummy output storage and must throw a
`std::runtime_error` identifying that MPI function. These dummy calls never
reach a native MPI implementation.

Native collective entries are tested with derived types as well. When
`Allgather` is implemented by the shared fallback, the test verifies that a
custom contiguous datatype is rejected clearly. This limitation belongs to the
fallback, not to native MPI or other operations on derived types.

`MPI_Init` and `MPI_Finalize` belong to the shared launcher: the rank body checks
initialization/finalization state and the caller checks finalization after the
launch. Calling either lifecycle function directly through the tMPI frontend
must throw because the launcher owns it. `MPI_Abort` is run by the after-hook in
a separate invocation of the same executable; the checker requires a diagnostic
proving the intended call was reached, a failing exit status, and no return from
the call. A timeout is a failure, never a passing abort test.

`check-coverage.py` reads the generator's actual argument metadata and requires a
real call in the test source or an explicit lifecycle/destructive owner from
`lifecycle-coverage.json` for every mapping. Comments do not satisfy coverage.
It also checks that the generated unavailable-call and capability checks match
the current inventory. Run it before a mapping change is committed:

```
python3 regtest/basic/rt-mpiwrap-api/check-coverage.py
```

For a deliberate inventory change, first add a semantic scenario and review the
capability contract, then regenerate the two local includes with `--write`.
This mechanical check prevents forgotten scenarios; semantic quality still
requires review. The test does not claim exhaustive parameter/error-path
coverage of the MPI standard.

An already built executable can be checked independently, including the
isolated abort case:

```
python3 check-isolated.py ./exe
python3 check-isolated.py mpiexec -n 4 ./exe
```

Native builds require MPI-2.2 or later for `MPI_C_BOOL` and `MPI_INT64_T`, as do
these particular datatype scenarios. They do not change the MPI version
requirements of the provider itself.

## Upstream nonblocking progress limitation

The bundled upstream tMPI `Test` and `Testany` implementations do not process
newly posted incoming envelopes. Four ranks that each post `Irecv` and `Isend`
and then repeatedly test their receives do not progress. Calling `Waitall`
afterward completes the transfers. This violates the [MPI progress guarantee
for repeated receive tests](https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/node75.htm),
which does not require the sender to complete its nonblocking send first.

Consequently the tMPI capability contract excludes `Test` and `Testany`;
`Testsome` is also unavailable because its shared replacement needs a valid
`Testany`. The ordinary API test verifies their unavailable exceptions and
uses `Wait`/`Waitany` for the corresponding communication scenarios. Native MPI
still executes the complete nonblocking-test scenarios. This does not restrict
`Isend` or ready-send emulation: completing these operations through blocking
waits remains valid.

`tmpi-progress.cpp.in` is the minimal native-tMPI-only reproducer, retained for
checking a future upstream fix. It is not part of the normal executable build
and does not link PLUMED's wrapper. It polls for a bounded number of attempts,
then completes outstanding requests before leaving the rank callback. Its exit
status is **nonzero when the defect is reproduced**; it is not an expected-fail
CI test. Compile it against a standalone `src/thread_mpi` CMake build:

```
c++ -std=c++17 -DSIZEOF_LONG_LONG_INT -D__PLUMED_HAS_THREAD_MPI \
  -I../../../src/thread_mpi -x c++ -c tmpi-progress.cpp.in -o progress.o
c++ progress.o /path/to/build/libplumed_thread_mpi.a -o progress
python3 -c 'import subprocess; subprocess.run(["./progress"], check=True, timeout=30)'
```

The test's `Waitsome` scenario accepts any nonempty valid completion subset,
rather than requiring all ready requests in one call. A production replacement
must additionally preserve the standard's [fairness requirement](https://www.mpi-forum.org/docs/mpi-5.0/mpi50-report/node78.htm);
a naive single `Waitany` call does not establish that property. The tMPI
capability contract therefore also excludes `Waitsome` until a standard-safe
implementation is available; its unavailable exception is checked. Native MPI
retains the full `Waitsome` scenario.
