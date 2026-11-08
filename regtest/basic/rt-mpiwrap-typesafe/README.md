# MPI pointer type checks

This serial test does not initialize or select an MPI provider. It checks local
integer/pointer/wrapper communicator types, const and untyped pointers, rejected
wrong-type arguments, provider pointers, exception translation, and
metadata forwarding through ordinary, CLTool and GREX commands. Raw legacy
communicator flags and unknown future communicator annotations remain accepted.
The native integer in a legacy provider-mode argument is deliberately smaller
than wrapper storage.

Both sides use only extra metadata in the existing safeptr flags. Old kernels
ignore the communicator exclusion bit and provider type code 6. New kernels
accept legacy void/unannotated communicator pointers. The new provider API
requires dedicated type code 6 when metadata is present, including for opaque
pointers; it does not need another exclusion bit. Native-MPI kernels ignore
injection entirely, including null, invalid and incorrectly typed descriptors.
C99's untyped interface is unchanged; C11 and C++ callers can carry the new
annotations. Communicators remain
pointer-only; values are rejected if the local type distinguishes them.
Run `bash regtest/basic/rt-mpiwrap-typesafe/check-callers.sh` from the repository
root for C11/C++98/11/17 transport probes and native kernel exception checks without an MPI installation.
The opaque-pointer probe deliberately leaves its implementation type incomplete.
The native probes execute the actual `TypesafePtr` checks for integer, opaque-pointer
and void-pointer handles, including nested-pointer rejection and legacy flags.
The probes also check type-only wrapper inclusion and verify that TypesafePtr.h
and Communicator.h do not pull in provider dispatch.
To annotate communicators,
the local MPI declaration must precede Plumed.h (native mpi.h, the wrapper
mpi.h or mpi-types.h, or GROMACS MPI bindings). Host-adapter mode includes native mpi.h itself.

`__PLUMED_WRAPPER_C_MPI_COMM_CHECK` defaults to 1 in the existing C11 path.
`__PLUMED_WRAPPER_CXX_MPI_COMM_CHECK` defaults to the anonymous-namespace
switch: shared public C++ wrappers are unchecked by default, while anonymous
wrappers classify arguments when MPI declarations are visible. Explicit shared
C++ checking requires matching MPI declarations and options in every translation
unit; missing declarations are diagnosed. Disabling these checks preserves the
provider type tag and all generic type checks.

Run `check-header-odr.sh` for mixed MPI-aware/unaware translation units at O0,
with reversed link order and O2/LTO when supported. `check-callers.sh` exercises
explicit enable/disable, anonymous defaults and shared defaults independently.
The API directory's `check-host-header-linkage.sh` accepts MPICC/MPICXX and checks
C99/C++98+ host adapters across both namespace and linkage modes.
