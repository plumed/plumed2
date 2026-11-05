# MPI provider wrapper

This module lets one PLUMED kernel use MPI through one selected runtime while
preserving the existing Communicator source interface. It has three mutually
exclusive configure modes: native MPI, internal thread-MPI, and fully serial.
The native build keeps PLUMED's existing direct MPI path. The thread-MPI build
uses the pinned GROMACS thread-MPI implementation. A fully serial build has
neither provider and retains the existing serial behavior.

## Configuration and capability reporting

Native MPI is discovered through --enable-mpi (enabled by default). Internal
thread-MPI defaults on when native MPI is unavailable; use
--disable-thread-mpi to turn it off. Explicitly requesting both available
providers is an error. --disable-mpi --disable-thread-mpi configures a fully
serial build.

__PLUMED_HAS_MPI means native MPI only. __PLUMED_HAS_THREAD_MPI means the
internal runtime is built. __PLUMED_MPIWRAP is always defined by configure:
availability of the wrapper frontend is separate from selecting a provider.
Native MPI still uses its actual handle types and direct calls. Config reports
has mpi when either MPI provider exists, and also reports real_mpi and
thread_mpi separately. The CLI provides --has-mpi for combined capability,
--has-real-mpi for a native build or a selected real-MPI provider, and
--has-thread-mpi for available internal thread-MPI (or a selected thread
provider). Both capabilities may be true. --mpi-ranks-are-threads describes
only the selected provider, independently of initialization; CHDIR safety uses
the separate initialized-provider runtime check. Queries do not initialize MPI.
The normal regression suite runs with direct MPI when compiled with MPI, or
with internal thread-MPI when compiled without MPI. Real-MPI initialization
and finalization belong to the executable or host application, never the
serial-built PLUMED kernel. Host injection is tested by a separate MPI program.
The regression runner preserves old test settings and native launcher overrides;
older executables without thread-MPI support retain their native MPI launcher.

Provider-mode Communicator construction does not evaluate MPI_COMM_SELF or
select a provider. Its default handle is explicitly unresolved until an MPI
operation; then its lazy accessor selects the provider and obtains SELF once
MPI is initialized. Inject before the first provider-dependent constant or
call. Provider selection and initialization are separate; serial operations
still use local fallbacks until MPI is initialized. A fully disabled build
retains its inert serial surrogate. Internal capability inquiries do not select a provider. Explicit CLI inquiries
(`--has-mpi`, `--has-real-mpi`, `--has-thread-mpi` and
`--mpi-ranks-are-threads`) do select the provider so they can report an
environment-selected implementation; they never initialize MPI.
The separate `mpiwrap/mpi.h` inquiry reports whether an installed, initialized
provider executes ranks as threads; the PLUMED action uses it for CHDIR safety.

## Handle and ABI boundary

In provider mode, each public MPI handle is a trivial, standard-layout POD
record. Its zero-filled storage starts at offset zero and holds the native
handle representation. This is a deliberate compatibility contract:
existing PLUMED callers can keep passing &Communicator::Get_comm() to nested
kernels without changing their source or creating export buffers. A compatible
receiver copies exactly the selected provider's reported native handle size
from that address into zero-filled wrapper storage. The storage capacity
is validated before use. Matching bytes do not make a handle valid across
unrelated MPI runtimes.

The wrapper and native adapter communicate through a versioned C-compatible
provider descriptor. Injection uses the existing command interface, advertised
by getApiVersion >= 11; the kernel symbol table remains version 4. The provider
descriptor has its own ABI version. Fixed entry layouts are frozen within a provider ABI major; new optional
capabilities use named metadata entries. Each injected implementation takes
one pointer per argument plus a result-storage pointer, following the draft;
there is no array of argument pointers. Function pointers use a generic
function-pointer carrier and are cast back to the exact generated per-function
transport type before invocation. IDs are hashes of function names used for
validation, never addresses or table offsets. Entries append contiguously in
array order, including conditional entries. PLUMED handles native storage,
array strides and status conversion; injected implementations call native MPI.
The nullable description is copied at installation for diagnostics only.
`ranks_are_threads` describes execution; absence of a function means unavailable.

PLUMED_MPIW_HANDLE_CAPACITY defaults to max(8, sizeof(void *));
PLUMED_MPIW_STATUS_CAPACITY defaults to 128 bytes. Rare overrides can be supplied
as compiler definitions, for example
-DPLUMED_MPIW_HANDLE_CAPACITY=24 -DPLUMED_MPIW_STATUS_CAPACITY=192.
Use the same definitions for the kernel and every client using wrapper types,
and rebuild them together. Installed headers retain the defaults; configure
does not select or propagate these overrides. Native MPI's direct
types are independent of these capacities, as is the thin native provider;
provider installation checks that its native values fit the receiving kernel. Handle storage is a zero-filled
byte-array POD. Provider ABI major 3 uses fixed per-function pointer transport
and rejects older incompatible provider layouts before reading their entries.
Providers must fit those configured representations,
and installation validates the capacities and native scalar metadata before
publishing function pointers. The frontend signatures
follow MPI-5 const correctness. The native adapter keeps any compatibility
cast for older native headers local to an input-only MPI call; callback types
and writable arguments retain their exact types.

Internal thread-MPI launching uses the ordinary command interface: configure
CLTool arguments, including `--tmpi-threads N`, and execute `CLTool run`.
The initial PLUMED object dispatches the launch; each rank constructs its own
tool and communicator state. No public kernel-handle type or separate launcher
function is needed. The kernel owns internal worker startup and finalization.
Nested/concurrent launches are rejected. Sequential reuse with different rank
counts must be covered by the maintained launcher check. An unrecoverable rank
exception terminates the process rather than unwinding through blocked peers;
recoverable errors use the existing command exception translation.

Standalone C++ programs can use the internal `mpiwrap_provider::launchRanks`
entry declared in `thread-mpi-launcher.h`. CLTools uses this same entry for
internal launch. Select `internal_thread_mpi` for embedded workers, or
`standalone_owned` for a native-MPI program launched by `mpiexec`; the latter
initializes MPI, checks the requested process count, executes one callback per
process, then finalizes. `host_owned` requires already initialized native MPI
and never initializes or finalizes it. The supplied native communicator bytes
are valid during the callback; construct rank-local objects there. This is an
internal C++ facility, not a new `Plumed.h` or symbol-table API.

C++ make regtests opt in with `plumed_mpi_harness=yes` and a positive
`mpiprocs` in `config`. The runner passes `-n N` to the executable and uses
`mpiexec -np N ./exe -n N` for direct MPI or `./exe -n N` for internal tMPI.
The executable parses `-n` and calls the shared launcher. Existing make tests
keep their previous behavior. `regtest/basic/rt-mpiwrap-launcher` is a four-rank
launcher smoke test, not the complete mapped-API test (which is still pending).

## Inject the host MPI from Plumed.h

When the host needs to inject the MPI implementation already used by its MD
code, define `PLUMED_MPIW_HOST_ADAPTER` in exactly one translation unit before
including `Plumed.h`. This embeds the host adapter independently of the wrapper
implementation, so the host can link normally to `libplumed`. Define
`__PLUMED_WRAPPER_IMPLEMENTATION` as well only when embedding the wrapper
implementation, for example for runtime loading. Include the host MPI bindings first. For GROMACS
thread-MPI, include `thread_mpi/mpi_bindings.h` first and ensure
`GMX_THREAD_MPI` has its numeric build value. Include `gromacs/version.h` so
`GMX_VERSION` is available, or explicitly define `PLUMED_MPIW_GROMACS_VERSION`
using the same numeric encoding (2026.0 is 20260000). Unknown versions and
versions earlier than 2026.0 are rejected when building the injected implementation. This opt-in header embeds the
same generated adapter body, so it does not need a separate provider source or
provider library.

The host initializes real MPI, then creates an ordinary object and injects
before any command that uses MPI:

```cpp
plumed p = plumed_create();
// Or: plumed p = plumed_create_dlopen(path);
plumed_cmd(p, "injectMPIImplementation", PLMD::plumed_mpi_implementation());
```

C++ callers use ordinary `Plumed` constructors and the same `cmd` command.
The command is available with kernel API version 11 (`getApiVersion`);
the symbol table remains version 4. Its descriptor argument is the descriptor
address itself, not the address of a pointer. Errors use the existing command
exception/nothrow translation machinery. Object creation retains its historical
invalid-object failure convention.

Injection is **library-wide**, applying to every object in that loaded kernel.
Inject into each newly created object: destroying the last object may unload
the kernel, so a later object is not guaranteed to retain an earlier injection.
Compatible repeated injections are allowed; replacing a selected provider is
not. Independent copies of the generated header may have different local
adapter and constant-accessor addresses. Repeated injection accepts those
copies only when native runtime identity, ABI/layout, function availability and
signatures, and constant widths/metadata agree. The original callback table
remains installed. Objects own their normal loader references independently. Applications
that need explicit lifetime control can retain their own `dlopen` handle and
use `plumed_create_dlsym`.

Constructing an object does not select a provider. Its default communicators
remain unresolved until actual MPI use and obtain SELF only after MPI is
initialized. `Get_comm()` can still be used to forward an unresolved
communicator before initialization. In an internal thread-MPI build, its
storage contains the provider's `MPI_COMM_NULL` sentinel, and wrapper MPI calls
reject that value. Other providers retain their native pre-initialization
behavior. The unresolved state is independent of native handle bytes, so a
native zero-valued communicator remains distinguishable from a placeholder.

The host must keep the injected adapter code and MPI runtime loaded while the
kernel can call them, finalize every PLUMED object before finalizing MPI, and
own all real `MPI_Init`/`MPI_Finalize` calls. Copying a descriptor does not keep
its executable code loaded. PLUMED does not initialize real MPI merely because
a provider was injected or loaded.

Existing `setMPIComm` callers continue passing the address of their native
communicator. PLUMED imports the provider-sized bytes into zero-filled wrapper
storage, whose native byte prefix has the same layout. This byte-copy contract
is the reason for the chosen layout: callers, including nested PLUMED actions,
need no export buffer or communicator-forwarding changes.
The generated adapter omits thread-MPI functions and constants unavailable in
the selected GROMACS headers; PLUMED reports calls to omitted functions as a
clear runtime error.

## Build a native provider

Use `plumed mkmpi --compiler /absolute/path/to/mpicc --out provider.so`
(or an MPI C++ compiler driver) to build the shared native adapter. The command
reuses configure-tested compiler/shared-library flags, relies on the `.cpp`
suffix for language selection, and adds no language-standard flag. The adapter
is C++98-compatible and does not link a PLUMED kernel.
Use a `.dylib` output name on macOS. Pass each additional compiler or linker
argument separately with repeated `--cflag ARG` or `--ldflag ARG` options;
arguments are never shell-evaluated. The selected compiler must support C++.
The command compiles a small generated stub against the installed public
`Plumed.h`; it does not need adapter source files or the original source
checkout. Failed builds preserve an existing output.

Building a provider does not initialize MPI: the host remains responsible
for its MPI lifecycle. Set `PLUMED_MPI_WRAPPER` to an absolute provider-library path
or an absolute executable MPI compiler path. Explicit pre-creation injection
takes precedence; otherwise the environment selection precedes the internal
provider. A loadable library must expose the provider factory. A non-loadable
executable compiler invokes the same `mkmpi` command in a private temporary
directory, with compiler diagnostics inherited by the process. Files with
shared-library suffixes are never executed as compiler commands. The kernel
retains the loaded provider library; temporary build files are removed after
loading. This selection does not initialize the host MPI runtime.

For troubleshooting, set `PLUMED_MPI_WRAPPER_FORCE=1` together with
`PLUMED_MPI_WRAPPER` before any provider is installed. This makes that
environment-selected provider take precedence over an incoming injection.
Without the force flag, injection retains priority. An absent provider path,
an invalid force value, or a late attempt to replace an installed provider
raises an error. Injected callback code remains host-owned, including when injection is bypassed
by the force setting. Already published callbacks are never replaced.

## Generated interface and runtime files

All module files live directly in src/mpiwrap. The machine-readable inputs
are support-matrix.tsv, argument-metadata.json, constant-metadata.json,
type-metadata.json, signature-versions.json, availability-metadata.json, and
thread-mpi-api.tsv. Signature revisions are independent of MPI standard versions.
Build the repository formatter once with make -C astyle. generate.py validates
these inputs and generates the wrapper declarations, provider ABI, constant
accessors, route tables, signature checks, and coverage reports in the same
directory. The checked-in generated files are reproducible
with:

    python3 src/mpiwrap/generate.py --check

The .cpp.in and .c.in files are compile-only tests and probes; the module's
wildcard source build does not compile them. The implementation unity source
mpiwrap.cpp supplies the registry and frontend in wrapper mode, including when
the embedded runtime is disabled. Its internal adapter is compiled only when
thread-MPI is configured. mpiwrap-launcher.cpp
includes the internal rank launcher. mpiwrap_implementation.cpp is an opt-in
translation unit for a provider DSO; it is not compiled into ordinary PLUMED
kernels. Imported thread-MPI sources remain a separate module at
src/thread_mpi.

## Current implementation and limits

The internal provider has a tested launch path and a usable initial MPI subset,
including communicator rank/size, collectives, point-to-point operations,
status/count handling, and selected handle/datatype operations. Functions with
no safe route remain declared and their availability inquiry reports false;
calling them reports the missing capability. The standard variadic
MPI_Pcontrol remains direct in native mode and unavailable through provider
forwarding. Custom tMPI error handlers remain excluded as specified by the
stable provider surface.

The native adapter can publish 113 of the 118 retained function entries;
actual publication also depends on the host MPI version and supported constants. The omitted entries are MPI_Init/Finalize
(lifecycle-owned), MPI_Pcontrol (variadic), MPI_Comm_create_errhandler and
MPI_Op_create (callback bridges).
Their frontend declarations remain available and their inquiries report the
missing provider capability.
Custom error-handler creation stays outside the provider scope: PLUMED leaves
the host application’s error-handler policy unchanged. Custom reduction support
is a required first-release follow-up, starting from `mpiw_old_draft`; the
current omission is not the intended final release contract. Review and document
its callback ABI assumptions and test operation lifetime before freezing the
injection calls. Direct native MPI calls remain available. The draft's callback-reference conversions are
useful precedent, not a substitute for that review. No callback registration
or retained-handler allocation occurs for these unavailable APIs. The mock-native probes compile against real
provider headers and exercise call-through without initializing MPI; they are
not real MPI runtime communication tests. Native MPI startup and multi-process
communication were blocked in the development sandbox by process-manager and
socket permissions, so those semantics still need verification outside that
environment. The pinned thread-MPI source has a one-rank
out-of-place Allreduce defect; the wrapper handles its size-one identity case
without relying on that defective operation.

The embedded and injected tMPI adapters use the same thin native-call renderer.
Shared PLUMED code supplies Allgather/Allgatherv using the provider's basic C
predefined datatypes and standard C sizeof values. Custom and compound types
throw clearly on these fallback paths; provider-native collectives retain their
normal datatype support. The tMPI metadata explicitly states that its
communicators are intracommunicators, since it lacks Comm_test_inter.

Ready sends Rsend/Irsend delegate to Send/Isend; callers must have posted the
matching receive before starting either ready send. Ssend/Issend remain
unavailable because standard sends can complete through eager buffering before
a receive starts. Bsend remains unavailable because it needs the separate
attached-buffer completion contract.

Wait uses the shared Waitany replacement because upstream Wait leaves a stale
request handle. Test and Testany remain unavailable: a bounded standalone
reproducer demonstrates that they do not progress matching nonblocking messages.
Testsome therefore lacks a compliant prerequisite, and Waitsome remains
unavailable because a Waitany-only replacement cannot guarantee its required
fairness. Shared some-function replacements remain available to providers that
supply the necessary compliant primitives. Testall remains unavailable because
upstream modifies requests while returning a false completion flag.
Errhandler_free is omitted because upstream can free a borrowed runtime handler;
no replacement silently assumes ownership.

## Verification

The focused checks can be run from the repository root:

    python3 src/mpiwrap/validate-manifest.py
    python3 src/mpiwrap/generate.py --check
    sh src/mpiwrap/check-provider-registry.sh
    bash src/mpiwrap/check-kernel-lifecycle.sh
    bash src/mpiwrap/check-regtest-runner.sh
    MPICH_CC_WRAPPER=/path/to/mpich/bin/mpicc \
      MPICH_CXX_WRAPPER=/path/to/mpich/bin/mpicxx \
      OMPI_CC_WRAPPER=/path/to/openmpi/bin/mpicc \
      OMPI_CXX_WRAPPER=/path/to/openmpi/bin/mpicxx \
      bash src/mpiwrap/check-native-adapter.sh

Header and generated signature checks for local MPICH, Open MPI, and pinned
tMPI use:

    MPICH_CXX=/usr/bin/clang++ OMPI_CXX=/usr/bin/clang++ \
      sh src/mpiwrap/check-generated.sh \
      /path/to/mpich/bin/mpicxx \
      /path/to/openmpi/bin/mpicxx \
      src/thread_mpi

Internal launch semantics are exercised by
src/mpiwrap/check-thread-mpi-cli.sh src/lib/plumed-static after a configured
thread-MPI build. The tests cover 1, 2, and 4 ranks, callback error handling,
provider lifecycle, raw communicator byte-prefix handling, native/legacy
runner selection, and CHDIR rejection while thread-MPI is active.

## Maintaining GROMACS compatibility

Injected GROMACS thread-MPI support starts at GROMACS 2026. The embedded
runtime remains pinned to stable 2026.4 in the separately licensed thread_mpi
module. Earlier versions lacked Type_free, retaining derived datatype storage
until finalization; 2025 also corrected Error_string's length-output type.

Audit historical headers from a local full GROMACS clone without network access:

    python3 src/mpiwrap/check-gromacs-history.py /path/to/gromacs --check

To review a future release, rerun with --through vYYYY.N without --check,
inspect the resulting gromacs-history.json diff, and update availability rules
only after reviewing the corresponding implementations. The inventory records
47 stable releases from 2020 through 2026.4; declarations alone do not establish
standard-compatible semantics. Updating the embedded importer does not replace
this review.

The audited development-main snapshot 73d6be5d778e adds UINT64_T and makes Wait
clear completed requests, relative to 2026.4. It introduces no further public
functions or signature changes. This is a recorded snapshot, not a promise
about a future 2027 release. Keep the 2026.4 request-clearing workaround and
omit newer constants from hosts that do not provide them.

Native MPI availability uses MPI_VERSION/MPI_SUBVERSION rather than testing
whether a function name is a macro. The MPI-2 functions Comm_f2c, Finalized,
Comm_get_errhandler and Comm_set_errhandler are omitted for older implementations;
MPI_IN_PLACE also requires MPI-2. Optional long-long predefined types are
published conservatively from MPI-2.1, when required by the standard. The same
guards protect implementations and table rows. Synthetic MPI-1.1 probes verify
that older headers need not declare omitted functions or constants.

The capacity-override installation check uses a configured checkout:

    ./configure --disable-mpi
    MPICH_CXX=/usr/bin/clang++ OMPI_CXX=/usr/bin/clang++ \
      bash src/mpiwrap/check-capacity-install.sh /path/to/mpich/mpicxx \
      /path/to/openmpi/mpicxx

It runs the public install-headers recipe and compiles installed headers
with explicit compiler definitions for 24-byte handles and 192-byte status
storage. The build/runtime matrix above remains separate validation.

A separate MPI host checks injection into an already-built serial kernel:

```sh
MPICH_CXX=/usr/bin/clang++ bash regtest/mpiwrap/check-mpi-host-injection.sh \
  /path/to/mpicxx /path/to/mpiexec /absolute/path/to/libplumedKernel.dylib
```

Use matching compiler/launcher paths from one MPI installation. The harness
builds only the host, initializes and finalizes MPI there, and checks the
`rt-mpi0` outputs with 1, 2 and 4 ranks using both default and explicit kernel
selection. It unsets provider-selection variables so successful communication
proves injection, rather than environment loading. Repeat with another MPI
installation against the same serial kernel. Existing C/Fortran/Python command
interfaces can launch internal thread-MPI without a new public handle type;
Python calls should release the GIL during the long-running command.

### MPI installed after a serial kernel build

The dedicated Linux injection job installs Open MPI after building PLUMED
serially. A small C++98 host reruns the suite through injection; focused checks
also exercise 1, 2 and 4 ranks and both creation paths. Debug and Intel jobs add
a posterior Open MPI injection pass after their normal suites. The focused
checks also build a provider DSO with `plumed mkmpi` and pass `--no-injection`
to the same host, selecting that DSO through `PLUMED_MPI_WRAPPER`. The same
checks can be repeated manually with another MPI implementation.
The same serial kernel is tested again through injection;
there is no PLUMED rebuild. The host owns MPI initialization
and finalization. `PLUMED_CLI_PROGRAM_NAME` selects it for CLTool tests only,
and the regression runner derives the installed library paths. External-engine
tests retain their own initialization and capability
checks, so MPI cases requiring a natively built external engine remain skipped.
No GROMACS build is required. The standalone embedded tMPI tests also run.

Run the same host test locally with an already built serial kernel:

```sh
bash regtest/mpiwrap/check-mpi-host-injection.sh /absolute/path/to/mpicxx \
  /absolute/path/to/mpiexec /absolute/path/to/libplumedKernel.so \
  /absolute/path/to/provider.so
```

The optional fourth argument enables environment-provider cases in addition to
injection cases. Omitting it retains the injection-only test. The host consumes
`--no-injection` as its first option and does not forward it to PLUMED. It owns
MPI initialization and finalization for parallel execution in both
modes; `--no-mpi` queries avoid MPI lifecycle calls.


## Integration with existing PLUMED code

The MPI provider implementation keeps existing callers and native-MPI initialization
unchanged. In particular `src/main/main.cpp`, kernel initializer/symbol-table
code and existing wrapper entry-point implementations are unchanged. There is
no new kernel symbol-table entry or public object type.

| Files | Reason for the feature change |
| --- | --- |
| `configure.ac`, generated `configure` | Discover native MPI as before; default to internal tMPI otherwise, with an explicit opt-out. Define the wrapper infrastructure independently of runtime selection. No capacity configure options or module `.on`/`.off` manipulation. |
| `CMakeLists.txt` | Enable the wrapper and optional embedded Win32/POSIX thread backend; retain registration objects with each platform's linker syntax. Existing Windows configuration remains supported. |
| `src/.gitignore` | Track the imported flat runtime sources and provenance record using the existing module-ignore conventions. |
| `src/{core,generic,tools}/Makefile` | Declare the new module/header dependencies for normal builds and installed headers. |
| `src/config/Makefile` | Export native-provider compiler/link flags for `mkmpi`; regenerate options when this Makefile changes. The provider uses a C ABI and does not link a C++ PLUMED runtime. |
| `src/maketools/update-config-txt.sh` | Separate configured real-MPI and internal-tMPI capabilities while retaining combined `has mpi`. |
| `src/core/CLToolMain.cpp` | Parse internal rank counts, invoke the shared launcher, construct rank-local commands and report selected-provider runtime capabilities. |
| `src/core/PlumedMain.cpp` | Expose injection through the existing command mechanism and advertise command API version 11. Construction still permits injection before provider-dependent use. |
| `src/tools/Communicator.{h,cpp}` | Replace compile-time serial branches with runtime initialization checks, preserve serial fallbacks, import native bytes through the single boundary helper and resolve default SELF before communication. An explicit flag distinguishes deferred SELF from a valid all-zero native handle. Private buffer descriptions use const references for fatal cppcheck diagnostics; receive buffers remain writable. |
| `src/tools/PlumedHandle.cpp` | Forward the selected implementation to a compatible nested kernel through the existing command API; preserve older-kernel behavior and release the created reference on failure. |
| `src/generic/Plumed.cpp` | Reject process-global CHDIR while a thread-MPI runtime is initialized. Communicator forwarding at existing call sites is unchanged. |
| `src/core/Value.h` | Suppress only the observed cppcheck false positive on a member-reference accessor, without changing its code or return type. |
| `src/maketools/codecheck` | Exclude the imported runtime like other vendor code, avoid opt-in adapter configurations in kernel analysis, and exempt only structural rules incompatible with generated C ABI headers. Semantic diagnostics remain enabled. |
| `src/wrapper/Plumed.h` | Add the opt-in generated host/native adapter and its one implementation accessor. Existing C/C++ wrapper implementations and exception handling are unchanged. |
| `scripts/mkmpi.sh` | Build a standalone provider from the already-installed public header, using the chosen MPI driver, C++98-compatible bodies and atomic output publication. No extra implementation sources are installed. |
| `plugins/pycv/FindPlumed.cmake` | Link external MPI only for a native-MPI build; embedded tMPI does not require external MPI headers or libraries. |
| `regtest/scripts/run`, `mpi*.cpp`, `mpi-test.h` | Select direct MPI, an injecting CLI host or internal ranks while keeping old `mpiprocs` values and old-kernel behavior. External MPI lifecycle remains owned by the executable. |
| Converted existing MPI regtests | Mechanical lifecycle/communicator changes only; original calculations, rank counts and reference files are retained. Four cases with process-global environment or independent serial processes remain native-only. |
| New API/launcher/provider regtests | Check the mapped API, unavailable calls, transport, lifecycle, nested kernels without compiling GROMACS. |
| `.ci/check-warnings`, `.codecov.yml` | Treat the imported thread-MPI runtime like other vendored libraries for warning ownership and coverage metrics; the PLUMED MPI adapter remains checked. |
| Linux/Windows workflows | Keep master's existing jobs, including sanitizer/no-OpenMP/no-dlopen checks. Add one serial-build/full-injection job, posterior OpenMPI passes for debug and Intel, focused provider tests and explicit Windows smoke checks. |
