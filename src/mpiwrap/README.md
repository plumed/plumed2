# MPI provider wrapper

This module lets one PLUMED kernel use MPI through an explicitly selected
runtime while preserving the existing Communicator source interface. Native
MPI builds keep PLUMED's existing direct MPI path. Serial-built kernels can
load an external MPI provider at runtime and optionally include embedded
thread-MPI. The embedded runtime uses the pinned GROMACS implementation and
starts only when the explicit internal launcher is used.

This separation also enables smaller package build matrices: one `plumed`
package can serve MPI-independent users and users of either external MPI
implementation, with small optional `plumed-mpi` provider packages. See
[Packaging and conda extensions](#packaging-and-conda-extensions).

## Configuration and capability reporting

Native MPI is discovered through --enable-mpi (enabled by default). Internal
thread-MPI defaults on when native MPI is unavailable; use
--disable-thread-mpi to turn it off. Explicitly requesting both available
providers is an error. Building embedded thread-MPI does not make it an
implicit serial fallback: launch it explicitly with --tmpi-threads N.
--disable-mpi --disable-thread-mpi builds a kernel with neither MPI runtime
linked in; injection and runtime provider loading remain available.

__PLUMED_HAS_MPI means native MPI only. __PLUMED_HAS_THREAD_MPI means the
internal runtime is built. __PLUMED_MPIWRAP is always defined by configure:
availability of the wrapper frontend is separate from selecting a provider.
Native MPI still uses its actual handle types and direct calls.
`plumed config has mpi` reports whether native or internal thread-MPI was enabled
at build time; `real_mpi` and `thread_mpi` report those build choices separately. Loading
a provider does not change this configuration record. The CLI provides
--has-mpi for combined runtime capability,
--has-real-mpi for a native build or a selected real-MPI provider, and
--has-thread-mpi for available internal thread-MPI (or a selected thread
provider). Both capabilities may be true. --mpi-ranks-are-threads describes
only the selected provider, independently of initialization; CHDIR safety uses
the separate initialized-provider runtime check. Capability queries do not
select embedded tMPI or initialize MPI. An explicitly configured environment
provider may be selected for a query. Real-MPI initialization
and finalization belong to the executable or host application, never the
serial-built PLUMED kernel. Host injection is tested by a separate MPI program.
The regression runner preserves old test settings and native launcher overrides;
older executables without thread-MPI support retain their native MPI launcher.

Provider-mode Communicator construction does not evaluate MPI_COMM_SELF or
select a provider. Its default handle is explicitly unresolved until an MPI
operation; SELF is resolved after the selected provider is initialized.
Provider selection and initialization are separate; serial operations still
use local fallbacks until MPI is initialized, including in builds without
embedded thread-MPI. Importing a native communicator without an
injected or environment-selected provider first looks beside the loaded
kernel for `libplumed-mpi-provider.so` (`libplumed-mpi-provider.dylib` on macOS).
If absent, it triggers `mkmpi` without a `--compiler` option to build the
native adapter. It searches `$MPICC`, `$MPICXX`, `mpicc`, then `mpicxx`.
This search must find a compiler compatible with the host MPI runtime.
The kernel never initializes or finalizes real MPI. The
`--mpi-ranks-are-threads` query succeeds only when a thread-based provider has
been selected.
The package `mpiwrap/mpiwrap.h` inquiry reports whether an installed, initialized
provider executes ranks as threads; the PLUMED action uses it for CHDIR safety.

## Namespaces

All C++ implementation types, functions and state live under `PLMD::mpiwrap`,
including `mpiwrap::provider_installed()`, the registry and
`mpiwrap::pointers`. Native MPI-name macros forward to this namespace. Standard
MPI type aliases remain available for unchanged callers. The C-compatible
adapter and descriptor use conditional namespace guards: C++ definitions are
in the module namespace, while C99 compilation and exported C symbol names
are preserved. The public `plumed_mpi_implementation` spelling remains available
through a using-declaration. Imported runtime definitions and its test helpers
are under `PLMD::thread_mpi`; executable entry points remain global.

Declaration-only/inline headers declare `PLUMED_HEADER_ONLY`, and shared
C headers declare `PLUMED_C_COMPATIBLE`, using the generic plumedcheck rules.
Neither property relaxes namespace or guard checks. No mpiwrap-specific
plumedcheck suppression is used.

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

Standalone C++ programs can use the internal `mpiwrap::launchRanks`
entry declared in `mpiwrap.h`. CLTools uses this same entry for
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
launcher smoke test; `regtest/basic/rt-mpiwrap-api` provides the mapped-API
checks.

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
plumed_cmd(p, "injectMPIImplementation", PLMD::mpiwrap::plumed_mpi_implementation());
```

C++ callers use ordinary `Plumed` constructors and the same `cmd` command.

The accessor returns a typed provider descriptor pointer. Safeptr type code 6
identifies it; older kernels ignore this reserved code. Provider-mode kernels
require this dedicated type whenever type metadata is supplied, including for
otherwise opaque object pointers. There is no legacy typed provider argument to
preserve. The untyped C interface remains unchecked, as for other commands.
The descriptor ABI and kernel symbol table are unchanged. Native-MPI kernels
ignore injection without reading or validating the supplied descriptor. Their
provider registry, dispatch and fallback implementations are compiled out;
standard MPI calls use the linked MPI library directly.

When the local MPI declaration precedes `Plumed.h`, communicator arguments are
checked against that translation unit's `MPI_Comm*`, preserving const pointers.
The reserved flag `0x20000000` means the argument is incompatible. New kernels
check it through `TypesafePtr::getMPIComm()` at the communicator boundary,
before discovery or byte copying. In native mode this explicit request is
necessary because MPI_Comm can alias an ordinary type such as int: normal
integer commands must ignore the communicator exclusion flag. Provider-mode
`get<const MPI_Comm*>()` can also check it because the wrapper class is distinct.
Old callers omit the flag and old kernels ignore it. In provider mode the
native size is not compared against `sizeof` the wrapper class.

Communicator commands retain the pointer-only interface: pass `&comm`.
An integer or distinct opaque-pointer handle passed by value is rejected when
the local type is known. Explicitly erased void pointers remain unchecked,
including when MPI_Comm is itself void*. C99's untyped interface is unchanged.
Host-adapter mode includes native `mpi.h` before setting up its type checks.
Provider checks use type code 6 alone. No additional provider exclusion bit is
needed; `0x40000000` and higher bits remain reserved.

`regtest/basic/rt-mpiwrap-typesafe` checks exceptions, forwarding and legacy
metadata. Its `check-callers.sh` checks C11 and C++98/11/17 transport against
synthetic integer, opaque-pointer and void-pointer MPI types, plus the actual
wrapper class, without MPI initialization.
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

Use `plumed mkmpi --out provider.so` to search for an executable compiler in
the order `$MPICC`, `$MPICXX`, `mpicc`, `mpicxx`. Environment values may be
paths or executable names on `PATH`; unavailable candidates are skipped.
Use `plumed mkmpi --compiler /absolute/path/to/mpicc --out provider.so`
(or an MPI C++ compiler driver) to build the shared native adapter. The command
reuses configure-tested compiler/shared-library flags, relies on the `.c`
suffix for language selection, and adds no language-standard flag. The adapter
is C99/C++98-compatible and does not link a PLUMED kernel.
Use a `.dylib` output name on macOS. Pass each additional compiler or linker
argument separately with repeated `--cflag ARG` or `--ldflag ARG` options;
arguments are never shell-evaluated. Either a C or C++ compiler driver works.
The command compiles a small generated stub against the installed public
`Plumed.h`; it does not need adapter source files or the original source
checkout. Failed builds preserve an existing output.

Building a provider does not initialize MPI: the host remains responsible
for its MPI lifecycle. Set `PLUMED_MPI_WRAPPER` to an absolute provider-library path
or an absolute executable MPI compiler path. Host injection
takes precedence; otherwise the environment provider can be selected
explicitly. Embedded tMPI is selected only by the internal launcher. A
loadable library must expose the provider factory. A non-loadable
executable compiler invokes the same `mkmpi` command in a private temporary
directory, with compiler diagnostics inherited by the process. Files with
shared-library suffixes are never executed as compiler commands. The kernel
retains the loaded provider library; temporary build files are removed after
loading. This selection does not initialize the host MPI runtime.

Packages can install a prebuilt adapter beside `libplumedKernel` as
`libplumed-mpi-provider.so` (`libplumed-mpi-provider.dylib` on macOS).
For a custom program name, replace `plumed` with that configured name; for
example, `plumed_custom` uses `libplumed_custom-mpi-provider.so`.
The kernel resolves its actual loaded library path, including symlinks,
rather than using the current Conda
environment or working directory. When a native communicator is supplied,
an installed adapter is tried after injection and explicit environment
selection, but before automatic compilation. A present adapter that cannot
load or fails provider validation raises an error with its path; it does not
silently fall back to a compiler. This discovery neither initializes MPI nor
selects a provider during ordinary construction or capability inquiries.

For troubleshooting, set `PLUMED_MPI_WRAPPER_FORCE=1` together with
`PLUMED_MPI_WRAPPER` before any provider is installed. This makes that
environment-selected provider take precedence over an incoming injection.
Without the force flag, injection retains priority. An absent provider path,
an invalid force value, or a late attempt to replace an installed provider
raises an error. Injected callback code remains host-owned, including when injection is bypassed
by the force setting. Already published callbacks are never replaced.

## Build an MPI executable

`mkmpi --executable` recompiles the exact `src/main/main.cpp` used for the normal
PLUMED executable. The source is also installed, so the original checkout is not
needed. Both modes use the same compiler search: `$MPICC`, `$MPICXX`, `mpicc`,
then `mpicxx`. Executable mode copies the source unchanged into a temporary
`main.c` so a C driver selects C without compiler-specific language flags.
The main uses the C wrapper interface and is C99/C++98-compatible. Neither
mode forces a language standard; either a C or C++ driver works.

```sh
plumed mkmpi --executable --compiler /absolute/path/to/mpicc --out plumed-mpi
mpiexec -n 8 ./plumed-mpi driver --plumed plumed.dat --ixyz trajectory.xyz
```

The default `--link shared` links the current `libplumed`. `--link runtime`
embeds the public wrapper's loader instead, with the current kernel as its
default; `PLUMED_KERNEL` can select a different kernel at execution time.
This allows a precompiled serial kernel to use MPI installed afterwards.
The executable and MPI/kernel dependencies must be available on the compute
nodes. The selected kernel must support injection for this serial-kernel use.

The ordinary main owns MPI initialization/finalization, destroys PLUMED before
finalizing MPI, and passes the native communicator to CLTools. It injects MPI
before other setup when the kernel advertises command API 11 or newer. Older
kernels retain the historical native-MPI path. `--no-mpi` skips MPI initialization,
communicator setup and finalization; `--mpi` restores them if it follows
`--no-mpi`. Neither option changes injection. No additional launcher options are
introduced. PLUMED command exceptions are transferred through the existing C
error interface, including nested messages. On failure the main prints the
diagnostic, releases the error and PLUMED object, and aborts, retaining fatal
exception behavior without entering MPI finalization while peers may be blocked.
Native-MPI kernels ignore injection, including when this same main is built as
the normal native-MPI executable.

To exercise environment loading while using the same executable:

```sh
export PLUMED_MPI_WRAPPER=/absolute/path/to/provider.so
export PLUMED_MPI_WRAPPER_FORCE=1
mpiexec -n 8 ./plumed-mpi driver --plumed plumed.dat --ixyz trajectory.xyz
```

The provider must use the host executable's MPI implementation. Injection is
still requested by the main; the kernel's existing force mechanism overrides it.
Provider-only `mkmpi` builds retain their existing compiler discovery and
C++98-compatible implementation; executable mode does not change them.

## Packaging and conda extensions

Runtime providers allow packaging MPI support separately from the main PLUMED
library. A distributor can build one serial kernel and reuse it with MPICH or
Open MPI, rather than recompiling the full kernel for every MPI variant.
Python wrappers still need their usual Python/platform builds, but do not
need an additional build matrix dimension for the MPI implementation. A
prebuilt provider also removes the end user's need for an MPI compiler or
PLUMED sources when a host passes a native communicator without injection.

PLUMED's own conda recipes prototype this arrangement:

| Package | Role |
| --- | --- |
| `plumed` | One MPI-independent kernel, ordinary CLI, headers and wrappers, with embedded thread-MPI. Uses non-MPI FFTW and has no external MPI runtime dependency. |
| `plumed-mpi`, MPICH variant | Small provider library and MPI-enabled executable, depending on the matching `plumed` version and MPICH. |
| `plumed-mpi`, Open MPI variant | The same two files built against Open MPI, depending on the matching `plumed` version and Open MPI. |
| `py-plumed` | Python wrapper depending on the same base `plumed`, independent of the MPI variant. |

The base keeps the name `plumed`; there is no `plumed-kernel` or `nompi`
add-on. Both MPI implementations are variants of the package `plumed-mpi`, so
only one provider package can be installed at a time. Each installs
`lib/libplumed-mpi-provider.so` (`.dylib` on macOS) next to the kernel and
`bin/plumed-mpi`, built with `plumed mkmpi` against the already installed base.
Installing, changing or removing this add-on does not require replacing the
kernel. The package name also leaves room for future, separately designed
extensions such as `plumed-cuda-extensions`; no CUDA extension is implemented
by this work.

The executable `plumed-mpi` owns external MPI startup and injects its provider.
An existing MD host using an older wrapper can instead supply its native
communicator, triggering discovery of the adjacent provider without an
environment variable or runtime compilation. Both paths must use the same
compatible MPI implementation as the host. The ordinary `plumed` executable
keeps its serial behavior and explicit `--tmpi-threads N` launcher; installing
the add-on does not make that executable call `MPI_Init`. Merely installing
the provider also need not make `plumed --has-real-mpi` succeed before provider
selection. Use the MPI executable or a communicating host to exercise it.

This is a prototype in PLUMED's own conda CI/channel recipes, not a claim that
conda-forge already distributes this layout. See
[`conda/README.md`](../../conda/README.md) for installation examples, the package
matrix and build commands. The package tests cover relocated installations,
injection, adjacent-provider discovery with runtime compilation disabled,
and explicit internal thread-MPI. Fresh-environment checks add and remove
each variant and compare the base kernel bytes. This tests the intended
benefit directly: optional MPI support around one unchanged PLUMED build.

## Source layout and generation

All files live directly in `src/mpiwrap`. The normal PLUMED rules determine
which files are installed and compiled: `.h` files are installed; `.cpp` files
are compiled and linked; templates, metadata and test fixtures are ignored by
the automatic source build. There are no `.hpp` or `.inc` files in this module.

### Installed headers

| File | Interface |
| --- | --- |
| `mpi.h` | MPI-facing umbrella for normal C++ callers. Includes `mpi-types.h` and `mpi-functions.h`; native builds retain the actual MPI interface. |
| `mpi-types.h` | Complete MPI handle/status types and standard `MPI_*` type aliases, without provider dispatch or predefined-value accessors. Used by common headers such as `TypesafePtr.h` and `Communicator.h`. |
| `mpi-functions.h` | MPI function mappings and predefined-value names; includes the type and package interfaces. |
| `types.h` | Shared namespaced transport and handle definitions used by both interfaces. Internal implementation header; does not introduce global MPI aliases. |
| `mpiwrap.h` | Package interface: provider ABI, typed pointers and inline dispatch, capability inquiries, native communicator import and rank launching. Includes `types.h`; does not add provider-mode MPI name mappings. |
| `internal.h` | Shared implementation declarations, principally the registry. Installed by the normal header rules but internal, with no supported external API contract. |

Normal MPI users include `"mpiwrap/mpi.h"` (or
`"plumed/mpiwrap/mpi.h"` from an installed tree). Code needing package services
includes `mpiwrap.h`. Header search paths should name the parent directory,
so the module's `mpi.h` does not shadow the native MPI implementation's
`<mpi.h>`.

Headers that need only MPI types should include `mpi-types.h`; their implementation
files include `mpi.h` or `mpiwrap.h` explicitly when they need MPI calls or package
services. This avoids repeatedly parsing provider tables and inline dispatch in
unrelated translation units, including during static analysis. It does not change
handle layout, runtime selection, or the type checks in `TypesafePtr`. The full
umbrella remains source compatible, and the type and package headers support
either inclusion order. All these headers are generated by `generate.py` and
installed by the normal module header rules.

The C++ package interface requires C++17. Its provider ABI section is also
C-compatible; `PLUMED_MPIW_ABI_ONLY` exposes that section without the C++
frontend, including for C++98 compatibility probes. Host injection and standalone
provider building continue to use the self-contained public `Plumed.h`, whose
adapter remains compatible with C99 and C++98. Complete handle definitions and
inline dispatch stay in headers to preserve type safety and direct call-through.

### Compiled sources

| File | Maintained/generated | Responsibility |
| --- | --- | --- |
| `mpiwrap.cpp` | Generated from `mpiwrap.cpp.in` | Registry validation, dispatch/constant bridges, pointer publication and frontend predefined-value accessors. |
| `loader.cpp` | Maintained | Provider selection, dynamic loading and automatic adapter compilation. |
| `fallbacks.cpp` | Maintained | Shared replacement algorithms. |
| `launcher.cpp` | Maintained | Direct-MPI and embedded thread-MPI launch paths. |
| `tmpi.cpp` | Generated | Thin embedded tMPI adapter and descriptor. |

Each source contains its own build guards. The embedded adapter is compiled
only when thread-MPI is configured; provider-mode registry/loader/fallback
implementations are excluded from native-MPI builds. No `.cpp` file includes
another implementation file. Imported tMPI remains a separate module at
`src/thread_mpi`.

### Generator inputs and checks

`mpiwrap.h.in` and `mpiwrap.cpp.in` contain maintained implementation scaffolding
with explicit generation markers. `provider-abi.h.in` defines the C-compatible
provider ABI. `native-adapter.c.in` contains the maintained native adapter body;
the generator fills its repetitive sections and embeds it directly in
`src/wrapper/Plumed.h`. The same source drives embedded tMPI adapter generation.
There is no separately installed adapter body or separately compiled native
provider source in the kernel module. `mkmpi` builds a small temporary source
against installed `Plumed.h`.

Machine-readable inputs are `support-matrix.tsv`, `argument-metadata.json`,
`constant-metadata.json`, `type-metadata.json`, `signature-versions.json`,
`availability-metadata.json` and `thread-mpi-api.tsv`. Signature revisions are
independent of MPI standard versions. Generated function lists and tables are
written directly into their complete header/source consumers, rather than
emitted as include fragments.

Build the repository formatter once with `make -C astyle`. Regenerate with:

    python3 src/mpiwrap/generate.py

Verify reproducibility without changing files with:

    python3 src/mpiwrap/generate.py --check

Verification scripts and compile-probe fixtures live in
[`rt-mpiwrap-api`](../../regtest/basic/rt-mpiwrap-api/README.md) and
[`rt-mpiwrap-launcher`](../../regtest/basic/rt-mpiwrap-launcher/README.md),
including generated probe `.cpp.in` files. They are compiled explicitly by their
check scripts and consume the same installed headers and public `Plumed.h` as
real callers. The module contains no test executables, scratch logs or generated
include fragments. `support-matrix.md` is the generated human-readable API
support report; validation and coverage checks use the source metadata directly.

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
the host application’s error-handler policy unchanged. Provider-mode custom
reductions are deferred to a later stage and are not a first-release requirement.
`MPI_Op_create` remains cleanly unavailable through providers in the initial
release; direct native MPI calls remain available. A future extension can start
from `mpiw_old_draft`, with a review of callback ABI assumptions and operation
lifetime before publishing callback-capable injection entries. The draft's
callback-reference conversions are useful precedent, not a substitute for that
review. No callback registration or retained-handler allocation occurs for these
unavailable APIs. The mock-native probes compile against real provider headers and exercise call-through without initializing MPI; they are
not real MPI runtime communication tests. Separate host tests and CI exercise
real MPI startup, multi-process communication, injection and environment-provider
loading. The conda package tests also exercise discovery of an installed
provider after relocation. These runtime tests require permission to start MPI
process managers and create their communication sockets; passing a compile or
mock-provider probe alone does not establish runtime correctness.
The pinned thread-MPI source has a one-rank
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
    sh regtest/basic/rt-mpiwrap-api/check-provider-registry.sh
    bash regtest/basic/rt-mpiwrap-api/check-provider-discovery.sh
    bash regtest/basic/rt-mpiwrap-api/check-kernel-lifecycle.sh
    bash regtest/basic/rt-mpiwrap-launcher/check-regtest-runner.sh
    MPICH_CC_WRAPPER=/path/to/mpich/bin/mpicc \
      MPICH_CXX_WRAPPER=/path/to/mpich/bin/mpicxx \
      OMPI_CC_WRAPPER=/path/to/openmpi/bin/mpicc \
      OMPI_CXX_WRAPPER=/path/to/openmpi/bin/mpicxx \
      bash regtest/basic/rt-mpiwrap-api/check-native-adapter.sh

The provider-discovery check uses a mock compiler and provider to exercise
adjacent adapter discovery after relocation and through a kernel symlink,
compiler fallback when absent, rejection of broken installed adapters,
native C/Fortran handle import, injection/environment precedence, and retry
after a compiler failure without starting real MPI.

Header and generated signature checks for local MPICH, Open MPI, and pinned
tMPI use:

    MPICH_CXX=/usr/bin/clang++ OMPI_CXX=/usr/bin/clang++ \
      sh regtest/basic/rt-mpiwrap-api/check-generated.sh \
      /path/to/mpich/bin/mpicxx \
      /path/to/openmpi/bin/mpicxx \
      src/thread_mpi

Internal launch semantics are exercised by
regtest/basic/rt-mpiwrap-launcher/check-thread-mpi-cli.sh src/lib/plumed-static after a configured
thread-MPI build. The tests cover 1, 2, and 4 ranks, callback error handling,
provider lifecycle, raw communicator byte-prefix handling, native/legacy
runner selection, and CHDIR rejection while thread-MPI is active.

## Extending the implementation

Change the maintained metadata/templates, regenerate, then review the complete
outputs. Do not edit generated `.h`/`.cpp` files or the generated regions of
`Plumed.h` directly. The module implements the MPI subset PLUMED needs; adding a
name does not require providing it on every runtime.

### Adding or enabling a function

1. Add the standard signature and scope to `support-matrix.tsv`, and a matching
   function record to `argument-metadata.json`. Use MPI-5 const correctness in
   the frontend. `MPI_Comm_rank` is a small existing example; array/status and
   callback entries show the additional metadata those cases need.
2. Classify every argument: direction, conversion, scalar/array extent,
   sentinels, and when output is valid on success or error. Specify callback
   signatures and lifetime separately. Native handles are provider-sized byte
   values: arrays require native strides, and output storage must remain
   zero-filled beyond the native value. Never cast wrapper arrays to native
   arrays or read native handles as wrapper objects.
3. Record the route for direct MPI, embedded tMPI and injected real MPI. A
   declaration alone does not establish correct tMPI semantics: review the
   implementation and add a reproducer. Keep unsupported/incompatible routes
   unavailable, with a clear exception. To publish a native adapter function,
   mark its `native_adapter` as `implemented` and keep the matching
   `availability-metadata.json` function lists in sync.
4. Set MPI-version requirements in `mpi_standard_functions` and conservative
   GROMACS availability in `availability-metadata.json`. Use MPI version
   numbers for real MPI, not `#ifdef MPI_Function`. For a new tMPI symbol, update
   `thread-mpi-api.tsv` from the audited headers and review its release history
   as described below. The generated body and function table must use identical
   guards.
5. Add a signature record in `signature-versions.json`. The fingerprint is
   calculated by `canonical_signature_fingerprint()` in `generate.py`; the
   validator rejects mismatches. Function IDs are derived deterministically
   from names, not table positions or addresses. A signature revision is not an
   MPI standard version. Changing an existing transport contract requires an
   explicit compatibility design; do not bump a revision merely to silence
   validation and assume already compiled adapters remain compatible.
6. Existing conversion kinds are implemented by
   `render_provider_inline_functions()`, `render_provider_dispatch_bridges()`
   and `render_native_adapter_function()` in `generate.py`. Extend these only
   if the new argument shape is not already supported. Keep native-version
   compatibility casts local to input-only arguments in the native adapter;
   do not weaken the frontend's const correctness or callback types.
7. If a missing function can be implemented using other standard primitives,
   place that algorithm in `fallbacks.cpp`, its needed declarations in
   `mpiwrap.h.in`, and its routing/availability logic in the generator. Select
   it only when the function is absent and all prerequisites are present.
   Preserve error, completion, ownership and communicator semantics; a similar
   prototype or matching numerical result is insufficient.
8. Add semantic calls to `regtest/basic/rt-mpiwrap-api/main.cpp`. Update
   `thread-mpi-contract.json` only after verifying the expected capability;
   lifecycle functions can instead name a real test owner in
   `lifecycle-coverage.json`. Run its `check-coverage.py --write` and review the
   generated unavailable-call and capability checks. Retain tests of clean
   failure for omitted functions, and add an isolated process test for calls
   such as Abort that cannot safely run inside the normal test process.
9. Run `validate-manifest.py`, regenerate, and run `generate.py --check`.
   Exercise the header/prototype checks against MPICH, Open MPI and the
   bundled tMPI, including older-MPI/GROMACS-version guards. Run the mapped-API
   regtest through native MPI, embedded tMPI and the injected/provider host
   paths where that route is supported. The commands in Verification and
   `regtest/basic/rt-mpiwrap-api/README.md` provide the starting points.

The transport remains a versioned, optional name/function table. A frozen host
adapter cannot acquire new functions from a newer kernel: absent entries must
still produce a justified kernel fallback or a clean unavailable exception.
Preserve old entry layouts and signatures when extending the table; incompatible
changes require explicit ABI versioning. Additional compatibility information
belongs in extensible metadata, rather than new mandatory descriptor fields.

### Adding constants, types or a provider-specific workaround

Add a constant to `constant-metadata.json` and `support-matrix.tsv`, specifying
its category, C/native representation and provider routes. Update the constant
availability lists/version overrides in `availability-metadata.json`. Native
values are obtained through accessors; do not hard-code their representation
or eagerly evaluate thread-local communicator constants. Test both successful
publication and omission on runtimes/versions that lack the constant.

Add a new handle/status/callback category to `type-metadata.json` only after
reviewing all scalar, array, import/export and callback conversions it requires.
Handle capacity overrides remain compiler definitions; native builds retain
their actual types, while provider installation validates the receiving storage.
For shared collective fallbacks, C predefined datatypes can use their standard
C types and `sizeof`; custom/compound types need a separate layout design and
currently throw on those fallback paths.

A provider workaround must follow the standard semantics and remain shared
between injected and embedded providers where possible. Record its required
capabilities/compatibility metadata and test the upstream defect independently.
Do not put a second replacement algorithm in the injected `Plumed.h` adapter.
Keep the thin adapter limited to native calls and necessary native-ABI handling.

### Callback extensions: starting point and unresolved choices

Direct native-MPI compilation already uses the native callback types and calls
MPI directly. Provider-mode `MPI_Op_create` and
`MPI_Comm_create_errhandler` remain unavailable today. Their argument metadata
and the generated callback transport/invocation helpers in `mpiwrap.h` provide
a starting point, but do not implement registration or ownership. The old
`mpiw_old_draft` is useful design material; it is not a build dependency or proof
that its callback conversions are portable.

Custom reductions are a possible later-stage extension, outside the initial
release scope. Custom error-handler callbacks
remain outside the intended PLUMED library use: preserve the host's error
handling policy. The pinned tMPI headers do not expose `tMPI_Op_create` or
`tMPI_Op_free`; implementing a real-MPI callback bridge alone will not provide
custom operations on that runtime. Keep that route unavailable until upstream
support or a separately reviewed compatible implementation exists.

There are two approaches worth evaluating before freezing a callback-capable
host adapter:

- **Native typed trampolines.** Compile a trampoline with the host's exact
  `MPI_User_function` signature. Copy its native datatype argument into a
  temporary zero-filled wrapper value, then invoke the saved callback through
  its actual wrapper type. The standard callback supplies no user context, so
  registration needs a way to identify the callback, for example a bounded set
  of explicitly generated trampoline slots. This avoids templates and calling
  through incompatible callback types, but needs slot allocation, concurrency,
  exhaustion and lifetime handling in the native adapter.
- **A documented ABI shortcut based on the draft.** Evaluate whether the
  callback pointer/reference conversions can be supported on an explicitly
  stated platform ABI. Distinguish a tested ABI assumption from C/C++ language
  portability, audit pointee sizes/aliasing as well as argument passing, and
  communicate compatibility through optional metadata. A function-pointer cast
  alone does not establish a portable bridge.

For either approach, resolve these points before publishing an entry:

- Track callback code, state, operation ownership and provider-library lifetime.
  Roll back failed registrations. Do not recycle state while an in-flight
  reduction can still call it, including after an operation is marked for
  freeing. Specify cleanup without asking a PLUMED kernel to finalize real MPI.
- Support concurrent rank threads and independent operations. A process-global
  “current callback” is insufficient. Convert callback-local native handle
  bytes using their reported sizes; never reinterpret a small native handle as
  a larger wrapper object.
- Prevent exceptions from unwinding through MPI's C callback boundary. Decide
  how a callback failure is reported without leaving other ranks blocked.
- Test reduction results with multiple datatypes/ranks, commutative and
  noncommutative operations, concurrent registrations, repeated create/free,
  resource exhaustion, failure rollback and adapter/kernel lifetime. Include
  an older adapter that omits the entry, which must still fail cleanly.

These are extension requirements and options, not promises of current callback
support. Decide the native registration/lifetime mechanism before changing the
frozen host-side surface.

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

GROMACS 2026's tMPI build exposes the long-long datatype constants only when
`SIZEOF_LONG_LONG_INT` is defined. The injected adapter applies that same
compile-time guard to those three constants, so a normal build without
long-long tMPI support remains valid. `check-native-adapter.sh` compiles the
host header in C99 and C++98 modes with this macro both present and absent.

### GROMACS 2026 runtime injection

The `gromacs-2026.0` patch uses the usual top-level `Plumed.h` link (or forwarding
header with `--include`). The integration and its test helper include that header;
GROMACS's bundled copy stays unchanged. Reverting removes the top-level patch
files through the standard patch machinery. Use a PLUMED 2.11 kernel with API
version 11 or newer. The PLUMED kernel can remain serial;
GROMACS supplies its initialized MPI implementation to PLUMED before setup.
The embedded PLUMED tMPI runtime is not involved. Earlier PLUMED kernels retain
the legacy native-MPI path and must be built against the matching MPI library;
multi-rank GROMACS thread-MPI runs are rejected with those older kernels.

For a runtime-linked PLUMED build, apply the patch before configuring GROMACS:

    plumed patch -p --runtime -e gromacs-2026.0
    export PLUMED_KERNEL=/path/to/libplumedKernel.so
    gmx mdrun -ntmpi 4 -ntomp 2 -plumed plumed.dat

On macOS the kernel library has the `.dylib` extension.

For native MPI, configure GROMACS with `-DGMX_THREAD_MPI=OFF -DGMX_MPI=ON` and
launch the MPI executable with its matching process manager:

    mpiexec -n 4 gmx_mpi mdrun -ntomp 2 -plumed plumed.dat

The integration smoke matrix was checked on 2026-11-07 using GROMACS tag
`v2026.4`, its own thread-MPI, MPICH 5.0.1 and Open MPI 5.0.10, all loading
the same serial PLUMED kernel. All three passed at 1, 2 and 4 ranks, plus
4 ranks with 2 OpenMP threads each. It
checked an analytic distance restraint, PLUMED generalized forces, atom force
increments against a no-PLUMED baseline, and rank-count consistency. This is a
focused integration check, not the full GROMACS upstream test suite.

Patch/revert checks also passed for tags `v2026.0` through `v2026.4`, both with
the default links and with `--include`, restoring every original source file
byte for byte. Fresh patched and unpatched `v2026.4` thread-MPI builds passed:
after reverting, single-rank PLUMED use still worked and the original
multi-rank thread-MPI rejection returned, confirming that injection was removed.

The capacity-override installation check uses a configured checkout:

    ./configure --disable-mpi
    MPICH_CXX=/usr/bin/clang++ OMPI_CXX=/usr/bin/clang++ \
      bash regtest/basic/rt-mpiwrap-api/check-capacity-install.sh /path/to/mpich/mpicxx \
      /path/to/openmpi/mpicxx

It runs the public install-headers recipe and compiles installed headers
with explicit compiler definitions for 24-byte handles and 192-byte status
storage. The build/runtime matrix above remains separate validation.

A separate MPI host checks injection into an already-built serial kernel:

```sh
MPICH_CXX=/usr/bin/clang++ bash regtest/basic/rt-mpiwrap-api/check-mpi-host-injection.sh \
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
serially. `plumed mkmpi --executable` builds the injecting host from the
installed standard main. Focused checks exercise 1, 2 and 4 ranks in shared and
runtime modes. Debug and Intel jobs add a posterior Open MPI injection pass
after their normal suites. The focused checks also build a provider DSO with
`plumed mkmpi` and select it through `PLUMED_MPI_WRAPPER` together with
`PLUMED_MPI_WRAPPER_FORCE=1`. The same checks can be repeated manually with
another MPI implementation.
The same serial kernel is tested again through injection;
there is no PLUMED rebuild. The host owns MPI initialization
and finalization. `PLUMED_CLI_PROGRAM_NAME` selects it for CLTool tests only,
and the regression runner derives the installed library paths. External-engine
tests retain their own initialization and capability
checks, so MPI cases requiring a natively built external engine remain skipped.
No GROMACS build is required. The standalone embedded tMPI tests also run.

Run the same host test locally with an already built serial kernel:

```sh
bash regtest/basic/rt-mpiwrap-api/check-mpi-host-injection.sh /absolute/path/to/mpicxx \
  /absolute/path/to/mpiexec /absolute/path/to/libplumedKernel.so \
  /absolute/path/to/provider.so
```

The optional fourth argument enables environment-provider cases in addition to
injection cases. Omitting it retains the injection-only test. Both paths use
unmodified standard-main behavior; `--no-mpi` queries avoid MPI lifecycle calls.


## Integration with existing PLUMED code

The MPI provider implementation keeps existing callers and native-MPI initialization
unchanged. The standard main now requests injection before CLTool setup and is
reused by `mkmpi --executable`. Kernel initializer/symbol-table code and existing
wrapper entry-point implementations are unchanged. There is no new kernel
symbol-table entry or public object type.

| Files | Reason for the feature change |
| --- | --- |
| `configure.ac`, generated `configure` | Discover native MPI as before; include embedded tMPI by default in serial builds, with an explicit opt-out. Run it only through the explicit launcher. Define the wrapper infrastructure independently of runtime selection. No capacity configure options or module `.on`/`.off` manipulation. |
| `CMakeLists.txt` | Enable the wrapper and optional embedded Win32/POSIX thread backend; retain registration objects with each platform's linker syntax. Existing Windows configuration remains supported. |
| `src/.gitignore`, `src/{mpiwrap,thread_mpi}/.gitignore` | Use the existing per-module ignore conventions to track sources, generator/import tools and metadata while ignoring build artifacts. |
| `src/{core,generic,tools}/Makefile` | Declare the new module/header dependencies for normal builds and installed headers. |
| `src/config/Makefile` | Export native-provider compiler/link flags for `mkmpi`; regenerate options when this Makefile changes. The provider uses a C ABI and does not link a C++ PLUMED runtime. |
| `src/maketools/update-config-txt.sh` | Separate configured real-MPI and internal-tMPI capabilities while retaining combined `has mpi`. |
| `src/core/CLToolMain.cpp` | Parse internal rank counts, invoke the shared launcher, construct rank-local commands and report selected-provider runtime capabilities. |
| `src/core/PlumedMain.cpp` | Expose injection through the existing command mechanism and advertise command API version 11. Construction still permits injection before provider-dependent use. |
| `src/tools/Communicator.{h,cpp}` | Replace compile-time serial branches with runtime initialization checks, preserve serial fallbacks, import native bytes through the single boundary helper and resolve default SELF before communication. An explicit flag distinguishes deferred SELF from a valid all-zero native handle. Private buffer descriptions use const references for fatal cppcheck diagnostics; receive buffers remain writable. The header includes only MPI types; the implementation includes the full MPI interface. |
| `src/tools/TypesafePtr.h` | Include the MPI type interface for communicator checks without importing provider tables or dispatch. The checking implementation is unchanged. |
| `src/tools/PlumedHandle.cpp` | Forward the selected implementation to a compatible nested kernel through the existing command API; preserve older-kernel behavior and release the created reference on failure. |
| `src/generic/Plumed.cpp` | Reject process-global CHDIR while a thread-MPI runtime is initialized. Communicator forwarding at existing call sites is unchanged. |
| `src/core/Value.h` | Suppress only the observed cppcheck false positive on a member-reference accessor, without changing its code or return type. |
| `src/maketools/codecheck` | Exclude the imported runtime like other vendor code, avoid opt-in adapter configurations in kernel analysis, and exempt only structural rules incompatible with generated C ABI headers. Semantic diagnostics remain enabled. |
| `src/wrapper/Plumed.h` | Add the opt-in generated host/native adapter and its one implementation accessor. Existing C/C++ wrapper implementations and exception handling are unchanged. |
| `scripts/mkmpi.sh` | Build a standalone provider from the public header or recompile the installed standard main into a shared/runtime MPI executable, with atomic output publication. |
| `src/main/main.cpp`, `src/lib/Makefile` | Request injection on MPI-enabled entry before CLTool setup, retaining older-kernel behavior. Install this exact source for `mkmpi --executable`; no second CLI host implementation. |
| `plugins/pycv/FindPlumed.cmake` | Link external MPI only for a native-MPI build; embedded tMPI does not require external MPI headers or libraries. |
| `regtest/scripts/run`, `mpi-test.h` | Select direct MPI, an injecting CLI host or internal ranks while keeping old `mpiprocs` values and old-kernel behavior. External MPI lifecycle remains owned by the executable. |
| Converted existing MPI regtests | Mechanical lifecycle/communicator changes only; original calculations, rank counts and reference files are retained. Four cases with process-global environment or independent serial processes remain native-only. |
| New API/launcher/provider regtests | Check the mapped API, unavailable calls, transport, lifecycle, nested kernels without compiling GROMACS. |
| `.ci/check-warnings`, `.codecov.yml` | Treat the imported thread-MPI runtime like other vendored libraries for warning ownership and coverage metrics; the PLUMED MPI adapter remains checked. |
| Linux/Windows workflows | Keep master's existing jobs, including sanitizer/no-OpenMP/no-dlopen checks. Add one serial-build/full-injection job, posterior OpenMPI passes for debug and Intel, focused provider tests and explicit Windows smoke checks. |
| `conda/` recipes and conda workflow | Build one MPI-independent `plumed`, optional `plumed-mpi` provider/executable variants, and the existing Python variants. Test relocation, automatic discovery, injection and add-on removal without changing the kernel. |
