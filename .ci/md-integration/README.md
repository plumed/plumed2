# MD integration CI

The `MD integration` workflow builds serial PLUMED once on Ubuntu 24.04 and
macOS 15 (Apple Silicon). Each installation is shared as an artifact with the
GROMACS MPI matrix and with two MPI plus one serial job for each of LAMMPS and
AmberTools sander. PLUMED, GROMACS and LAMMPS have OpenMP enabled; sander uses
one thread per rank or runs as a serial executable.
Artifacts preserve installed symlinks and permissions and use the same absolute
installation prefix in producer and consumer jobs. They are not shared between
operating systems or architectures.
Linux jobs add the private installation's `lib` directory to `LD_LIBRARY_PATH`
before using the installed executable, including after downloading the artifact.
OpenMP flags are passed both when compiling and linking the kernel. Before
uploading an artifact, `check-kernel.sh` loads it from a plain C executable
without MPI/OpenMP link flags. This catches missing runtime dependencies that
an OpenMP-enabled MD host could otherwise hide.

All source versions are pinned in [`versions.json`](versions.json). The installed PLUMED patcher applies,
reverts and reapplies the distributed `gromacs-2026.0` patch. The revert must
restore the tracked GROMACS sources exactly and remove the top-level wrapper
files. The patched build uses the installed current `Plumed.h`. GROMACS's bundled
wrapper remains untouched. Each job installs/selects its own host MPI:

- GROMACS thread-MPI;
- Open MPI;
- MPICH.

Ubuntu 24.04's packaged MPICH uses external PMIx and cannot form a shared MPI
world with its bundled Hydra launcher ([MPICH issue #7064](https://github.com/pmodels/mpich/issues/7064)).
The Linux MPICH job therefore uses Open MPI's PMIx-capable launcher, while
GROMACS and the injected provider still compile/link against MPICH. UCX uses
TCP, shared memory and self transports, avoiding unavailable hosted-runner RDMA
devices. This reuses the workaround already validated in PLUMED's injection
tests. On macOS, MPICH uses its own launcher. Every real-MPI job first runs the
existing native MPI probe with four ranks, without GROMACS or PLUMED, so a
broken MPI installation fails before the expensive GROMACS build.

The workflow is CPU-only and uses GROMACS's bundled FFTPACK to avoid installing
an FFTW package that can bring in a different MPI. This input uses cutoff
electrostatics, so FFT performance is irrelevant. Full GROMACS compilation is
intentional; GPU builds and GROMACS's own regression suite are outside this test.

## Running the workflow

MD builds are optional and run only through `workflow_dispatch`. In GitHub
Actions, select **MD integration**, choose **Run workflow**, and select the
branch to test. Ordinary push, pull-request and scheduled CI remain enabled
independently; the MD workflow does not suppress or replace any existing job.

The complete matrix has 20 jobs: two PLUMED artifact builds, six GROMACS jobs,
six LAMMPS jobs and six AmberTools jobs. All 20 jobs passed on Linux and macOS
in [the integration validation run](https://github.com/GiovanniBussi/plumed2/actions/runs/42368828147),
including all four serial AmberTools/LAMMPS jobs. Both serial engines also
passed local Apple Silicon builds and numerical validation, including LAMMPS
with one and two OpenMP threads.

The reusable `md-gromacs.yml`, `md-lammps.yml` and `md-amber.yml` are called by the parent after
the PLUMED artifacts are ready. Further engines can have separate reusable
workflows sharing those artifacts. LAMMPS and AmberTools injection/header
edits stay local to these tests, rather than adding entries to the distributed
`plumed patch` collection.

## Updating engine versions

Edit **only `versions.json` to change a source pin**. The workflows read it
through `version.py`; no commit SHA, download URL or archive checksum is duplicated
in workflow files. Git engines record a readable release/tag and its exact
40-character commit SHA. AmberTools records its official archive URL and SHA256.
The GROMACS entry also selects the compatible distributed PLUMED patch.

For GROMACS or LAMMPS, fetch the new release in a local engine checkout, then
use `git rev-parse 'RELEASE_TAG^{commit}'` to obtain its commit (including for
annotated tags). Update `version`, `tag` and `ref` together. For AmberTools,
download the new archive and use `shasum -a 256 ARCHIVE` to update `url`,
`version` and `sha256`. Keep exact hashes rather than moving branch references.

Inspect an entry locally with, for example:

```sh
python3 .ci/md-integration/version.py gromacs
```

Before accepting an update, run the MD integration workflow and inspect every
backend's numerical summary. The jobs apply/reverse/reapply the relevant patch,
so source changes that require adapting a hook fail visibly. Update the small
test-only LAMMPS/Amber patch if its injection point moved; a GROMACS major release
may need a new distributed patch. Also check any fixture paths inside Amber's
archive and keep its C-only injection hook minimal. Uploads retain logs and
results for diagnosis. Version updates intentionally do not regenerate PLUMED's
embedded thread-MPI module: that has its own importer and version policy.

## Numerical and injection checks

`gromacs/` contains the topology, a 20-step MD input, PLUMED input, and a portable
Python runner. Before running `grompp`, the runner copies the unchanged 5832-atom
argon configuration from `src/testutils/simulationdatabase/argon5832.gro` in the
pinned GROMACS checkout into its results directory. The large coordinate file
and precompiled TPR are not duplicated in PLUMED: each backend prepares its own
TPR. Outputs remain in the requested results directory, including the exact
input files, command logs, trajectories and `summary.json`.

For every backend the runner tests 1, 2 and 4 ranks with one OpenMP thread each,
then **4 ranks with 2 OpenMP threads per rank**. The hybrid case explicitly sets
`-ntomp 2`, `OMP_NUM_THREADS=2` and `PLUMED_NUM_THREADS=2`; it is not left to CPU
count heuristics. The runner verifies the rank/thread counts in GROMACS's log
and that PLUMED receives those counts. GROMACS owns rank startup and MPI
initialization; PLUMED's internal thread-MPI launcher is not used.
The Open MPI launcher receives `--oversubscribe --bind-to none` explicitly, so
four ranks and hybrid threads can run on smaller hosted runners without relying
on version-dependent MPI environment settings. Other backends receive no Open
MPI-specific launcher options unless explicitly using Open MPI's launcher for
the PMIx workaround above.

Every biased run sets `PLUMED_MPI_WRAPPER` to a nonexistent path, disables forced
provider selection and sets `MPICC`/`MPICXX` to `/usr/bin/false`. Success therefore
requires the injected host implementation to take precedence over discovery;
the test cannot silently pass by compiling or loading another provider.

Checks cover all 21 output frames:

- finite output and complete, ordered time steps;
- analytical harmonic bias `50 * (distance - 0.3)^2` and generalized force
  `-100 * (distance - 0.3)`;
- rank/hybrid consistency of distances, bias, potential energy and positions;
- actual atomic force increments, from the first trajectory frame compared
  with a separate four-rank simulation without PLUMED.

The existing numerical tolerances from the local tests are retained: `2e-8`
for analytical bias/generalized force, `2e-6` nm for distance, `2e-4` kJ/mol for
bias, `2e-5` nm for positions, and `5e-6` relative to the initial potential energy
for energy reductions. The atomic-force comparison allows `0.02` kJ/mol/nm for
the text precision of `gmx dump`. Differences are reported rather than compared
against platform-specific reference files. Logs/results are uploaded even on
failure, with a seven-day retention period.

## Running locally

Use a serial PLUMED kernel and a GROMACS executable patched/built with the
current header. The output directory must not exist yet. From the PLUMED root:

```sh
python3 .ci/md-integration/gromacs/check.py --backend tmpi \
  --gromacs-source /path/to/gromacs-source \
  --gmx /path/to/gromacs-build/bin/gmx \
  --kernel /path/to/libplumedKernel.dylib --results /tmp/gromacs-tmpi-results
```

For Open MPI or MPICH, select the corresponding backend, use `gmx_mpi`, and add
`--launcher /path/to/the/matching/mpiexec`. Ensure the MPI's runtime libraries
and process-manager helpers are accessible in the environment. On Linux use
`libplumedKernel.so`. `--check-only` rechecks an existing result directory without
rerunning the simulations (it still invokes `gmx dump` to read stored forces).
When launching Ubuntu's packaged MPICH with Open MPI, also pass
`--launcher-kind openmpi` and set `UCX_TLS=tcp,sm,self`, as the workflow does.

`build-gromacs.sh SOURCE BUILD` builds an already patched tree. Set `BACKEND`
to `tmpi`, `openmpi` or `mpich`, and optionally `CC`, `CXX`, `MD_BUILD_JOBS`.
Real MPI additionally requires `MD_MPICC` and `MD_MPICXX`; Apple Clang OpenMP
requires `MD_LIBOMP` pointing to a prefix containing `include/omp.h` and
`lib/libomp.dylib`. The shared `setup-mpi.sh` is the GitHub-runner dependency
setup and writes selected paths to `GITHUB_ENV`/`GITHUB_PATH`. `check-mpi.sh`
compiles/runs the native MPI probe before either engine is built.

## LAMMPS

LAMMPS's stable release is selected by the manifest. The matrix runs Open MPI,
MPICH, and a true non-MPI serial build on both operating systems, using the
same serial PLUMED artifacts as GROMACS. LAMMPS does not have a GROMACS-style
thread-MPI backend.

The test copies the installed current `Plumed.h` to the LAMMPS source root.
The small test-only `lammps/injection.patch` changes `fix_plumed.cpp` to include
that copy, enables `PLUMED_MPIW_HOST_ADAPTER`, and calls `injectMPIImplementation`
after creating the PLUMED object and querying its API. Older kernels receive
one main-rank warning to use the same MPI library as LAMMPS or upgrade PLUMED;
their original communicator setup continues. No other LAMMPS code is changed.
CI verifies that reversing the patch restores all tracked sources. CMake uses
the installed PLUMED runtime wrapper (`PLUMED_MODE=runtime`), disables downloading
another PLUMED, and enables LAMMPS's `OPENMP` and `PLUMED` packages.

The runner generates a 256-atom periodic FCC argon crystal with fixed IDs and
zero initial velocities. IDs are assigned before decomposition so they do not
depend on how many ranks read the system. A 20-step trajectory uses the threaded
`lj/cut/omp` pair style, selected explicitly with `-sf omp -pk omp N`. A harmonic
restraint between atoms 1 and 2 has target 3.5 angstrom and force constant
0.1 eV/angstrom squared. PLUMED uses matching angstrom/eV/ps units.

As for GROMACS, each job tests 1, 2 and 4 MPI ranks with one OpenMP thread, then
4 ranks with two threads each, and a four-rank baseline without PLUMED. The
MPI compiler overrides and nonexistent provider path rule out automatic provider
discovery as a substitute for injection. After every simulation the runner
checks LAMMPS's processor grid and MPI/OpenMP counts. It also checks the rank
and thread counts received by PLUMED and verifies that the threaded pair style
was active. Singleton worlds therefore cannot silently pass.

Checks cover all 21 frames: ordered time steps and atom IDs, finite values,
distance measured independently from the native trajectory, analytical bias
and generalized force, and rank/hybrid consistency of all positions and forces.
The first frame is compared with the no-PLUMED baseline to verify the actual
force increment on both restrained atoms and zero increments on every other
atom. `ENERGY` is compared with LAMMPS's native potential; `fix_modify plumed
energy no` keeps this native potential independent of the bias. Tolerances are
`2e-10` for distance/bias/generalized force, `2e-9` eV for native energy and
rank comparisons of COLVAR, `2e-8` for rank comparisons of positions/forces,
and `2e-9` eV/angstrom for the force increments. Inputs, exact commands, native
logs, atom trajectories and `summary.json` are retained and uploaded.

To build locally, first copy the current header and apply the test patch:

```sh
cp /path/to/plumed/include/plumed/wrapper/Plumed.h /path/to/lammps/Plumed.h
patch -d /path/to/lammps -p1 < .ci/md-integration/lammps/injection.patch
BACKEND=mpich MD_MPICXX=/path/to/mpicxx \
  bash .ci/md-integration/build-lammps.sh /path/to/lammps /path/to/build /path/to/plumed
python3 .ci/md-integration/lammps/check.py --backend mpich \
  --lmp /path/to/build/lmp --launcher /path/to/mpiexec \
  --kernel /path/to/libplumedKernel.dylib --results /tmp/lammps-mpich-results
```

Use `--backend openmpi` for Open MPI and `.so` for the kernel on Linux. The
same Apple Clang OpenMP settings (`MD_LIBOMP`) and Ubuntu MPICH launcher
workaround (`--launcher-kind openmpi`, `UCX_TLS=tcp,sm,self`) apply as described
above. `--check-only` validates existing outputs without launching MPI again.

The complete local matrix passed on Apple Silicon with MPICH 5.0.1 and Open MPI
5.0.10 against a serial OpenMP-enabled PLUMED kernel, including the four-rank,
two-thread case and native force/energy comparisons.

The serial build sets `MD_ENGINE_MODE=serial` and `BUILD_MPI=OFF`. Its runner
executes the binary directly, tests one and two OpenMP threads, and checks the
serial PLUMED setup, analytic bias/forces, native energy, and no-PLUMED baseline.
It does not launch MPI or inject an MPI provider. LAMMPS's serial MPI stubs
still report one task, which the runner checks alongside PLUMED's serial setup.

```sh
MD_ENGINE_MODE=serial bash .ci/md-integration/build-lammps.sh \
  /path/to/lammps /path/to/build /path/to/plumed
python3 .ci/md-integration/lammps/check.py --backend serial --mode serial \
  --lmp /path/to/build/lmp --kernel /path/to/libplumedKernel.dylib \
  --results /tmp/lammps-serial-results
```

## AmberTools sander

The AmberTools source archive and SHA256 are selected by `versions.json`.
`md-amber.yml` builds `sander.MPI` with Open MPI and MPICH, plus the ordinary
non-MPI `sander` target in serial jobs, on both operating systems. All use the
same serial PLUMED artifacts as the other engines. Python, CUDA, automatic
source updates and Amber's own test installation are disabled. CMake 3 is used
for Amber's bundled dependencies. MPI jobs use one thread per rank; serial jobs
run without an MPI launcher or synthetic rank count.

The test separately replaces `AmberTools/src/sander/Plumed.h` with the installed
header. The small `amber/injection.patch` contains only a C helper, its Fortran
declaration and a call after global-object creation. Existing `plumed=1` handling
guards that initialization, and the helper checks that the global object is
valid and supports API 11 before injecting. Older kernels retain their original
communicator setup. There are **no Amber CMake changes**.

The C source enables injection when its compiler can find `mpi.h`; using `cc`
with matching include/link flags works as well as using `mpicc`. If header
detection is unavailable or `mpi.h` is absent, compilation emits a warning and
the helper returns zero without accessing PLUMED. Fortran then writes a warning
to Amber's output on the master rank and proceeds with the original interface.
This path does not inject an MPI implementation: PLUMED must use the same MPI
library as sander through its existing native/provider configuration. Compilers
without `__has_include` can explicitly define `PLUMED_MPIW_HOST_ADAPTER` when
the matching C MPI headers and libraries are available. An MPI-enabled adapter
still needs the native MPI C symbols at link time; the tested MPI Fortran
libraries already supply the matching `libmpi` dependency.
Separate helper statuses distinguish missing C adapter support, an older
kernel API, and an invalid global object. Amber's master rank reports the
corresponding warning; the older-API warning recommends upgrading the kernel
rather than recompiling sander. A successful injection produces no warning.

`NMA.prmtop` and `NMA.inpcrd` are read directly from the pinned archive's
`AmberTools/test/paramfit/simplex_perfect_fit_NMA/prmtop_creation` fixture. The small
12-atom NMA system runs five steps in vacuum with a distance restraint on atoms
1 and 2 (`AT=0.2` nm, `KAPPA=1000` kJ/mol/nm squared). Inputs write ASCII
coordinates, forces and restart files, requiring no additional Python packages.
Sander evaluates the initial PLUMED step twice, so COLVAR has six rows with
times 0, 0, 0.001, 0.002, 0.003 and 0.004 ps.

For MPI jobs at 1, 2 and 4 ranks, the runner checks native Amber and PLUMED rank
counts, finite/complete output, analytic harmonic bias and squared generalized
force, and rank consistency of COLVAR, positions and forces. Serial jobs invoke
`sander` directly, verify the absence of the AMBER/MPI banner, and require
PLUMED's serial setup and analytic CV/bias/force checks. Both modes run separate
`plumed=0` baselines, which use an invalid kernel path and must not load PLUMED
or create COLVAR.
The first biased/native force difference is checked on every atom against the
analytic restraint force, with conversion from PLUMED's kJ/mol/nm to Amber's
kcal/mol/angstrom. The force tolerance, `0.0011` kcal/mol/angstrom, covers two
native dumps rounded to three decimals. A nonzero final coordinate change
also checks that the bias affects dynamics. Injection runs block provider
discovery exactly as the other engines do. Inputs, native logs, force and
coordinate trajectories, command lines and numerical summaries are uploaded.

For a local build, copy the new header and apply the small test patch, then use
ordinary compilers and the matching MPI wrappers:

```sh
cp /path/to/plumed/include/plumed/wrapper/Plumed.h /path/to/amber/AmberTools/src/sander/Plumed.h
patch -d /path/to/amber -p1 < .ci/md-integration/amber/injection.patch
CC=cc CXX=c++ FC=gfortran MD_MPICC=/path/to/mpicc \
  MD_MPICXX=/path/to/mpicxx MD_MPIFC=/path/to/mpifort \
  bash .ci/md-integration/build-amber.sh /path/to/amber /path/to/build
python3 .ci/md-integration/amber/check.py --launcher-kind mpich \
  --amber-source /path/to/amber \
  --launcher /path/to/mpiexec \
  --sander /path/to/build/AmberTools/src/sander/sander.MPI \
  --kernel /path/to/libplumedKernel.dylib --results /tmp/amber-mpich-results
```

Use `--launcher-kind openmpi` for Open MPI (also for Ubuntu's MPICH/PMIx
workaround), and `.so` for the kernel on Linux. The default matrix tests all
three rank counts; `--ranks 1` limits a local smoke run to one rank. Real MPI
launchers must run outside any sandbox that restricts process-manager sockets.

For a serial build, omit MPI wrappers and launchers:

```sh
MD_ENGINE_MODE=serial CC=cc CXX=c++ FC=gfortran \
  bash .ci/md-integration/build-amber.sh /path/to/amber /path/to/build
python3 .ci/md-integration/amber/check.py --mode serial \
  --amber-source /path/to/amber \
  --sander /path/to/build/AmberTools/src/sander/sander \
  --kernel /path/to/libplumedKernel.dylib --results /tmp/amber-serial-results
```
