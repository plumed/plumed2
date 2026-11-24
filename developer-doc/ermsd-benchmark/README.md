# ERMSD optimization benchmark

This experiment compares the original ERMSD, a streamed all-pairs implementation,
and linked-cell and buffered-list implementations on 1S72-derived RNA structures.
Results and machine details are recorded in `results.md`.

Python 3.10 or newer is required. Generated structures, trajectories, logs and
kernel snapshots should be kept outside the checkout. No Python packages are
required. `generate.py` downloads 1S72 if absent, selects complete C2/C4/C6 RNA
base triplets in the documented purine/pyrimidine order, and writes deterministic
120-frame XYZ trajectories in nm and corresponding PDB references. Sizes are
nested prefixes of 2871 complete bases; larger systems contain additional copies
translated 25 nm along x. Smooth excursions of amplitude 0.08 nm per component
are shared within each triplet, with additional 0.004 nm atomic distortions.
The default 0.2 nm neighbor skin is therefore crossed repeatedly.

## Local build and kernels

The measured upstream base is `08a6de70c0820d5dfb7b3bc48a8b2f9d959a9684`.
Build the baseline, the first streamed optimization commit, and the final commit
with identical compilers, options and module settings. Keep the resulting
`src/lib/libplumedKernel.dylib` (macOS) or `.so` (Linux) under the output directory
as `kernels/baseline.dylib`, `kernels/direct.dylib`, and `kernels/spatial.dylib`.
The snapshot filenames are accepted by the dynamic loader on either platform.
The first streamed implementation is commit `7c53afefe` on this task branch.

The local configuration used:

```sh
./configure --disable-mpi --enable-thread-mpi --disable-libsearch \
  --enable-modules=reset --prefix="$PWD/../.local-thread" \
  CXX=/opt/local/bin/clang++-mp-18 CC=/opt/local/bin/clang-mp-18 \
  CXXFLAGS='-O3 -g'
make -j8
```

OpenMP must be detected (`-fopenmp` here). PLUMED's thread count is controlled by
`PLUMED_NUM_THREADS`; setting only `OMP_NUM_THREADS` does not enable multiple
PLUMED threads. Both variables are set by the runner.

This upstream version can use an injected MPI provider with the same kernel:

```sh
plumed mkmpi --executable --link runtime \
  --compiler /opt/local/bin/mpicc-openmpi-mp --out "$ERMSD_BENCHMARK_DIR/plumed-openmpi"
```

Use matching Open MPI compilers and launchers; the measured configuration used
Open MPI 5.0.11. The runner uses `--host localhost:N --bind-to none -n N`.
Embedded thread-MPI is launched with `--tmpi-threads N`.

## Generate, validate and benchmark

Set output and executable locations explicitly, for example:

```sh
export ERMSD_BENCHMARK_DIR=/absolute/path/to/ermsd-data
export PLUMED_BENCHMARK_EXE=/absolute/path/to/plumed/src/lib/plumed-runtime
export MPIEXEC=/opt/local/bin/mpiexec-openmpi-mp
export PLUMED_OPENMPI_EXE="$ERMSD_BENCHMARK_DIR/plumed-openmpi"
python3 developer-doc/ermsd-benchmark/generate.py
python3 developer-doc/ermsd-benchmark/validate.py
python3 developer-doc/ermsd-benchmark/run.py --suite serial --repeats 3
python3 developer-doc/ermsd-benchmark/run.py --suite parallel --repeats 3
```

Run without other builds/tests or competing benchmarks. The runner performs
sequential trials, saves complete logs and exact commands, and records the
external benchmark B2+B3 elapsed times divided by their combined cycle counts.
These exclude initialization, first step, warm-up and trajectory parsing, but
include sharing, ERMSD, restraint-force application and other per-step PLUMED
overhead. Three independent process trials are summarized by their median.

Each normal trial has a three-second wall-clock budget. The original kernel at
1000/3000 bases instead runs 120 frames (maximum 120 seconds), ensuring complete
trajectory coverage despite large allocation costs. The 6000-base baseline is
intentionally skipped: its temporary derivative allocation alone needs 31.1 GB,
exceeding this machine's 24 GiB physical RAM.

The validation compares values and full atomic/box derivatives against the
original kernel over all 120 frames, including selected `PAIRS`, OpenMP,
thread-MPI, Open MPI and hybrid runs. It also checks analytical derivatives
against PLUMED's numerical differentiation on six bases.

Additional checks and full-trajectory repetitions used in the report:

```sh
python3 developer-doc/ermsd-benchmark/validate_extra.py
python3 developer-doc/ermsd-benchmark/run.py --suite serial --label coverage \
  --sizes 3000 6000 --methods streamed DIRECT --configs serial \
  --nsteps 120 --maxtime 120 --repeats 3
python3 developer-doc/ermsd-benchmark/run.py --suite parallel --label domain \
  --sizes 3000 6000 --methods VERLET \
  --configs tmpi2 tmpi4 mpi2 mpi4 tmpi2x2 mpi2x2 \
  --domain-decomposition --repeats 3
python3 developer-doc/ermsd-benchmark/summarize.py
```

`validate_extra.py` checks cutoff crossings, contact disappearance and coincident
centers against the original kernel or finite differences, then compares all
search methods and MPI/OpenMP on 6000 nucleotides. The large derivative dumps
use a sparse output stride to keep validation I/O bounded. The existing basic
regressions also cover the default make-whole behavior and `NOPBC`.

The `coverage` repeats replace the short serial all-pairs trials at 3000 and
6000 nucleotides in the summary, ensuring all 120 frames contribute. The
`domain` trials simulate distributed atom ownership and atom shuffling, in
addition to the default replicated-coordinate layout in the main scaling sweep.
Raw logs and trial data remain in the output directory; `timings.csv` records
medians, extrema and sample counts across the three process trials.
