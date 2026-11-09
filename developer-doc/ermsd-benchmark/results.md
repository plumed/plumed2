# ERMSD optimization: measured results

2026-11-09. All times below are **milliseconds per PLUMED step**, median of three
independent process trials. The original upstream base is
`08a6de70c0820d5dfb7b3bc48a8b2f9d959a9684`; the first streamed all-pairs stage is
`7c53afefe`. Final source, regression cases and reproduction scripts are on
`agents/ermsd/optimization`. See [reproduction instructions](README.md).

At 3000 nucleotides, the buffered implementation with four OpenMP threads is
approximately **549× faster** than the original serial ERMSD on this input.
The main gain comes from eliminating quadratic derivative storage and then
replacing all-pairs search with an exact sparse contact calculation.

## Serial implementation comparison

| Nucleotides | Original | Streamed all-pairs | Sparse all-pairs | Linked cells | Buffered list | Buffered, OMP4 |
|---:|---:|---:|---:|---:|---:|---:|
| 100 | 0.1504 | 0.08967 | 0.08883 | 0.08457 | 0.06614 | not measured |
| 300 | 3.638 | 0.5221 | 0.4867 | 0.2962 | 0.2138 | not measured |
| 1000 | 62.19 | 4.606 | 4.118 | 1.224 | 0.8334 | 0.3442 |
| 3000 | 577.9 | 38.66 | 33.46 | 4.237 | 2.82 | 1.052 |
| 6000 | not run | 150.2 | 130.7 | 8.507 | 5.678 | 2.097 |

The original 6000-base calculation was intentionally not run: its temporary
`(3N)^2` array of 4×3 derivative tensors alone requires 31.1 GB, exceeding this
machine's 24 GiB physical RAM. No timings are extrapolated.

“Streamed all-pairs” is the first optimization stage, retaining the original
all-pairs algorithm and dense reference matrix, while accumulating derivatives
immediately. “Sparse all-pairs” is final `PAIR_SEARCH=DIRECT`, with sparse
reference storage and early rejection. `CELLS` rebuilds PLUMED's linked-cell
NeighborList every evaluation; `VERLET` buffers that list and updates it by
center displacement. Serial CELLS/VERLET and OMP4 entries were remeasured after
the final source rebuild; slow serial all-pairs cases use all 120 input frames.

## OpenMP, buffered list

| Nucleotides | 1 thread | 2 threads | 4 threads | 8 threads |
|---:|---:|---:|---:|---:|
| 1000 | 0.8334 | 0.5353 | 0.3442 | 0.4029 |
| 3000 | 2.82 | 1.667 | 1.052 | 1.096 |
| 6000 | 5.678 | 3.328 | 2.097 | 2.078 |

Four threads are the most consistently effective setting on this host. The
8-thread gain at 6000 bases is small and within the observed trial spread.
`PLUMED_NUM_THREADS` and `OMP_NUM_THREADS` were both set explicitly.

## MPI, buffered list, replicated coordinates

| Nucleotides | thread-MPI 2 | thread-MPI 4 | thread-MPI 8 | Open MPI 2 | Open MPI 4 | Open MPI 8 |
|---:|---:|---:|---:|---:|---:|---:|
| 1000 | 0.6044 | 0.79 | 3.907 | 0.5714 | 0.6064 | 2.693 |
| 3000 | 1.92 | 2.069 | 7.944 | 1.907 | 1.724 | 7.284 |
| 6000 | 3.88 | 3.861 | 12.09 | 3.92 | 3.55 | 13.8 |

These runs use one OpenMP thread per rank. The embedded provider is launched by
`--tmpi-threads`; native MPI uses the matching Open MPI launcher and an executable
created with `plumed mkmpi`. Both run the same kernel built with the same compiler
and options. Native launches use `--host localhost:N --bind-to none -n N`.
Eight ranks add substantial overhead on this eight-logical-CPU machine.

## Hybrid ranks and threads, buffered list

| Nucleotides | thread-MPI 2 × OMP2 | Open MPI 2 × OMP2 |
|---:|---:|---:|
| 1000 | 0.4801 | 0.4359 |
| 3000 | 1.33 | 1.324 |
| 6000 | 2.651 | 2.713 |

## Simulated domain decomposition, buffered list

| Nucleotides | thread-MPI 2 | thread-MPI 4 | Open MPI 2 | Open MPI 4 | thread-MPI 2 × OMP2 | Open MPI 2 × OMP2 |
|---:|---:|---:|---:|---:|---:|---:|
| 3000 | 1.919 | 1.9 | 1.965 | 1.723 | 1.318 | 1.34 |
| 6000 | 3.883 | 3.969 | 3.983 | 3.698 | 2.629 | 2.744 |

These extra trials use `benchmark --domain-decomposition`, which distributes
atom ownership and shuffles the mapping. They include atom sharing and force
application, in addition to the ERMSD calculation. Main-sweep MPI timings above
use the tool's default replicated-coordinate layout.

## Correctness

- Full values and atomic/box derivatives over 120 perturbed frames of 100 bases
  agree with the original kernel; selected `PAIRS` are covered.
- At 6000 bases, serial DIRECT/CELLS/VERLET agree exactly on the checked outputs;
  OpenMP and both MPI providers differ only by floating-point summation order.
  The maximum value difference over all differential checks is 6.99e-15;
  maximum derivative difference is 8.05e-15.
- Explicit benchmark comparisons with shuffled, distributed atom ownership also
  agree between serial, embedded thread-MPI4 and Open MPI4 (<8e-15).
- Analytical derivatives also pass PLUMED finite differences. The maximum error,
  including numerical box derivatives, is 8.95e-06. Atomic-coordinate errors
  in the six-base check are about 3.5e-7; coincident-center checks are about 2.3e-7.
  Numerical box tests use a finite 10 nm box to avoid a singular deformation.
- `rt-ermsd`, `rt-ermsd2`, `rt-ermsd-nopbc`, `rt-ermsd-neighbors`, and
  `rt-colvar-input-errors` pass. New cases cover contact entry/exit, displacement
  thresholds, rotating bases, reference-only contacts, an empty instantaneous
  contact map, exactly zero ERMSD, coincident centers and invalid options.
- The saved generator was rerun independently and reproduces the measured
  trajectory bytes and reference coordinates at 100 and 6000 bases.

## Exactness and implementation choices

The enclosing sphere has radius `CUTOFF/2`: the ellipsoid semiaxes are
`CUTOFF/2`, `CUTOFF/2`, and `0.3*CUTOFF` in the internal length representation.
For the usual cutoff of 2.4 these correspond to 1.2, 1.2 and 0.72 nm.
The default skin is one sixth of the sphere radius (0.2 nm); the list radius is
therefore 1.4 nm. Every evaluation checks center displacements and rebuilds when
any center has moved more than half the skin. No stride is supplied. The bound
remains valid under arbitrary base rotations.

Candidate rows include both buffered instantaneous neighbors and every nonzero
reference contact. This preserves contributions from contacts that disappear
and accumulates positive squared differences directly, avoiding cancellation
between large reference and current sums near zero ERMSD. Each ordered pair is
computed in its source base's frame, preserving the original normalization.

The existing `tools::NeighborList` linked-cell implementation used by
CoordinationBase is reused for both reference and current contacts. Only buffer
invalidation is ERMSD-specific. The dense `calcMat` toolbox interface remains
available for source callers; the CV does not allocate its quadratic arrays.
At fixed density, runtime storage instead follows the number of bases/contacts.

Pair work is distributed over the action communicator and OpenMP workers;
workers accumulate private atomic derivative buffers, followed by reductions.
Local frame construction is also threaded. The zero-distance mapping uses its
analytic limit, and exactly zero ERMSD avoids division by zero.

Batched periodic distance evaluation was not added: ERMSD already makes its
selected molecule whole (unless NOPBC), then deliberately ignores PBC. The
periodic image calculation that motivates CoordinationBase's batching is absent.
Linked cells remove the dominant unnecessary distance work directly.

## Input, timing method and scope

The source is [1S72](https://files.rcsb.org/download/1S72.pdb). The generator
extracts 2871 complete canonical RNA bases, using three ring atoms per base in
the proper purine/pyrimidine order. Sizes 100/300/1000 are nested prefixes;
3000/6000 add translated copies 25 nm apart. No protein or solvent atoms are
included. There are 120 periodic, smoothly displaced frames: 0.08 nm per-component
triplet motion plus 0.004 nm independent atom distortion, seed 1729. These are
synthetic structural perturbations, rather than an equilibrated MD trajectory.
All ERMSD benchmarks use NOPBC and a restraint so derivatives are applied.

The machine reports eight logical CPUs and 24 GiB RAM, arm64 macOS. Compilation
uses MacPorts clang 18.1.8, `-O3 -g -fopenmp`, default modules and embedded
thread-MPI. The native provider is MacPorts Open MPI 5.0.11. No fast-math flags
are used. The same kernel supports both MPI providers by runtime injection.

The reported statistic is `(B2 time + B3 time)/(B2 cycles + B3 cycles)` from
PLUMED's benchmark tool. Initialization, first step, warm-up and trajectory
parsing are excluded; per-step sharing, calculation and force application are
included. Normal trials run for three seconds. Slow serial all-pairs cases use
120 frames so even the largest calculation samples the complete trajectory.
The first 20% is warm-up. Other experiment builds/tests were not run concurrently
with timed trials. These results characterize this trajectory and host, through
6000 nucleotides; they are not production-MD throughput measurements.

[Detailed timing statistics](timings.csv) include minima/maxima and sample
counts. [Distributed-atom statistics](domain-timings.csv),
[validation results](validation.json), and [environment/kernel hashes](environment.json)
are also retained. Complete trial logs, exact commands, kernels, references and
trajectories remain in the project `benchmarks/` directory.
