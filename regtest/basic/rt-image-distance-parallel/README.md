IMAGE_DISTANCE OpenMP/MPI validation and benchmark
=================================================

The MPI harness compares SERIAL and default modes bitwise for values, all
atomic forces and all nine virial components. It covers ranks with no pairs,
one-atom box-only behavior, extended coordinates, a skew shortest lattice
vector, 513-atom triclinic groups, coherent translations and exact ties whose
candidate pairs span multiple ranks/threads. An invalid pair owned only by rank
zero checks that errors are propagated collectively before any rank throws.
The original oracle/finite-difference fixture remains independent validation of
the geometry; this fixture specifically tests parallel partition/reduction.

The default test uses four ranks and the configured PLUMED_NUM_THREADS. Repeat
with 1, 2 and 4 OpenMP threads in separate processes to test team-size invariance.
The harness also supports `-n N` for internal thread-MPI and native MPI; use the
matching native Open MPI launcher for native builds. All tests print one shared
output file from rank zero; each rank has its own diagnostic log.

The optional benchmark uses fixed seeds and identical compact-group coordinates
in a 4x5x6 cell, three warmup frames, then three independently timed batches.
All modes must reproduce the serial reference distance within 1e-12.
Reported CSV fields are mode,ranks,threads,atoms,frames,median_ms_per_frame,distance.
Threads is the configured count; the action uses one thread when local work
is below 4096 pairs per requested thread.
Each batch uses the maximum elapsed time over ranks, excluding initialization.
Sizes are 100/1000/5000 with 1000/100/8 frames per batch. Timing includes the
public MD interface, force/virial propagation and MPI winner communication.

From a configured/built checkout, compile/run this fixture first. Its tmp/exe
accepts `--benchmark parallel` (default action), `--benchmark serial` (SERIAL
keyword), or `--benchmark reference` (the pre-parallelization action compiled
into a dedicated executable by benchmark.py). benchmark.py builds that reference
from an explicitly supplied Git revision, then runs a scaling matrix.
The reference differs only in its class/action name and public header includes;
it uses the same current kernel, compiler flags, input and timing harness.
The original serial implementation was preserved by rebased commit
03d741ee2d916ac2c52be6e83a09e989c057366b.

Example using the project-local paths documented in the project notes:

```
python3 benchmark.py --exe /absolute/path/to/rt-image-distance-parallel/tmp/exe \
  --reference-ref 03d741ee2d916ac2c52be6e83a09e989c057366b \
  --work /absolute/path/to/project/.tmp/image-distance-benchmark \
  --launcher /opt/local/bin/mpiexec-openmpi-mp --openmpi-local
```

Compiler wrapper overrides, if needed, must match those used to build the
kernel (e.g. OMPI_CXX/OMPI_CC); the script inherits them. It creates only files
inside --work. Timings are measurements, not regression thresholds. Large groups
remain quadratic in total pair work; MPI/thread overhead can outweigh gains for
small selections. No image or pair matrix is materialized.
