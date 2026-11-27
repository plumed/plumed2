IMAGE_DISTANCE spatial / OpenMP / MPI validation and benchmark
============================================================

The MPI harness compares SERIAL BRUTE_FORCE (the exhaustive pair/image scan),
SERIAL (one-worker spatial search) and default spatial work sharing bitwise for
values, every atomic force and all nine virial components. It covers empty
ranks, one-atom self images, extended coordinates, skew lattices, 513-atom
triclinic groups, exact distributed ties, coherent translations, permutations,
12 random spatial-bound cases, changing cells/coordinates on reused actions,
large translations, many-image fallback and collective worker errors.
The original rt-image-distance fixture independently checks the geometry with
wide image enumeration and atomic/box finite differences, including the tree
path on a 65-atom group and a replicated-group independent oracle.

The default fixture uses four ranks and the configured PLUMED_NUM_THREADS.
Repeat with 1/2/4 ranks and 1/2/4 threads in separate processes. The executable
accepts `-n N` for internal thread-MPI and native MPI; native builds need the
matching launcher too. Rank zero writes output; diagnostic logs are per rank.

The tree is rebuilt each frame. BRUTE_FORCE disables spatial pruning; SERIAL
only disables MPI/OpenMP work sharing. The tree falls back to the exhaustive
scan for fewer than 64 atoms, more than 4096 globally bounded shifts, or
more than max(32,N/8) shifts (a conservative query-cost guard).
There is no skin, cutoff approximation, stale list or cached previous winner.

The optional benchmark uses fixed seeds and identical inputs for the old Git
reference action and the current action in the SAME executable/compiler/kernel.
Three warmup frames precede three timed batches; the median maximum-rank elapsed
time is reported, excluding action initialization and including tree rebuilding,
the public MD interface, atomic/box force propagation and MPI communication.
Sizes are 100/1000/5000 with 1000/100/8 frames per batch. CSV columns are
mode,ranks,threads,atoms,frames,median_ms,distance,profile. Threads is configured;
the action caps the team by available rows and avoids threads below 4096 local
potential pairs per requested thread. All returned distances must agree within
1e-12 with the reference for that profile/size. Timings have no test thresholds.

Profiles (coordinates already whole):
- compact: uniform [-0.7,0.7]^3, 4x5x6 orthorhombic cell (previous benchmark).
- triclinic: same cloud, cell rows (4,0,0), (1.3,3.7,0), (-0.8,0.7,5.1).
- extended: compact cloud stretched eightfold along x, 4x5x6 cell.
- wide: compact cloud stretched 1000-fold along x; exercises the query-cost
  fallback rather than enumerating many globally overlapping images.
- density: cubic cloud at 100 atoms/nm^3 with a cubic cell 1 nm larger per side;
  molecular extent and cell grow with atom count. This is a synthetic spatial
  distribution, not an equilibrated molecular structure.
- coincident: all atoms at the origin, 4x5x6 cell; a degenerate stress case in
  which many tied candidates survive. It is optional and can remain quadratic.

Compile/run this fixture in a configured checkout first. benchmark.py reuses
its tmp/Makefile to build a separate executable containing the reference source
from --reference-ref, renamed to IMAGE_DISTANCE_REFERENCE. Only class/action
names and public header includes are changed. A separately linked plugin would
register in a different kernel registry for static test executables.

The spatial comparison uses the previous parallel scan commit
`ef3a0980a4ce17f174cfe0d35e8950cb05107842`. Each profile runs matched reference /
current cases at 1x1, 1x4, 1x8, 4x1 and 4x2 ranks x configured threads.
The older serial reference remains available at 03d741ee2, but it cannot measure
MPI/OpenMP scaling; choose the reference revision explicitly.

```
python3 benchmark.py --exe /absolute/path/to/rt-image-distance-parallel/tmp/exe \
  --reference-ref ef3a0980a4ce17f174cfe0d35e8950cb05107842 \
  --work /absolute/path/to/project/.tmp/image-distance-spatial-benchmark \
  --launcher /opt/local/bin/mpiexec-openmpi-mp --openmpi-local \
  --profiles compact triclinic extended density wide coincident
```

Compiler wrapper overrides, if needed, must match the kernel (e.g. OMPI_CXX /
OMPI_CC); the script inherits them. All scratch files go under --work. Compact/triclinic/extended/density
profiles are the default; wide and coincident must be requested explicitly. Direct
`tmp/exe -n N --benchmark parallel PROFILE` or `--benchmark serial PROFILE`
measure the current action; reference mode requires the executable built by
the script. Dense/tied cases can retain quadratic work. Tree construction is
replicated on MPI ranks and can dominate the much shorter accelerated query.
