# Hybrid OPES checkpoint restart

This repeats `rt-checkpoint-opes-mpi` with two requested OpenMP threads per
rank. Both tests use four MPI ranks, arranged as two replicas with two ranks
each, and exercise MPI walkers and adaptive kernels. Without OpenMP support,
PLUMED uses its usual single-thread fallback.

The original test forces one OpenMP thread and requires byte-identical COLVAR
and kernel files after restart. This test allows absolute/relative differences
of `1e-12` in COLVAR values because the order of OpenMP reductions is not fixed.
It still requires identical headers, row/column counts, finite values, and
byte-identical kernel files. Both tests also verify checkpoint restoration and
output truncation messages in separate replica logs.

The comparator is local to this test and uses only the Python standard library.
