This module provides [REPLEX_BIASVALUE](REPLEX_BIASVALUE.md), a replica exchange
implementation for multi-simulation runs in which the replicas exchange their
biases instead of their coordinates. The biases are linear combinations of
PLUMED values, and can be arranged on a multidimensional ladder. The trajectories
remain continuous in time, and the bias indices written by the action allow one
to recover the ensemble of each bias.

The module was originally developed for the paper cited below.

## Installation

This module is compiled by default. It can be disabled with
`--enable-modules=-replex_biasvalue` in the `./configure` command used to build
PLUMED.

## Usage

The exchanges are done entirely by PLUMED, so the MD code only needs to run
PLUMED on multiple replicas (e.g. `mdrun -multidir` in GROMACS, or
`plumed driver --multi`), without its own replica exchange. Since only the
biases differ between replicas, this is a form of Hamiltonian replica exchange
that requires no change to the MD code, and that also allows ladders in more
than one dimension.

A typical input computes the CVs and the bias functions as usual, and passes them
to a single REPLEX_BIASVALUE action, identical on all replicas, which adds the
bias of the current replica. The bias index of each replica, written to a file,
is then used to sort the trajectories by bias, e.g. with `gmx trjcat -demux`.
See [REPLEX_BIASVALUE](REPLEX_BIASVALUE.md) for examples, the output files and
restarts.
