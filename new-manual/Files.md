# Files

We tried to design PLUMED in such a manner that input/output is done consistently
irrespective of the file type. Most of the files written or read by PLUMED thus follow
the conventions discussed below. 

## Restart

Whenever the [RESTART](RESTART.md) option is used, all the files written by PLUMED are appended.
This makes it easy to analyze results of simulations performed as a chain of several sub-runs.
Notice that most of the PLUMED textual files have a header. The header is repeated at every
restart, unless the file was trimmed (see below) and ends with the same header. Additionally, several files have time in the first column. PLUMED just takes the value
of the physical time from the MD engine. As such, you could have that time starts again from zero
upon restart or not.

An exception from this behavior is given by files which are not growing as the simulation proceeds.
For example, grids written with [METAD](METAD.md) with `GRID_WFILE` are overwritten by default during the simulation.
As such, when restarting, there is no point in appending the file. Internally, PLUMED opens the file in append
mode but then rewinds it every time a new grid is dumped.

## Trimming the files at the restart step

When the MD engine tells PLUMED the step from which it restarts (the GROMACS patches do it after `mdrun -cpi`),
PLUMED removes from the files it appends to the lines written at or after that step by the interrupted run,
so that the files continue exactly where the restarted simulation starts, without steps written twice. The
removed lines are moved to `bck.trimmed.X.filename`. Only the last lines with a time not earlier than the one
of the restart step are removed, so the earlier sub-runs whose time restarted from zero are kept. Files
without a time field are appended as they are. The trimming can be switched off with `RESTART NO_TRIM`, or
for the files of a single action with the `NO_TRIM` flag. See [RESTART](RESTART.md) for the details.

## Exact restart with a checkpoint

Appending (also after trimming) leaves in the files the lines that the interrupted run wrote after the checkpoint
of the MD engine, so that after the restart the files contain the same steps twice, and
the actions that keep an internal state (e.g. the hills of [METAD](METAD.md)) are not
restarted exactly. With the [CHECKPOINT](CHECKPOINT.md) directive PLUMED writes its own
checkpoint file whenever the MD engine writes its own: at restart every file is truncated
to the size it had at the checkpoint, the state of the actions is restored and the checkpointing
step is computed again, so that the restarted run continues the interrupted one bitwise,
without repeating the header of the files.

## Backup

Whenever the [RESTART](RESTART.md) option is not used, PLUMED tries to write new files. If an old file
is found in the way, PLUMED takes a backup named "bck.X.filename" where X is a progressive number.
Notice that by default PLUMED only allows a maximum of 100 backup copies for a file.
This behavior can be changed by setting the environment variable `PLUMED_MAXBACKUP` to the desired number
of copies. E.g. export `PLUMED_MAXBACKUP=10` will fail after 10 copies. `PLUMED_MAXBACKUP=-1` will never fail - be careful
since your disk might fill up quickly with this setting.

## HDF5 tables

The files that PLUMED writes field by field, i.e. those with a `#! FIELDS` header such as the ones of
[PRINT](PRINT.md), the hills of [METAD](METAD.md) or the kernels of [OPES_METAD](OPES_METAD.md), are written in
[HDF5](https://www.hdfgroup.org/solutions/hdf5/) format when their name ends with `.h5`:

```plumed
d: DISTANCE ATOMS=1,2
PRINT ARG=d FILE=COLVAR.h5
METAD ARG=d SIGMA=0.1 HEIGHT=0.1 PACE=100 FILE=HILLS.h5
```

The values are stored with full precision (FMT is ignored), so that a restart from an HDF5 file of hills or
kernels gives the same bias as the run that wrote it. The files can be read with any HDF5 tool, e.g. in Python:

```python
import h5py

with h5py.File("COLVAR.h5") as f:
    time = f["block_0/time"][:]
    d = f["block_0/d"][:]
```

Every field is a dataset of the group `block_0`. A new group (`block_1`, ...) starts when the fields change (where
a text file would have a new `#! FIELDS` line), e.g. at a restart with a different input; the constant fields
(`#! SET` in text files) are attributes of the group, and the attributes `plumed_columns` and `plumed_constants` list
the fields in order. PLUMED reads the HDF5 tables it writes wherever it reads text files (restarts, [READ](READ.md),
...), and the [h5convert](h5convert.md) tool converts them to text files and back, so that every tool that reads
text files can be used. Restarts (appending, trimming at the restart step, [CHECKPOINT](CHECKPOINT.md)) work as with
text files: the removed rows are moved to an HDF5 file `bck.trimmed.X.filename.h5`.

HDF5 tables require PLUMED to be compiled with libhdf5 (`./configure` uses it when it finds it, see
`--enable-libhdf5`; the conda package of PLUMED uses the hdf5 package of conda). Only the files written field by
field can be HDF5 tables: the files written as free text (e.g. those of [DUMPATOMS](DUMPATOMS.md) and
[PRINT_NDX](PRINT_NDX.md)) keep their text formats. The hills files of multiple walkers that share files
(WALKERS_N of [METAD](METAD.md)) are text files too, since every walker reads the files of the others while they are
written (with WALKERS_MPI the hills can be an HDF5 table).

Grids (e.g. GRID_WFILE of [METAD](METAD.md), [DUMPGRID](DUMPGRID.md)) and state files (STATE_WFILE of
[OPES_METAD](OPES_METAD.md)) can be HDF5 tables too, with a row for every point of the grid as in the text files.
With their full precision, a restart from a grid (GRID_RFILE) or from a state (STATE_RFILE) gives the same bias as
the run that wrote them. When they are rewritten, the previous one is kept as `bck.last.filename`, as for text files.
For a COLVAR file with a few columns an HDF5 table is
not smaller or faster than a text file: HDF5 is useful for large files and when the full precision is needed.

# Replica-Suffix Replica suffix

When running with multiple replicas (e.g., with GROMACS, -multi option) PLUMED adds the replica index as a suffix to
all the files. The following command
will thus print files named COLVAR.0, COLVAR.1, etc for the different replicas.

```plumed
#SETTINGS NREPLICAS=2
d: DISTANCE ATOMS=1,2
PRINT ARG=d FILE=COLVAR
```

When reading a file, PLUMED will try to add the suffix. If the file is not found, it will fall back to
the name without suffix. The most important case is the reading of the plumed input file.
If you provide a file for each replica (e.g. plumed.0.dat, plumed.1.dat, etc) you will be able to
setup plumed differently on each replica. 
On the other hand, using a single plumed.dat will make all the replicas read the same file.

>[!CAUTION]
>This rule is true for almost all the files read by PLUMED. As of
>PLUMED version 2.4, the only exception is PDB files, where the replica suffix is not added.

Notice that when PLUMED adds the replica suffix, it recognizes the file extension and add the suffix _before_ the
extension. Before PLUMED 2.2, the only recognized suffix was ".gz". Since 2.2, any suffix with length
less or equal to five letters is recognized.

This means that if you use the following input in a multi-replica context 

```plumed
#SETTINGS NREPLICAS=2
d: DISTANCE ATOMS=1,2
PRINT ARG=d FILE=COLVAR.gz
METAD ARG=d FILE=test.HILLS SIGMA=0.1 HEIGHT=0.1 PACE=100
```

PLUMED will write files named COLVAR.0.gz, COLVAR.1.gz, test.0.HILLS, test.1.HILLS, etc
etc. This is useful since the preserved extension makes it easy
to process the files later.

