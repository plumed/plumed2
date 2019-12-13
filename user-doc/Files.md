\page Files Files

We tried to design PLUMED in such a manner that input/output is done consistently
irrespective of the file type. Most of the files written or read by PLUMED thus follow
the very same conventions discussed below. 

\section Restart

Whenever the \ref RESTART option is used, all the files written by PLUMED are appended.
This makes it easy to analyze results of simulations performed as a chain of several sub-runs.
Notice that most of the PLUMED textual files have a header. The header is repeated at every
restart. Additionally, several files have time in the first column. PLUMED just takes the value
of the physical time from the MD engine. As such, you could have that time starts again from zero
upon restart or not.

An exception from this behavior is given by files which are not growing as the simulation proceeds.
For example, grids written with \ref METAD with GRID_WFILE are overwritten by default during the simulation.
As such, when restarting, there is no point in appending the file. Internally, PLUMED opens the file in append
mode but then rewinds it every time a new grid is dumped.

\section Backup

Whenever the \ref RESTART option is not used, PLUMED tries to write new files. If an old file
is found in the way, PLUMED takes a backup named "bck.X.filename" where X is a progressive number.
Notice that by default PLUMED only allows a maximum of 100 backup copies for a file.
This behavior can be changed by setting the environment variable PLUMED_MAXBACKUP to the desired number
of copies. E.g. export PLUMED_MAXBACKUP=10 will fail after 10 copies. PLUMED_MAXBACKUP=-1 will never fail - be careful
since your disk might fill up quickly with this setting.

\section Replica-Suffix Replica suffix

When running with multiple replicas (e.g., with GROMACS, -multi option) PLUMED adds the replica index as a suffix to
all the files. The following command
will thus print files named COLVAR.0, COLVAR.1, etc for the different replicas.
\plumedfile
d: DISTANCE ATOMS=1,2
PRINT ARG=d FILE=COLVAR
\endplumedfile

When reading a file, PLUMED will try to add the suffix. If the file is not found, it will fall back to
the name without suffix. The most important case is the reading of the plumed input file.
If you provide a file for each replica (e.g. plumed.0.dat, plumed.1.dat, etc) you will be able to
setup plumed differently on each replica. 
On the other hand, using a single plumed.dat will make all the replicas read the same file.

\warning This rule is true for almost all the files read by PLUMED. As of
  PLUMED version 2.4, the only exception is PDB files, where the replica suffix is not added.

Notice that when PLUMED adds the replica suffix, it recognizes the file extension and add the suffix _before_ the
extension. Before PLUMED 2.2, the only recognized suffix was ".gz". Since 2.2, any suffix with length
less or equal to five letters is recognized.

This means that using in a multi-replica context an input such as
\plumedfile
d: DISTANCE ATOMS=1,2
PRINT ARG=d FILE=COLVAR.gz
METAD ARG=d FILE=test.HILLS SIGMA=0.1 HEIGHT=0.1 PACE=100
\endplumedfile
PLUMED will write files named COLVAR.0.gz, COLVAR.1.gz, test.0.HILLS, test.1.HILLS, etc
etc. This is useful since the preserved extension makes it easy
to process the files later.

