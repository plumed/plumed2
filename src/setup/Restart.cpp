/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "core/ActionSetup.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Checkpointing.h"
#include "tools/Exception.h"

namespace PLMD {
namespace setup {

//+PLUMEDOC GENERIC RESTART
/*
Activate restart.

Many MD calculations are run on supercomputer facilities that limit the amount of
time an individual calculation can run for.  As MD simulations take a long time to
complete, you often have to restart your calculations and build up your trajectories by
running a chain of sub-runs. Some MD codes tell PLUMED if the MD has restarted
the calculation from a checkpoint file while others do not. To avoid problems we recommend
using the RESTART commands and keywords in your input files whenever you are generating
trajetories by restarting calculations from a checkpoint file.

RESTART is a Setup directive and, as such, should appear
at the beginning of the input file. When you include the RESTART command
in your plumed.dat file you are telling PLUMED that you want to restart all
the actions. Restarting these actions ensures that, as discussed [here](Files.md), PLUMED appends
all the files that are open for writing. Appending to these files makes it
easier to analyze results from simulations that have been performed as a
chain of several sub-runs.

To consider what this means in practise consider the following input:

```plumed
d: DISTANCE ATOMS=1,2
PRINT ARG=d FILE=out
```

If you run a plumed calculation in a directory which already contains a file called
`out` that file called `out` will be renamed `bck.0.out` so that the data from the new
calculation can be output to a new version of the file `out` and so that the data in the
old version of the file is not lost.

If by contrast you run the calculation in that same directory with the following input:

```plumed
RESTART
d: DISTANCE ATOMS=1,2
PRINT ARG=d FILE=out
```

The data from the new simulation simulation will be appended to the file called `out` that
was already present in the data.

Notice that it is also possible to enable or disable restart on a per-action
basis using the RESTART keyword as illustrated below:

```plumed
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=1,2
PRINT ARG=d1 FILE=out1
PRINT ARG=d2 FILE=out2 RESTART=YES
```

If this input is used any files called `out1` are backed up to `bck.0.out1`.  Data from the
new calculation is then appended to any file called `out2`  that is present in the working directory.

Notice that the RESTART keyword is assigned one of three values:

- `RESTART=AUTO` is the default and means that global settings are used
- `RESTART=YES` means that the action is restarted
- `RESTART=NO` means that the action is not restarted even if there is a global restart command

The following input is thus equivalent to the first input that introduced the RESTART keyword.

```plumed
RESTART
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=1,2
PRINT ARG=d1 FILE=out1 RESTART=NO
PRINT ARG=d2 FILE=out2
```

[!CAUTION]
The RESTART directive can have also other important side effects. For methods such as [METAD](METAD.md)
or [PBMETAD](PBMETAD.md) when the action is told to restart all the hills from the previous run are read
in at the start of the calculation so that the bias potential at the start of the restarted simulation is
identical to the bias potential at the end of the previous run.

If you are absolutely certain that you do not want to restart any of the actions in your input you
can use an input like the example below:

```plumed
RESTART NO
d1: DISTANCE ATOMS=1,2
PRINT ARG=d1 FILE=out1
```

The `RESTART NO` command in this input tells PLUMED that no actions are to be restarted even if the MD
code has told PLUMED that the initial configuration came from a checkpoint file and that the trajectory
is being restarted.

## Trimming the files at the restart step

When an MD code restarts from its checkpoint, the simulation continues from the step S at which the
checkpoint was written, while the files written by PLUMED may contain lines written after it, by a run
that was interrupted (or that went on) after writing its checkpoint. If the MD code tells PLUMED the step
S from which it restarts (the GROMACS patches do this after `-cpi`, and [driver](driver.md) does it with
`--restart --initial-step S`), PLUMED removes these lines from the files it appends to, before reading or
writing them, so that they continue exactly where the restarted simulation starts, as GROMACS does with its
own files. The lines of step S are removed too, since step S is computed again. In this way the COLVAR
files contain every time once and, for instance, [METAD](METAD.md) reads only the hills deposited before
step S (the hill of step S is added again) and [OPES_METAD](OPES_METAD.md) only the kernels deposited before it.

The lines removed are the last lines of the file with a time equal to or later than the time of step S:
the earlier lines are kept, including those of earlier runs whose time restarted from zero. The trimming
applies to the files that have a `time` field (PRINT, the hills of METAD, the kernels of OPES_METAD, ...);
the other files (e.g. those of [DUMPATOMS](DUMPATOMS.md)) and the files where the times after S are not
increasing are appended as they are, with a line in the log. The lines removed from a file are not lost:
they are moved to a file `bck.trimmed.N.name` in the same directory (N is the first number not used yet).
Compressed files (`.gz`, written in append mode as a sequence of gzip streams) are trimmed too: since they
cannot be truncated, the lines that are kept are compressed again. Files that are rewritten at
every write (e.g. grids and state files) are not trimmed. When the bias is read from a grid
([METAD](METAD.md) with GRID_RFILE) or from a state file ([OPES_METAD](OPES_METAD.md) with STATE_RFILE),
which are written at the checkpoint of the MD code and already contain the hill or kernel of step S,
the lines of step S are kept. With the multiple walkers of METAD that share files (WALKERS_N) the hills
files are not trimmed. The log lists the files that have been trimmed and where their removed lines are. When a trimmed file ends with the same field header that the restarted run prints, the header is not printed again, so that the file is the same as the one of an uninterrupted run.

The trimming can be switched off for all the files, which are then appended as they are, with

```plumed
RESTART NO_TRIM
d1: DISTANCE ATOMS=1,2
PRINT ARG=d1 FILE=out1
```

or for the files of a single action with the NO_TRIM flag:

```plumed
d1: DISTANCE ATOMS=1,2
PRINT ARG=d1 FILE=out1 NO_TRIM
```

Notice that `RESTART NO_TRIM`, as RESTART, also switches restart on. With [CHECKPOINT](CHECKPOINT.md) the files
are instead truncated to the sizes they had when the checkpoint was written, and the restart is exact.

*/
//+ENDPLUMEDOC

class Restart :
  public virtual ActionSetup {
public:
  static void registerKeywords( Keywords& keys );
  explicit Restart(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(Restart,"RESTART")

void Restart::registerKeywords( Keywords& keys ) {
  ActionSetup::registerKeywords(keys);
  keys.addFlag("NO",false,"switch off restart - can be used to override the behavior of the MD engine");
  // the flag reserved for all the actions, with the meaning it has here (read by Action, see getNoTrim)
  keys.remove("NO_TRIM");
  keys.addFlag("NO_TRIM",false,"append to the files as they are, without removing the lines written at or after the step from which the MD engine restarts");
}

Restart::Restart(const ActionOptions&ao):
  Action(ao),
  ActionSetup(ao) {
  bool no=false;
  parseFlag("NO",no);
  if(getNoTrim()) {
    plumed.checkpointing.disableTrimming();
    log<<"  the files are not trimmed at the restart step\n";
  }
  bool md=plumed.getRestart();
  log<<"  MD code "<<(md?"did":"didn't")<<" require restart\n";
  if(no) {
    if(plumed.checkpointing.isRestartingFromCheckpoint()) {
      error("RESTART NO cannot be used when restarting from a checkpoint (see CHECKPOINT)");
    }
    if(md) {
      log<<"  Switching off restart\n";
    }
    plumed.setRestart(false);
    log<<"  Not restarting simulation: files will be backed up\n";
  } else {
    if(!md) {
      log<<"  Switching on restart\n";
    }
    plumed.setRestart(true);
    log<<"  Restarting simulation: files will be appended\n";
  }
}

}
}
