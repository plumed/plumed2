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
}

Restart::Restart(const ActionOptions&ao):
  Action(ao),
  ActionSetup(ao) {
  bool no=false;
  parseFlag("NO",no);
  bool md=plumed.getRestart();
  log<<"  MD code "<<(md?"did":"didn't")<<" require restart\n";
  if(no) {
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
