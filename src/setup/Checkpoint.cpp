/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2026 The plumed team
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
#include <filesystem>

namespace PLMD {
namespace setup {

//+PLUMEDOC GENERIC CHECKPOINT
/*
Write a PLUMED checkpoint file for an exact restart.

With this directive PLUMED writes its own checkpoint file, separate from the
one of the MD engine, whenever the MD engine writes its checkpoint (the MD
engine tells PLUMED through the `doCheckPoint` command of the interface; the
GROMACS patch does so) or, with STRIDE, every STRIDE steps for the engines
that do not send this signal.

The checkpoint stores the internal state of the actions (hills, accumulators,
random number generators, ...) at the start of the checkpointing step, together
with the sizes of all the files that PLUMED has written up to that moment.
When the simulation is restarted from it, every output file is truncated to the
size it had at the checkpoint, so that the lines written after it by the
interrupted run are removed, and the checkpointing step is computed again:
this requires the MD engine to restart from its own checkpoint of the same
step, as GROMACS does with `-cpi`. The result is an exact continuation of the
interrupted run, at variance with a plain [RESTART](RESTART.md) where files
are only appended.

CHECKPOINT is a setup directive and must come before any action that writes a
file, so put it at the top of the input. The file is written atomically and the
previous checkpoint is kept with the suffix `_prev` (`plumed_prev.cpt` for
`plumed.cpt`); with multiple replicas every replica writes its own file with the
usual replica suffix.

```plumed
CHECKPOINT FILE=plumed.cpt
d: DISTANCE ATOMS=1,2
PRINT ARG=d FILE=colvar
```

At restart, RESTART_FROM selects the checkpoint to read. The default, AUTO,
reads FILE when the MD engine tells PLUMED that the run is restarted (or when
[RESTART](RESTART.md) comes before CHECKPOINT) and starts a new run otherwise.
A file name reads that checkpoint unconditionally, NO never reads one (the
checkpoint is only written). Restarting from a checkpoint implies RESTART.

```plumed
CHECKPOINT FILE=plumed.cpt STRIDE=1000 RESTART_FROM=AUTO
d: DISTANCE ATOMS=1,2
PRINT ARG=d FILE=colvar
```

The first step presented by the MD engine after the restart must be the one
stored in the checkpoint, otherwise PLUMED stops with an error: this happens
when the PLUMED checkpoint and the restart file of the MD engine do not come
from the same moment of the simulation. With STRIDE, choose a stride that
matches the checkpointing of the MD engine.

An action with the per-action keyword `RESTART=NO` starts from scratch at the
restart step: its state is not restored from the checkpoint and its files,
trimmed at the checkpoint like all the others, are backed up and rewritten
instead of being appended. This is the way to add a new action with an
internal state to a continuation, or to reset one (e.g. a [METAD](METAD.md)
that starts depositing again from zero while everything else continues
exactly). The global `RESTART NO` directive after CHECKPOINT is an error,
since it contradicts the restart from the checkpoint.

The actions that store their state in the checkpoint include [METAD](METAD.md),
[PBMETAD](PBMETAD.md), [OPES_METAD](OPES_METAD.md), [OPES_METAD_EXPLORE](OPES_METAD_EXPLORE.md),
[METAINFERENCE](METAINFERENCE.md) (and the ISDB collective variables with DOSCORE),
[EDS](EDS.md), [DRR](DRR.md), [MOVINGRESTRAINT](MOVINGRESTRAINT.md), [ABMD](ABMD.md),
[EXTENDED_LAGRANGIAN](EXTENDED_LAGRANGIAN.md), [ACCUMULATE](ACCUMULATE.md), [COLLECT](COLLECT.md),
[COMMITTOR](COMMITTOR.md), [MAXENT](MAXENT.md), [READ](READ.md), [RANDOM_EXCHANGES](RANDOM_EXCHANGES.md),
[LOGMFD](LOGMFD.md), [FISST](FISST.md), [RESCALE](RESCALE.md), [EMMI](EMMI.md), [EMMIVOX](EMMIVOX.md) and
[EFFECTIVE_ENERGY_DRIFT](EFFECTIVE_ENERGY_DRIFT.md) (only with a single process): they
restore their state from it, so that their hills, kernels or status files are not read
(the growing ones are truncated).
Actions that keep an internal state but do not store it in the checkpoint
yet are listed in a warning in the log: for them the restart is the usual
file-based one and an exact continuation is not guaranteed.

The content of a checkpoint can be inspected with `plumed checkpoint --dump`.
Compressed (`.gz`) output files cannot be truncated and are not allowed
together with CHECKPOINT.

*/
//+ENDPLUMEDOC

class Checkpoint :
  public virtual ActionSetup {
public:
  static void registerKeywords( Keywords& keys );
  explicit Checkpoint(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(Checkpoint,"CHECKPOINT")

void Checkpoint::registerKeywords( Keywords& keys ) {
  ActionSetup::registerKeywords(keys);
  keys.add("compulsory","FILE","plumed.cpt","the name of the checkpoint file");
  keys.add("compulsory","STRIDE","0","write a checkpoint every STRIDE steps, in addition to those requested by the MD engine. 0 means only on request of the MD engine");
  keys.add("compulsory","RESTART_FROM","AUTO","the checkpoint to restart from: AUTO (FILE, if the MD engine is restarting), NO (never restart from a checkpoint), or the name of a checkpoint file");
}

Checkpoint::Checkpoint(const ActionOptions&ao):
  Action(ao),
  ActionSetup(ao) {
  std::string file;
  parse("FILE",file);
  long long int stride=0;
  parse("STRIDE",stride);
  if(stride<0) {
    error("STRIDE should be positive");
  }
  std::string restartFrom;
  parse("RESTART_FROM",restartFrom);
  checkRead();
  Checkpointing& cpt=plumed.checkpointing;
  if(restartFrom=="AUTO") {
    if(plumed.getRestart()) {
      restartFrom=file;
      log<<"  the MD engine is restarting: reading the checkpoint\n";
    } else {
      restartFrom="NO";
      log<<"  new run: not reading any checkpoint\n";
    }
  }
  if(restartFrom!="NO") {
    const std::string fname=cpt.suffixedPath(restartFrom);
    if(!std::filesystem::exists(fname)) {
      error("checkpoint file "+fname+" not found. Use RESTART_FROM=NO to restart without a checkpoint");
    }
    cpt.restartFrom(restartFrom);
  }
  cpt.setup(file,stride);
  log<<"  checkpoint file: "<<cpt.getPath()<<"\n";
  if(stride>0) {
    log<<"  written every "<<stride<<" steps and when the MD engine writes its checkpoint\n";
  } else {
    log<<"  written when the MD engine writes its checkpoint\n";
  }
}

}
}
