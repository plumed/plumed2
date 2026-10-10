/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "core/ActionRegister.h"
#include "core/ActionPilot.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"

namespace PLMD {
namespace generic {

//+PLUMEDOC GENERIC DEBUG
/*
Set some debug options.

This action can be used while debugging or optimizing plumed. The example
input below demonstates two useful ways you can use the DEBUG command

```plumed
# print detailed (action-by-action) timers at the end of simulation
a: DEBUG DETAILED_TIMERS
# dump every two steps which are the atoms required from the MD code
b: DEBUG logRequestedAtoms STRIDE=2
```


For diagnosing non-finite values in an action graph, set the environment
variable `PLUMED_NONFINITE_ACTION_TRACE` before starting the program. It
accepts comma- or space-separated action labels, or `*` for all labels.
The selection is read once per process and is disabled when unset. This
option is independent of the DEBUG action.

Selected actions are checked after their forward calculation for non-finite
values and scalar derivatives. During force propagation, selected action
forces and the atomic-position and cell force buffers are also checked.
A failure reports the step, forward/backward phase, triggering action,
target component and element. The diagnostic can add substantial overhead.
All spatial MPI ranks must use the same selection and action graph. The first failing rank is reported to every rank
before an exception is raised. Separate multi-simulation communicators are
not synchronized by this diagnostic; the MD code must terminate the other
replicas if one replica raises an exception.
It does not check all internal action state or grid-derivative layouts, and
cannot detect errors introduced later by the molecular dynamics engine.
A finite result is not evidence of a valid physical simulation.

*/
//+ENDPLUMEDOC
class Debug:
  public ActionPilot {
  OFile ofile;
  bool logActivity;
  bool logRequestedAtoms;
  bool novirial;
  bool detailedTimers;
public:
  explicit Debug(const ActionOptions&ao);
/// Register all the relevant keywords for the action
  static void registerKeywords( Keywords& keys );
  void calculate() override {}
  void apply() override;
};

PLUMED_REGISTER_ACTION(Debug,"DEBUG")

void Debug::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords(keys);
  keys.add("compulsory","STRIDE","1","the frequency with which this action is to be performed");
  keys.addFlag("LOGACTIVITY",false,"write in the log which actions are inactive and which are inactive,"
               " in PLUMED<2.6 and this flag can be activated only as 'logActivity'");
  keys.addFlag("LOGREQUESTEDATOMS",false,"write in the log which atoms have been requested at a given time,"
               " in PLUMED<2.6 and this flag can be activated only as 'logRequestedAtoms'");
  keys.addFlag("NOVIRIAL",false,"switch off the virial contribution for the entirety of the simulation");
  keys.addFlag("DETAILED_TIMERS",false,"switch on detailed timers");
  keys.add("optional","FILE","the name of the file on which to output these quantities");
}

Debug::Debug(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  logActivity(false),
  logRequestedAtoms(false),
  novirial(false) {
  parseFlag("LOGACTIVITY",logActivity);
  if(logActivity) {
    log.printf("  logging activity\n");
  }
  parseFlag("LOGREQUESTEDATOMS",logRequestedAtoms);
  if(logRequestedAtoms) {
    log.printf("  logging requested atoms\n");
  }
  parseFlag("NOVIRIAL",novirial);
  if(novirial) {
    log.printf("  Switching off virial contribution\n");
  }
  if(novirial) {
    plumed.novirial=true;
  }
  parseFlag("DETAILED_TIMERS",detailedTimers);
  if(detailedTimers) {
    log.printf("  Detailed timing on\n");
    plumed.detailedTimers=true;
  }
  ofile.link(*this);
  std::string file;
  parse("FILE",file);
  if(file.length()>0) {
    ofile.open(file);
    log.printf("  on file %s\n",file.c_str());
  } else {
    log.printf("  on plumed log file\n");
    ofile.link(log);
  }
  checkRead();
}

void Debug::apply() {
  if(logActivity) {
    const ActionSet&actionSet(plumed.getActionSet());
    int a=0;
    for(const auto & p : actionSet) {
      if(dynamic_cast<Debug*>(p.get())) {
        continue;
      }
      if(p->isActive()) {
        a++;
      }
    };
    if(a>0) {
      ofile<<"activity at step "<<getStep()<<": ";
      for(const auto & p : actionSet) {
        if(dynamic_cast<Debug*>(p.get())) {
          continue;
        }
        if(p->isActive()) {
          ofile.printf("+");
        } else {
          ofile.printf("-");
        }
      };
      ofile.printf("\n");
    };
  };
  if(logRequestedAtoms) {
    ofile<<"requested atoms at step "<<getStep()<<": ";
    const int* l;
    int n;
    plumed.cmd("createFullList",&n);
    plumed.cmd("getFullList",&l);
    for(int i=0; i<n; i++) {
      ofile.printf(" %d",l[i]);
    }
    ofile.printf("\n");
    plumed.cmd("clearFullList");
  }

}

}
}

