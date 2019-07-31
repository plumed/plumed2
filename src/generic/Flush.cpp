/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2019 The plumed team
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
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

namespace PLMD {
namespace generic {

using namespace std;

//+PLUMEDOC GENERIC FLUSH
/*
This command instructs plumed to flush all the open files with a user specified frequency.
Notice that all files are flushed anyway every 10000 steps.

This
is useful for preventing data loss that would otherwise arise as a consequence of the code
storing data for printing in the buffers. Notice that wherever it is written in the
plumed input file, it will flush all the open files.

\par Examples

A command like this in the input will instruct plumed to flush all the output files every 100 steps
\plumedfile
d1: DISTANCE ATOMS=1,10
PRINT ARG=d1 STRIDE=5 FILE=colvar1

FLUSH STRIDE=100

d2: DISTANCE ATOMS=2,11
# also this print is flushed every 100 steps:
PRINT ARG=d2 STRIDE=10 FILE=colvar2
\endplumedfile
(see also \ref DISTANCE and \ref PRINT).
*/
//+ENDPLUMEDOC

class Flush:
  public ActionPilot
{
public:
  explicit Flush(const ActionOptions&ao):
    Action(ao),
    ActionPilot(ao)
  {
    checkRead();
  }
  static void registerKeywords( Keywords& keys );
  void calculate() override {}
  void apply() override {}
  void update() override {
    plumed.fflush();
    log.flush();
    const ActionSet & actionSet(plumed.getActionSet());
    for(const auto & p : actionSet)
      p->fflush();
  }
};

PLUMED_REGISTER_ACTION(Flush,"FLUSH")

void Flush::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  keys.add("compulsory","STRIDE","the frequency with which all the open files should be flushed");
  keys.remove("LABEL");
}

}
}


