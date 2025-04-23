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
#include "PbcAction.h"
#include "DomainDecomposition.h"
#include "tools/Pbc.h"
#include "PlumedMain.h"
#include "ActionSet.h"
#include "ActionRegister.h"

//+PLUMEDOC ANALYSIS PBC
/*
Pass the cell vectors into PLUMED.

This command is used to pass the cell vectors from an MD code to PLUMED. The command is a specialised
version of the [PUT](PUT.md) command that accepts a $3\times 3$ matrix.  The additional specialisation is
required because the cell vectors play a special role when you do calculations with PLUMED.

Notice that when you create a [DOMAIN_DECOMPOSITION](DOMAIN_DECOMPOSITION.md) action a PBC action is created
automatically.

You would only use the PBC command if you were calling PLUMED from python or an MD code. There is no reason for using this command in a conventional PLUMED
input file like this:

```plumed
pbc: PBC
```

Instead this command would be called from python like this:

```python
import plumed

# Create a PLUMED object
p = plumed.Plumed()
p.cmd("readInputLine", "pbc: PBC")
```

*/
//+ENDPLUMEDOC

namespace PLMD {

PLUMED_REGISTER_ACTION(PbcAction,"PBC")

void PbcAction::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  keys.setValueDescription("matrix","a matrix containing the cell vectors that were passed from the MD code");
}

PbcAction::PbcAction(const ActionOptions&ao):
  Action(ao),
  ActionToPutData(ao),
  interface(NULL) {
  std::vector<std::size_t> shape(2);
  shape[0]=shape[1]=3;
  addValue( shape );
  setNotPeriodic();
  setUnit( "length", "energy" );
  getPntrToValue()->reshapeMatrixStore(3);
}


void PbcAction::setPbc() {
  if( !interface ) {
    std::vector<DomainDecomposition*> allput=plumed.getActionSet().select<DomainDecomposition*>();
    if( allput.size()>1 ) {
      warning("found more than one interface so don't know how to broadcast cell");
    }
    interface = allput[0];
  }
  Tensor box;
  if( interface ) {
    interface->broadcastToDomains( getPntrToValue() );
  }
  for(unsigned i=0; i<3; ++i)
    for(unsigned j=0; j<3; ++j) {
      box(i,j) = getPntrToValue()->get(3*i+j);
    }
  pbc.setBox(box);
}

void PbcAction::wait() {
  ActionToPutData::wait();
  setPbc();
}

void PbcAction::readBinary(std::istream&i) {
  ActionToPutData::readBinary(i);
  setPbc();
}

}



