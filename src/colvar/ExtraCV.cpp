/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2023 The plumed team
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
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR EXTRACV
/*
Allow PLUMED to use collective variables computed in the MD engine.

This feature requires the MD engine to use special instructions to pass to PLUMED the value of
some pre-computed collective variable. Check the documentation for the MD code to find out which
collective variables can be computed and passed to PLUMED. These variables can then be accessed by
name using the EXTRACV action.

This example takes the lambda variable pre-computed in GROMACS and applies a restraint to keep
it close to the value 3.

```plumed
l: EXTRACV NAME=lambda
RESTRAINT ARG=l KAPPA=10 AT=3
```


*/
//+ENDPLUMEDOC


class ExtraCV : public ActionShortcut {
public:
  explicit ExtraCV(const ActionOptions&);
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(ExtraCV,"EXTRACV")

ExtraCV::ExtraCV(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  std::vector<std::string> argn(1);
  parse("NAME",argn[0]);
  readInputLine( argn[0] + ": PUT UNIT=number SHAPE=0 MUTABLE PERIODIC=NO");
  if( getShortcutLabel()!=argn[0] ) {
    readInputLine( getShortcutLabel() + ": COMBINE ARG=" + argn[0] + " PERIODIC=NO");
  }
}

void ExtraCV::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","NAME","name of the CV as computed by the MD engine");
  keys.setValueDescription("scalar","the value of the CV that was passed from the MD code to PLUMED");
  keys.needsAction("PUT");
  keys.needsAction("COMBINE");
}

}
}



