/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2023 The plumed team
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

//+PLUMEDOC MCOLVAR MFILTER_MORE
/*
Apply one minus a switching function to the input vector.

This action has been depracated as it is equivalent to [MORE_THAN](MORE_THAN.md).

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

class MFilterMore : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit MFilterMore(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(MFilterMore,"MFILTER_MORE")

void MFilterMore::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.setDeprecated("MORE_THAN");
  keys.add("compulsory","DATA","the vector you wish to transform");
  keys.add("compulsory","SWITCH","the switching function that transform");
  keys.addDeprecatedFlag("LOWMEM","");
  keys.addFlag("HIGHEST",false,"this flag allows you to recover the highest of these variables.");
  keys.addOutputComponent("highest","HIGHEST","scalar","the largest of the colvars");
  keys.needsAction("CUSTOM");
  keys.needsAction("GROUP");
  keys.needsAction("MORE_THAN");
  keys.needsAction("HIGHEST");
}

MFilterMore::MFilterMore(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  warning("This action has been depracated.  Look at the log to see how the same result is achieved with the new syntax");
  bool lowmem;
  parseFlag("LOWMEM",lowmem);
  if( lowmem ) {
    warning("LOWMEM flag is deprecated and is no longer required for this action");
  }
  std::string dd;
  parse("DATA",dd);
  std::string swit;
  parse("SWITCH",swit);
  readInputLine( getShortcutLabel() + "_grp: GROUP ATOMS=" + dd + "_grp");
  readInputLine( getShortcutLabel() + ": MORE_THAN ARG=" + dd + " SWITCH={" + swit + "}");
  bool highest;
  parseFlag("HIGHEST",highest);
  if( highest ) {
    readInputLine( getShortcutLabel() + "_filtered: CUSTOM ARG=" + dd + "," + getShortcutLabel() + " FUNC=x*y PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_highest: HIGHEST ARG=" + getShortcutLabel() + "_filtered");
  }
}

}
}
