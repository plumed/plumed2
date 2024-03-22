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
Basically equivalent to LESS_THAN.

This action has been depracated

\par Examples


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
  keys.add("compulsory","DATA","the vector you wish to transform");
  keys.add("compulsory","SWITCH","the switching function that transform");
  keys.needsAction("GROUP"); keys.needsAction("MORE_THAN");
}

MFilterMore::MFilterMore(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao)
{
  warning("This action has been depracated.  Look at the log to see how the same result is achieved with the new syntax");
  std::string dd; parse("DATA",dd);
  std::string swit; parse("SWITCH",swit);
  readInputLine( getShortcutLabel() + "_grp: GROUP ATOMS=" + dd + "_grp");
  readInputLine( getShortcutLabel() + ": MORE_THAN ARG=" + dd + " SWITCH={" + swit + "}");
}

}
}
