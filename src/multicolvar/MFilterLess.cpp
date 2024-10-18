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

//+PLUMEDOC MCOLVAR MFILTER_LESS
/*
Apply a switching function to the input vector.

This action has been depracated as it is equivalent to [LESS_THAN](LESS_THAN.md).

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

class MFilterLess : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit MFilterLess(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(MFilterLess,"MFILTER_LESS")

void MFilterLess::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.setDeprecated("LESS_THAN");
  keys.add("compulsory","DATA","the vector you wish to transform");
  keys.add("compulsory","SWITCH","the switching function that transform");
  keys.setValueDescription("vector","a vector that has the same dimension as the input vector with elements equal to one if the corresponding component of the vector is less than a tolerance and zero otherwise");
  keys.needsAction("GROUP");
  keys.needsAction("LESS_THAN");
}

MFilterLess::MFilterLess(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  warning("This action has been depracated.  Look at the log to see how the same result is achieved with the new syntax");
  std::string dd;
  parse("DATA",dd);
  std::string swit;
  parse("SWITCH",swit);
  readInputLine( getShortcutLabel() + "_grp: GROUP ATOMS=" + dd + "_grp");
  readInputLine( getShortcutLabel() + ": LESS_THAN ARG=" + dd + " SWITCH={" + swit + "}");
}

}
}
