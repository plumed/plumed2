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
#include "ActionForInterface.h"
#include "ActionRegister.h"
#include "PlumedMain.h"
#include "ActionSet.h"

namespace PLMD {

void ActionForInterface::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys);
  ActionWithValue::registerKeywords( keys );
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  keys.add("optional","ROLE","Get the role this value plays in the code can be x/y/z/m/q to signify that this is x, y, z positions of atoms or masses or charges of atoms");
}

ActionForInterface::ActionForInterface(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  firststep(true),
  wasscaled(false),
  wasset(false) {
  if( keywords.exists("ROLE") && getName()!="DOMAIN_DECOMPOSITION") {
    parse("ROLE",role);
  }
}

std::string ActionForInterface::getRole() const {
  return role;
}

}
