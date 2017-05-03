/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#include "ActionSetup.h"
#include "PlumedMain.h"
#include "ActionSet.h"
#include "tools/Exception.h"

namespace PLMD {

ActionSetup::ActionSetup(const ActionOptions&ao):
  Action(ao)
{
  const ActionSet& actionset(plumed.getActionSet());
  for(const auto & p : actionset) {
// check that all the preceeding actions are ActionSetup
    if( !dynamic_cast<ActionSetup*>(p) ) error("Action " + getLabel() + " is a setup action, and should be only preceeded by other setup actions");
  }
}

void ActionSetup::registerKeywords( Keywords& keys ) {
  Action::registerKeywords(keys);
  keys.remove("LABEL");
}

}
