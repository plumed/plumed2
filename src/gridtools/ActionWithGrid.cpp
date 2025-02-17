/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#include "ActionWithGrid.h"

namespace PLMD {
namespace gridtools {

ActionWithGrid* ActionWithGrid::getInputActionWithGrid( Action* action ) {
  ActionWithGrid* ag = dynamic_cast<ActionWithGrid*>( action );
  if( !ag && action->getName()=="ACCUMULATE" ) {
    ActionWithArguments* aa=dynamic_cast<ActionWithArguments*>( action );
    plumed_assert( aa );
    ag = dynamic_cast<ActionWithGrid*>( (aa->getPntrToArgument(0))->getPntrToAction() );
  }
  plumed_assert( ag );
  return ag;
}

void ActionWithGrid::registerKeywords( Keywords& keys ) {
  ActionWithVector::registerKeywords( keys );
}

ActionWithGrid::ActionWithGrid(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  firststep(true) {
}

void ActionWithGrid::calculate() {
  plumed_assert( !actionInChain() );
  if( firststep ) {
    setupOnFirstStep( true );
    firststep=false;
  }

  runAllTasks();
}

}
}
