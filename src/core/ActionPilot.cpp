/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#include "ActionPilot.h"

using namespace std;
namespace PLMD {

void ActionPilot::registerKeywords(Keywords& keys) {
}

ActionPilot::ActionPilot(const ActionOptions&ao):
  Action(ao),
  stride(1)
{
  if( keywords.exists("STRIDE") ) {
    parse("STRIDE",stride);
    if( !keywords.style("STRIDE","hidden") ) log.printf("  with stride %d\n",stride);
  } else {
    stride=0;
  }
}

bool ActionPilot::onStep()const {
  if( stride>0 ) return getStep()%stride==0;
  return false;
}

int ActionPilot::getStride()const {
  return stride;
}

void ActionPilot::setStride( const int& n ) {
  stride=n;
}

}


