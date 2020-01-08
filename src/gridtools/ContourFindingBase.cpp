/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2020 The plumed team
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
#include "ContourFindingBase.h"

namespace PLMD {
namespace gridtools {

void ContourFindingBase::registerKeywords( Keywords& keys ) {
  ActionWithInputGrid::registerKeywords( keys );
  keys.add("compulsory","CONTOUR","the value we would like to draw the contour at in the space");
  keys.remove("KERNEL"); keys.remove("BANDWIDTH");
}

ContourFindingBase::ContourFindingBase(const ActionOptions&ao):
  Action(ao),
  ActionWithInputGrid(ao),
  mymin(this)
{
  if( ingrid->noDerivatives() ) error("cannot find contours if input grid has no derivatives");
  parse("CONTOUR",contour);
  log.printf("  calculating dividing surface along which function equals %f \n", contour);
}

}
}
