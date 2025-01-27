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
namespace contour {

void ContourFindingBase::registerKeywords( Keywords& keys ) {
  gridtools::ActionWithGrid::registerKeywords( keys );
  keys.addInputKeyword("compulsory","ARG","grid","the labels of the grid in which the contour will be found");
  keys.add("compulsory","CONTOUR","the value we would like to draw the contour at in the space");
  gridtools::EvaluateGridFunction gg;
  gg.registerKeywords(keys);
}

ContourFindingBase::ContourFindingBase(const ActionOptions&ao):
  Action(ao),
  ActionWithGrid(ao),
  mymin(this) {
  parse("CONTOUR",contour);
  function.read( this );
  log.printf("  calculating dividing surface along which function equals %f \n", contour);
}

void ContourFindingBase::setupOnFirstStep( const bool incalc ) {
  function.setup( this );
  setupValuesOnFirstStep();
}

}
}
