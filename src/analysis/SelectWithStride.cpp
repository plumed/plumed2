/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2019 The plumed team
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
#include "LandmarkSelectionBase.h"
#include "core/ActionRegister.h"

//+PLUMEDOC LANDMARKS LANDMARK_SELECT_STRIDE
/*
Select every \f$k\f$th landmark from the trajectory.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace analysis {

class SelectWithStride : public LandmarkSelectionBase {
public:
  static void registerKeywords( Keywords& keys );
  explicit SelectWithStride( const ActionOptions& ao );
  void selectLandmarks();
};

PLUMED_REGISTER_ACTION(SelectWithStride,"LANDMARK_SELECT_STRIDE")

void SelectWithStride::registerKeywords( Keywords& keys ) {
  LandmarkSelectionBase::registerKeywords( keys );
}

SelectWithStride::SelectWithStride( const ActionOptions& ao ):
  Action(ao),
  LandmarkSelectionBase(ao)
{
}

void SelectWithStride::selectLandmarks() {
  unsigned stride = std::floor( my_input_data->getNumberOfDataPoints() / getNumberOfDataPoints() ), max=stride*getNumberOfDataPoints();
  for(unsigned i=0; i<max; i+=stride) selectFrame( i );
}

}
}
