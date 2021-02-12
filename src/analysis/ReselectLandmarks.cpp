/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "tools/Random.h"

//+PLUMEDOC LANDMARKS RESELECT_LANDMARKS
/*
This allows us to use one measure in landmark selection and a different measure in dimensionality reduction

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace analysis {

class ReselectLandmarks : public LandmarkSelectionBase {
private:
  LandmarkSelectionBase* mylandmarks;
public:
  static void registerKeywords( Keywords& keys );
  explicit ReselectLandmarks( const ActionOptions& ao );
  void selectLandmarks() override;
};

PLUMED_REGISTER_ACTION(ReselectLandmarks,"RESELECT_LANDMARKS")

void ReselectLandmarks::registerKeywords( Keywords& keys ) {
  LandmarkSelectionBase::registerKeywords(keys);
  keys.remove("NLANDMARKS");
  keys.add("compulsory","LANDMARKS","the action that selects the landmarks that you want to reselect using this action");
}

ReselectLandmarks::ReselectLandmarks( const ActionOptions& ao ):
  Action(ao),
  LandmarkSelectionBase(ao)
{
  std::string datastr; parse("LANDMARKS",datastr);
  mylandmarks = plumed.getActionSet().selectWithLabel<LandmarkSelectionBase*>( datastr );
  if( !mylandmarks ) error("input to LANDMARKS is not a landmark selection action");
  nlandmarks = mylandmarks->nlandmarks;

  if( (mylandmarks->my_input_data)->getNumberOfDataPoints()!=my_input_data->getNumberOfDataPoints() ) error("mismatch between amount of landmark class and base class");
}

void ReselectLandmarks::selectLandmarks() {
  for(unsigned i=0; i<mylandmarks->getNumberOfDataPoints(); ++i) selectFrame( mylandmarks->getDataPointIndexInBase(i) );
}

}
}
