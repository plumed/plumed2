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
#include "tools/Random.h"

//+PLUMEDOC LANDMARKS LANDMARK_SELECT_RANDOM
/*
Select a random set of landmarks from a large set of configurations.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace analysis {

class SelectRandomFrames : public LandmarkSelectionBase {
private:
  unsigned seed;
public:
  static void registerKeywords( Keywords& keys );
  explicit SelectRandomFrames( const ActionOptions& ao );
  void selectLandmarks() override;
};

PLUMED_REGISTER_ACTION(SelectRandomFrames,"LANDMARK_SELECT_RANDOM")

void SelectRandomFrames::registerKeywords( Keywords& keys ) {
  LandmarkSelectionBase::registerKeywords(keys);
  keys.add("compulsory","SEED","1234","a random number seed");
}

SelectRandomFrames::SelectRandomFrames( const ActionOptions& ao ):
  Action(ao),
  LandmarkSelectionBase(ao)
{
  parse("SEED",seed);
}

void SelectRandomFrames::selectLandmarks() {
  Random r; r.setSeed(-seed);
  unsigned nframe=my_input_data->getNumberOfDataPoints();
  unsigned nland=getNumberOfDataPoints();

  std::vector<bool> selected( nframe, false );

  unsigned fcount=0;
  while (fcount<nland) {
    unsigned iframe = std::floor( r.U01()*nframe );
    if (!selected[iframe]) {
      selected[iframe]=true;
      selectFrame( iframe );
      ++fcount;
    }
  }
}

}
}
