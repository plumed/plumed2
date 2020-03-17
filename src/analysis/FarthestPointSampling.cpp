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
#include "tools/Random.h"

//+PLUMEDOC LANDMARKS LANDMARK_SELECT_FPS
/*
Select a set of landmarks using farthest point sampling.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace analysis {

class FarthestPointSampling : public LandmarkSelectionBase {
private:
  unsigned seed;
  Value* distance_matrix;
  Value* outmatrix;
public:
  static void registerKeywords( Keywords& keys );
  explicit FarthestPointSampling( const ActionOptions& ao );
  void selectLandmarks() override;
};

PLUMED_REGISTER_ACTION(FarthestPointSampling,"LANDMARK_SELECT_FPS")

void FarthestPointSampling::registerKeywords( Keywords& keys ) {
  LandmarkSelectionBase::registerKeywords(keys);
  keys.add("compulsory","SEED","1234","a random number seed");
}

FarthestPointSampling::FarthestPointSampling( const ActionOptions& ao ):
  Action(ao),
  LandmarkSelectionBase(ao),
  distance_matrix(NULL),
  outmatrix(NULL)
{
  for(unsigned i=0;i<getNumberOfArguments();++i) {
      if( !distance_matrix && getPntrToArgument(i)->getRank()==2 ) { distance_matrix = getPntrToArgument(i); outmatrix = getPntrToOutput(i); }
      else if( distance_matrix && getPntrToArgument(i)->getRank()==2 ) error("should only be one distance matrix in input arguments");
  }
  if( !distance_matrix ) error("could not find distance matrix in input arguments");
  parse("SEED",seed); plumed_assert( outmatrix->getRank()==2 );
}

void FarthestPointSampling::selectLandmarks() {
  // Select first point at random
  Random random; random.setSeed(-seed); selectFrame( std::floor( nvals*random.RandU01() )  );

  // Now find all other landmarks
  for(unsigned i=1; i<nlandmarks; ++i) {
    // Find point that has the largest minimum distance from the landmarks selected thus far
    double maxd=0; unsigned lll=0;
    for(unsigned j=0; j<nvals; ++j) {
      double mind=outmatrix->get( j );
      for(unsigned k=1; k<i; ++k) {
        double tdist = outmatrix->get( k*nvals + j );
        if( tdist<mind ) { mind=tdist; }
      }
      if( mind>maxd ) { maxd=mind; lll=j; }
    }
    selectFrame( lll );
  }
}

}
}
