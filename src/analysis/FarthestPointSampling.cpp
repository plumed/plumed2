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
public:
  static void registerKeywords( Keywords& keys );
  FarthestPointSampling( const ActionOptions& ao );
  void selectLandmarks();
};

PLUMED_REGISTER_ACTION(FarthestPointSampling,"LANDMARK_SELECT_FPS")

void FarthestPointSampling::registerKeywords( Keywords& keys ) {
  LandmarkSelectionBase::registerKeywords(keys);
  keys.add("compulsory","SEED","1234","a random number seed");
}

FarthestPointSampling::FarthestPointSampling( const ActionOptions& ao ):
  Action(ao),
  LandmarkSelectionBase(ao)
{
  if( !dissimilaritiesWereSet() ) error("dissimilarities have not been calcualted in input action");
  parse("SEED",seed);
}

void FarthestPointSampling::selectLandmarks() {
  std::vector<unsigned> landmarks( getNumberOfDataPoints() );

  // Select first point at random
  Random random; random.setSeed(-seed); double rand=random.RandU01();
  landmarks[0] = std::floor( my_input_data->getNumberOfDataPoints()*rand );
  selectFrame( landmarks[0] );

  // Now find distance to all other points (N.B. We can use squared distances here for speed)
  Matrix<double> distances( getNumberOfDataPoints(), my_input_data->getNumberOfDataPoints() );
  for(unsigned i=0; i<my_input_data->getNumberOfDataPoints(); ++i) distances(0,i) = my_input_data->getDissimilarity( landmarks[0], i );

  // Now find all other landmarks
  for(unsigned i=1; i<getNumberOfDataPoints(); ++i) {
    // Find point that has the largest minimum distance from the landmarks selected thus far
    double maxd=0;
    for(unsigned j=0; j<my_input_data->getNumberOfDataPoints(); ++j) {
      double mind=distances(0,j);
      for(unsigned k=1; k<i; ++k) {
        if( distances(k,j)<mind ) { mind=distances(k,j); }
      }
      if( mind>maxd ) { maxd=mind; landmarks[i]=j; }
    }
    selectFrame( landmarks[i] );
    for(unsigned k=0; k<my_input_data->getNumberOfDataPoints(); ++k) distances(i,k) = my_input_data->getDissimilarity( landmarks[i], k );
  }
}

}
}
