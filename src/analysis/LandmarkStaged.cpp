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
#include <iostream>

//+PLUMEDOC LANDMARKS LANDMARK_SELECT_STAGED
/*
Select a set of landmarks using the staged algorithm.

\par Examples

*/
//+ENDPLUMEDOC


namespace PLMD {
namespace analysis {

class LandmarkStaged : public LandmarkSelectionBase {
private:
  unsigned seed;
  double gamma;
public:
  static void registerKeywords( Keywords& keys );
  explicit LandmarkStaged( const ActionOptions& ao );
  void selectLandmarks() override;
};

PLUMED_REGISTER_ACTION(LandmarkStaged,"LANDMARK_SELECT_STAGED")

void LandmarkStaged::registerKeywords( Keywords& keys ) {
  LandmarkSelectionBase::registerKeywords(keys);
  keys.add("compulsory","GAMMA","the gamma parameter to be used in weights");
  keys.add("compulsory","SEED","1234","a random number seed");
}

LandmarkStaged::LandmarkStaged( const ActionOptions& ao ):
  Action(ao),
  LandmarkSelectionBase(ao)
{
  parse("SEED",seed); parse("GAMMA",gamma);
  log.printf("  probability of selecting voronoi polyhedra equal to exp(-weight/%f) \n", gamma );
}

void LandmarkStaged::selectLandmarks() {
  unsigned int n = getNumberOfDataPoints(); // The number of landmarks to pick
  unsigned int N = my_input_data->getNumberOfDataPoints();  // The total number of frames we can choose from
  unsigned int m = static_cast<int>( sqrt(n*N) );
  std::vector<unsigned> fpslandmarks(m);
  // Select first point at random
  Random random; random.setSeed(-seed); double rand=random.RandU01();
  fpslandmarks[0] = std::floor( N*rand );

  // using FPS we want to find m landmarks where m = sqrt(nN)
  // Now find distance to all other points
  Matrix<double> distances( m, N );
  for(unsigned int i=0; i<N; ++i) {
    distances(0,i) = my_input_data->getDissimilarity( fpslandmarks[0], i );
  }

  // Now find all other landmarks
  for(unsigned i=1; i<m; ++i) {
    // Find point that has the largest minimum distance from the landmarks selected thus far
    double maxd=0;
    for(unsigned j=0; j<N; ++j) {
      double mind=distances(0,j);
      for(unsigned k=1; k<i; ++k) {
        if( distances(k,j)<mind ) { mind=distances(k,j); }
      }
      if( mind>maxd ) { maxd=mind; fpslandmarks[i]=j; }
    }
    for(unsigned k=0; k<N; ++k) distances(i,k) = my_input_data->getDissimilarity( fpslandmarks[i], k );
  }

  // Initial FPS selection of m landmarks completed
  // Now find voronoi weights of these m points
  std::vector<unsigned> poly_assign( N );
  std::vector<double> weights( m, 0 );
  voronoiAnalysis( fpslandmarks, weights, poly_assign );

  //Calulate total weight of voronoi polyhedras
  double vweight=0; for(unsigned i=0; i<m; i++) vweight += exp( -weights[i] / gamma );

  std::vector<bool> selected(N, false); unsigned ncount=0;
  while ( ncount<n) {
//  generate random number and check which point it belongs to. select only it was not selected before
    double rand = vweight*random.RandU01();
    double running_vweight=0;
    for(unsigned jpoly=0; jpoly<m; ++jpoly) {
      running_vweight+=exp( -weights[jpoly] / gamma );
      if( running_vweight>=rand ) {
        double tweight=0;
        for(unsigned i=0; i<poly_assign.size(); ++i) {
          if( poly_assign[i]==jpoly ) tweight += getWeight( i );
        }
        double rand_poly = tweight*random.RandU01();
        double running_tweight=0;
        for(unsigned i=0; i<N; ++i) {
          if( poly_assign[i]==jpoly ) {
            running_tweight += getWeight( i );
            if( running_tweight>=rand_poly && !selected[i] ) {
              selectFrame(i); selected[i]=true; ncount++; break;
            } else if( running_tweight>=rand_poly ) {
              break;
            }
          }
        }
      }
    }
  }
}

}
}
