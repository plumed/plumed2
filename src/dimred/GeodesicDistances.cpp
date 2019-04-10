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
#include "analysis/AnalysisBase.h"
#include "tools/PDB.h"
#include "core/ActionRegister.h"

//+PLUMEDOC ANALYSIS GEODESIC_DISTANCES
/*
Calculate the matrix of geodesic distances between a trajectory of atomic configurations.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace dimred {

class GeodesicDistances : public analysis::AnalysisBase {
private:
  double eps;
  unsigned nneigh;
  Matrix<double> dissimilarities;
public:
  static void registerKeywords( Keywords& keys );
  GeodesicDistances( const ActionOptions& ao );
/// Do the analysis
  void performAnalysis();
/// This ensures that classes that use this data know that dissimilarities were set
  bool dissimilaritiesWereSet() const { return true; }
/// Get the squared dissimilarity between two reference configurations
  double getDissimilarity( const unsigned& i, const unsigned& j );
/// This is just to deal with ActionWithVessel
  void performTask( const unsigned&, const unsigned&, MultiValue& ) const { plumed_error(); }
};

PLUMED_REGISTER_ACTION(GeodesicDistances,"GEODESIC_DISTANCES")

void GeodesicDistances::registerKeywords( Keywords& keys ) {
  AnalysisBase::registerKeywords( keys ); 
  keys.add("optional","EPSILON","points that are separated by less than this value are connected in the graph");
  keys.add("compulsory","NNEIGHBORS","0","connect each point to its k nearest neighbors");
}

GeodesicDistances::GeodesicDistances( const ActionOptions& ao ):
  Action(ao),
  AnalysisBase(ao),
  eps(0)
{
  if(!dissimilaritiesWereSet() ) error("dissimilarities have not been calcualted in input actions");
  log.printf("  calculating geodesic distances based on dissimilarity matrix \n");
  nneigh=0; parse("NNEIGHBORS",nneigh);
  if( nneigh==0 ) {
      parse("EPSILON",eps);
      if( eps<epsilon ) error("did not find NNEIGHBORS or EPSILON keyword in input");
      log.printf("  costructing graph by connecting points that are %f appart \n", eps );
  } else {
      log.printf("  constructing graph by connecting each point to its %u nearest neighbors \n", nneigh );
  }
}

void GeodesicDistances::performAnalysis() {
  // Resize dissimilarities matrix and set all elements to zero
  if( !usingLowMem() ) {
    dissimilarities.resize( getNumberOfDataPoints(), getNumberOfDataPoints() ); dissimilarities=0;
  }
  // Calculate the neighbors of each of the points
  if( nneigh>0 ) {}
  else {}
}

double GeodesicDistances::getDissimilarity( const unsigned& iframe, const unsigned& jframe ) {
  plumed_dbg_assert( iframe<getNumberOfDataPoints() && jframe<getNumberOfDataPoints() );
  if( !usingLowMem() ) {
      dissimilarities(iframe,jframe)=dissimilarities(jframe,iframe)=my_input_data->getDissimilarity( iframe, jframe );
      return dissimilarities(iframe,jframe);
  }
  // An implementation of Djikestras algorithm
  return my_input_data->getDissimilarity( iframe, jframe );
}

}
}
