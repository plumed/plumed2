/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013,2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "AnalysisWithDataCollection.h"
#include "core/ActionRegister.h"
#include "reference/ReferenceConfiguration.h"

namespace PLMD {
namespace analysis {

class EuclideanDissimilarityMatrix : public AnalysisWithDataCollection {
private:
  Matrix<double> dissimilarities;
public:
  static void registerKeywords( Keywords& keys );
  EuclideanDissimilarityMatrix( const ActionOptions& ao );
/// Do the analysis
  void performAnalysis();
/// This ensures that classes that use this data know that dissimilarities were set
  bool dissimilaritiesWereSet() const { return true; }
/// Get the squared dissimilarity between two reference configurations
  double getDissimilarity( const unsigned& i, const unsigned& j );
/// This is just to deal with ActionWithVessel
  void performTask( const unsigned& , const unsigned& , MultiValue& ) const { plumed_error(); }
};

PLUMED_REGISTER_ACTION(EuclideanDissimilarityMatrix,"EUCLIDEAN_DISSIMILARITIES")

void EuclideanDissimilarityMatrix::registerKeywords( Keywords& keys ){
  AnalysisWithDataCollection::registerKeywords( keys );
}

EuclideanDissimilarityMatrix::EuclideanDissimilarityMatrix( const ActionOptions& ao ):
Action(ao),
AnalysisWithDataCollection(ao)
{
}

void EuclideanDissimilarityMatrix::performAnalysis(){
  // Resize dissimilarities matrix and set all elements to zero
  dissimilarities.resize( getNumberOfDataPoints(), getNumberOfDataPoints() ); dissimilarities=0;
}

double EuclideanDissimilarityMatrix::getDissimilarity( const unsigned& iframe, const unsigned& jframe ){
  plumed_dbg_assert( iframe<dissimilarities.nrows() && jframe<dissimilarities.ncols() );
  if( dissimilarities(iframe,jframe)>0. ){ return dissimilarities(iframe,jframe); }
  if( iframe!=jframe ){ 
     bool ignore;
     dissimilarities(iframe,jframe) = dissimilarities(jframe,iframe) = distance( getPbc(), getArguments(), getReferenceConfiguration(iframe,ignore), getReferenceConfiguration(jframe,ignore), true ); 
     return dissimilarities(iframe,jframe);
  }
  return 0.0;
}

}
}
