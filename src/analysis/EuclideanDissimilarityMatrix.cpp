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
#include "Analysis.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace analysis {

class EuclideanDissimilarityMatrix : public Analysis {
public:
  static void registerKeywords( Keywords& keys );
  EuclideanDissimilarityMatrix( const ActionOptions& ao );
  void calcDissimilarity( const unsigned& , const unsigned& );
  void performAnalysis(){}
  void performTask(){ plumed_error(); }
  bool dissimilaritiesWereSet() const { return true; }
};

PLUMED_REGISTER_ACTION(EuclideanDissimilarityMatrix,"EUCLIDEAN_DISSIMILARITIES")

void EuclideanDissimilarityMatrix::registerKeywords( Keywords& keys ){
  Analysis::registerKeywords( keys );
}

EuclideanDissimilarityMatrix::EuclideanDissimilarityMatrix( const ActionOptions& ao ):
Action(ao),
Analysis(ao)
{
}

void EuclideanDissimilarityMatrix::calcDissimilarity( const unsigned& iframe, const unsigned& jframe ){
  double d = getDistanceBetweenFrames( iframe, jframe, true );
  setDissimilarityMatrixElement( iframe, jframe, d );
}

}
}
