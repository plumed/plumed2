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
#include "VectorProductMatrix.h"
#include "core/PlumedMain.h"
#include "core/ActionRegister.h"

//+PLUMEDOC ANALYSIS DISSIMILARITIES
/*
Calculate the matrix of dissimilarities between a trajectory of atomic configurations.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class DissimilarityMatrix : public VectorProductMatrix { 
private:
  bool squared;
public:
  static void registerKeywords( Keywords& keys );
  DissimilarityMatrix( const ActionOptions& ao );
  double computeVectorProduct( const unsigned& index1, const unsigned& index2,
                               const std::vector<double>& vec1, const std::vector<double>& vec2,
                               std::vector<double>& dvec1, std::vector<double>& dvec2, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(DissimilarityMatrix,"DISSIMILARITIES")

void DissimilarityMatrix::registerKeywords( Keywords& keys ) {
  VectorProductMatrix::registerKeywords( keys );
  keys.addFlag("SQUARED",false,"calculate the square of the dissimilarity matrix");
}

DissimilarityMatrix::DissimilarityMatrix( const ActionOptions& ao ):
  Action(ao),
  VectorProductMatrix(ao)
{
  parseFlag("SQUARED",squared);
  if( squared ) log.printf("  computing the square of the dissimilarity matrix\n");
  forcesToApply.resize( getNumberOfDerivatives() );
  setNotPeriodic();
}

double DissimilarityMatrix::computeVectorProduct( const unsigned& index1, const unsigned& index2,
    const std::vector<double>& vec1, const std::vector<double>& vec2,
    std::vector<double>& dvec1, std::vector<double>& dvec2, MultiValue& myvals ) const {
  double dist = 0; 
  for(unsigned i=0;i<vec1.size();++i) {
       double tmp = getPntrToArgument(arg_ends[i])->difference(vec2[i], vec1[i]); dist += tmp*tmp;
       dvec1[i] = 2*tmp; dvec2[i] = -2*tmp;
  }
  if( squared ) return dist ;
  dist = sqrt( dist );
  for(unsigned i=0;i<vec1.size();++i) { dvec1[i] /= 2*dist; dvec2[i] /= 2*dist; }
  return dist;
}

}
}
