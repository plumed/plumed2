/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2017 The plumed team
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
#include "MatrixProductBase.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace adjmat {

class Dot : public MatrixProductBase {
public:
  static void registerKeywords( Keywords& keys );
  explicit Dot(const ActionOptions&);
  void setupForTask( const unsigned& current, MultiValue& myvals, std::vector<unsigned> & indices, std::vector<Vector>& atoms ) const override;
  double computeVectorProduct( const unsigned& index1, const unsigned& index2,
                               const std::vector<double>& vec1, const std::vector<double>& vec2,
                               std::vector<double>& dvec1, std::vector<double>& dvec2, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(Dot,"DOT")

void Dot::registerKeywords( Keywords& keys ) {
  MatrixProductBase::registerKeywords( keys );
}

Dot::Dot(const ActionOptions& ao):
  Action(ao),
  MatrixProductBase(ao)
{
  setNotPeriodic();
}

void Dot::setupForTask( const unsigned& current, MultiValue& myvals, std::vector<unsigned> & indices, std::vector<Vector>& atoms ) const {
  if( getPntrToArgument(1)->getRank()==1 ) {
      unsigned nr = 0, nvals = getPntrToArgument(1)->getShape()[0]; 
      for(unsigned i=0;i<nvals;++i) {
          if( skip_ieqj && myvals.getTaskIndex()==i ) continue;
          if( fabs( getPntrToArgument(1)->get(i) )>epsilon ) nr++;
      }
      if( indices.size()!=nr ) indices.resize( 1+nr );
      indices[0]=current; unsigned k = 1;
      for(unsigned i=0;i<nvals;++i) {
          if( skip_ieqj && myvals.getTaskIndex()==i ) continue;
          if( fabs( getPntrToArgument(1)->get(i) )>epsilon ) { indices[k]=i; k++; }
      }
      myvals.setSplitIndex( indices.size() ); myvals.setNumberOfIndices( indices.size() ); 
      return;
  } 
  MatrixProductBase::setupForTask( current, myvals, indices, atoms );
}

double Dot::computeVectorProduct( const unsigned& index1, const unsigned& index2,
    const std::vector<double>& vec1, const std::vector<double>& vec2,
    std::vector<double>& dvec1, std::vector<double>& dvec2, MultiValue& myvals ) const {
  double val=0;
  for(unsigned i=0; i<vec1.size(); ++i) {
    val += vec1[i]*vec2[i];
    dvec1[i]=vec2[i]; dvec2[i]=vec1[i];
  }
  return val;
}

}
}
