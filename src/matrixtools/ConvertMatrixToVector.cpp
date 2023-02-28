/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2017 The plumed team
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
#include "ActionWithInputMatrices.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace matrixtools {

class ConvertMatrixToVector : public ActionWithInputMatrices {
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit ConvertMatrixToVector(const ActionOptions&);
///
  void buildTaskListFromArgumentRequests( const unsigned& ntasks, bool& reduce, std::set<AtomNumber>& otasks ) override ;
/// Do the calculation
  void completeMatrixOperations() override;
///
  void apply() override;
///
  std::vector<unsigned> getValueShapeFromArguments() override;
///
  double getForceOnMatrixElement( const unsigned& imat, const unsigned& jrow, const unsigned& kcol ) const override;
};

PLUMED_REGISTER_ACTION(ConvertMatrixToVector,"FLATTEN")

void ConvertMatrixToVector::registerKeywords( Keywords& keys ) {
  ActionWithInputMatrices::registerKeywords( keys );
}

ConvertMatrixToVector::ConvertMatrixToVector(const ActionOptions& ao):
  Action(ao),
  ActionWithInputMatrices(ao)
{
  if( getNumberOfArguments()!=1 ) error("should only be one argument for this action");
  std::vector<unsigned> shape( getValueShapeFromArguments() ); addValue( shape );
}

std::vector<unsigned> ConvertMatrixToVector::getValueShapeFromArguments() {
  std::vector<unsigned> shape( getPntrToArgument(0)->getShape() ), nshape(1); 
  if( shape.size()!=2 ) error("input should be a matrix");
  nshape[0] = shape[0]*shape[1]; return nshape; 
}

void ConvertMatrixToVector::buildTaskListFromArgumentRequests( const unsigned& ntasks, bool& reduce, std::set<AtomNumber>& otasks ) {
  // By calling this function here we ensure that if we are multiplying if we are calculating an lq6 in a particular volume
  // the calculation runs fast.  Only the parts of the matrix that we need to compute the dot products in the volume of interest
  // are computed.  There are lots of zeros in the matrix because we don't need all the elements in this case
  propegateTaskListsForValue( 0, ntasks, reduce, otasks );
}

void ConvertMatrixToVector::completeMatrixOperations() {
  unsigned nedge=0; retrieveEdgeList( 0, nedge ); std::vector<unsigned> shape( getPntrToArgument(0)->getShape() );
  for(unsigned i=0; i<nedge;++i) getPntrToOutput(0)->set( pairs[i].first*shape[1] + pairs[i].second, vals[i] ); 
  if( getPntrToArgument(0)->isSymmetric() ) {
      for(unsigned i=0; i<nedge;++i) getPntrToOutput(0)->set( pairs[i].second*shape[1] + pairs[i].first, vals[i] );
  }
}

void ConvertMatrixToVector::apply() {
  // Apply force on the matrix
  applyForceOnMatrix(0);
}

double ConvertMatrixToVector::getForceOnMatrixElement( const unsigned& imat, const unsigned& jrow, const unsigned& kcol ) const {
  return getPntrToOutput(0)->getForce(jrow*getPntrToArgument(0)->getShape()[1]+kcol);
}


}
}
