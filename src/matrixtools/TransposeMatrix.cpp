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

class TransposeMatrix : public ActionWithInputMatrices {
private:
  std::vector<std::string> aLabelsInChain;
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit TransposeMatrix(const ActionOptions&);
///
  unsigned getNumberOfColumns() const override;
///
  void buildTaskListFromArgumentRequests( const unsigned& ntasks, bool& reduce, std::set<unsigned>& otasks ) override ;
/// Do the calculation
  void completeMatrixOperations() override;
///
  void apply() override;
///
  std::vector<unsigned> getValueShapeFromArguments() override;
///
  double getForceOnMatrixElement( const unsigned& imat, const unsigned& jrow, const unsigned& kcol ) const override;
};

PLUMED_REGISTER_ACTION(TransposeMatrix,"TRANSPOSE")

void TransposeMatrix::registerKeywords( Keywords& keys ) {
  ActionWithInputMatrices::registerKeywords( keys );
}

TransposeMatrix::TransposeMatrix(const ActionOptions& ao):
  Action(ao),
  ActionWithInputMatrices(ao)
{
  if( getNumberOfArguments()!=1 ) error("should only be one argument for this action");
  if( getPntrToArgument(0)->isSymmetric() ) warning("input matrix is symmetric.  Transposing will achieve nothing!");
  std::vector<unsigned> shape( getValueShapeFromArguments() ); addValue( shape ); 
}

std::vector<unsigned> TransposeMatrix::getValueShapeFromArguments() {
  std::vector<unsigned> shape(2);
  if( getPntrToArgument(0)->getRank()==0 ) {
     shape.resize(0);
  } else if( getPntrToArgument(0)->getRank()==1 ) {
     shape[0]=1;
     shape[1]=getPntrToArgument(0)->getShape()[0]; 
  } else {
     if( getPntrToArgument(0)->getShape()[0]==1 ) { 
         shape.resize(1); shape[0] = getPntrToArgument(0)->getShape()[1]; 
     } else {
         shape[0]=getPntrToArgument(0)->getShape()[1]; 
         shape[1]=getPntrToArgument(0)->getShape()[0];
     }
  }
  return shape;
}

void TransposeMatrix::buildTaskListFromArgumentRequests( const unsigned& ntasks, bool& reduce, std::set<unsigned>& otasks ) {
  // By calling this function here we ensure that if we are multiplying if we are calculating an lq6 in a particular volume
  // the calculation runs fast.  Only the parts of the matrix that we need to compute the dot products in the volume of interest
  // are computed.  There are lots of zeros in the matrix because we don't need all the elements in this case
  propegateTaskListsForValue( 0, ntasks, reduce, otasks );
}

unsigned TransposeMatrix::getNumberOfColumns() const {
  return getPntrToArgument(0)->getShape()[0];
}

void TransposeMatrix::completeMatrixOperations() {
  // Retrieve the non-zero pairs
  if( getPntrToArgument(0)->getRank()<=1 || getPntrToOutput(0)->getRank()==1 ) {
      unsigned nv=getPntrToArgument(0)->getNumberOfValues();
      for(unsigned i=0; i<nv; ++i) getPntrToOutput(0)->set( i, getPntrToArgument(0)->get(i) );
  } else {
      unsigned nedge=0; retrieveEdgeList( 0, nedge ); std::vector<unsigned> shape( getPntrToOutput(0)->getShape() );
      for(unsigned i=0; i<nedge;++i) getPntrToOutput(0)->set( pairs[i].second*shape[1] + pairs[i].first, vals[i] ); 
      if( getPntrToArgument(0)->isSymmetric() ) {
          for(unsigned i=0; i<nedge;++i) getPntrToOutput(0)->set( pairs[i].first*shape[1] + pairs[i].second, vals[i] );
      }
  }
}

void TransposeMatrix::apply() {
  // Apply force on the matrix
  if( getPntrToOutput(0)->forcesWereAdded() ) {
      if( getPntrToArgument(0)->getRank()<=1 || getPntrToOutput(0)->getRank()==1 ) {
          unsigned nv=getPntrToArgument(0)->getNumberOfValues();
          for(unsigned i=0; i<nv; ++i) getPntrToArgument(0)->addForce( i, getPntrToOutput(0)->getForce(i) );
      } else applyForceOnMatrix(0);
  }
}

double TransposeMatrix::getForceOnMatrixElement( const unsigned& imat, const unsigned& jrow, const unsigned& kcol ) const {
  return getPntrToOutput(0)->getForce(kcol*getPntrToOutput(0)->getShape()[1]+jrow);
}


}
}
