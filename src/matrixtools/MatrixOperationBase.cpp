/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "MatrixOperationBase.h"

namespace PLMD {
namespace matrixtools {

void MatrixOperationBase::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  keys.addInputKeyword("compulsory","ARG","matrix","the input matrix");
  keys.remove("NUMERICAL_DERIVATIVES");
  keys.addDeprecatedKeyword("MATRIX","ARG");
}

MatrixOperationBase::MatrixOperationBase(const ActionOptions&ao):
  Action(ao),
  ActionWithArguments(ao),
  ActionWithValue(ao) {
  if( getNumberOfArguments()==0 ) {
    std::vector<Value*> args;
    parseArgumentList("MATRIX",args);
    requestArguments(args);
  }
  if( getNumberOfArguments()!=1 ) {
    error("should only be one argument to this action");
  }
  if( getPntrToArgument(0)->getRank()!=2 || getPntrToArgument(0)->hasDerivatives() ) {
    if( getName()=="TRANSPOSE" ) {
      if (getPntrToArgument(0)->getRank()!=1 || getPntrToArgument(0)->hasDerivatives() ) {
        error("input to this argument should be a matrix or vector");
      }
    } else {
      error("input to this argument should be a matrix");
    }
  }
}

void MatrixOperationBase::retrieveFullMatrix( Matrix<double>& mymatrix ) {
  if( mymatrix.nrows()!=getPntrToArgument(0)->getShape()[0] || mymatrix.nrows()!=getPntrToArgument(0)->getShape()[1] ) {
    mymatrix.resize( getPntrToArgument(0)->getShape()[0], getPntrToArgument(0)->getShape()[1] );
  }
  unsigned nedge;
  getPntrToArgument(0)->retrieveEdgeList( nedge, MOBpairs, MOBvals  );
  mymatrix=0;
  bool symmetric = getPntrToArgument(0)->isSymmetric();
  for(unsigned i=0; i<nedge; ++i ) {
    mymatrix( MOBpairs[i].first, MOBpairs[i].second ) = MOBvals[i];
    if( symmetric ) {
      mymatrix( MOBpairs[i].second, MOBpairs[i].first ) = MOBvals[i];
    }
  }
}


void MatrixOperationBase::apply() {
  if( doNotCalculateDerivatives() ) {
    return;
  }

  bool forces=false;
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    if( getPntrToComponent(i)->forcesWereAdded() ) {
      forces=true;
      break;
    }
  }
  if( !forces ) {
    return;
  }

  Value* mat = getPntrToArgument(0);
  unsigned ncols=mat->getNumberOfColumns();
  for(unsigned i=0; i<mat->getShape()[0]; ++i) {
    unsigned ncol = mat->getRowLength(i);
    for(unsigned j=0; j<ncol; ++j) {
      mat->addForce( i*ncols+j, getForceOnMatrixElement( i, mat->getRowIndex(i,j) ), false );
    }
  }
}


}
}
