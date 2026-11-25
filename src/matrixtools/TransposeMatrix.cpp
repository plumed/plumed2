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
#include "MatrixOperationBase.h"
#include "core/ActionRegister.h"

//+PLUMEDOC MCOLVAR TRANSPOSE
/*
Calculate the transpose of a matrix

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace matrixtools {

class TransposeMatrix : public MatrixOperationBase {
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit TransposeMatrix(const ActionOptions&);
///
  unsigned getNumberOfDerivatives() override {
    return 0;
  }
///
  void prepare() override ;
///
  void calculate() override ;
///
  void apply() override ;
///
  double getForceOnMatrixElement( const unsigned& jrow, const unsigned& krow ) const override;
};

PLUMED_REGISTER_ACTION(TransposeMatrix,"TRANSPOSE")

void TransposeMatrix::registerKeywords( Keywords& keys ) {
  MatrixOperationBase::registerKeywords( keys );
  keys.setValueDescription("the transpose of the input matrix");
}

TransposeMatrix::TransposeMatrix(const ActionOptions& ao):
  Action(ao),
  MatrixOperationBase(ao) {
  if( getPntrToArgument(0)->isSymmetric() ) {
    error("input matrix is symmetric.  Transposing will achieve nothing!");
  }
  std::vector<unsigned> shape;
  if( getPntrToArgument(0)->getRank()==0 ) {
    error("transposing a scalar?");
  } else if( getPntrToArgument(0)->getRank()==1 ) {
    shape.resize(2);
    shape[0]=1;
    shape[1]=getPntrToArgument(0)->getShape()[0];
  } else if( getPntrToArgument(0)->getShape()[0]==1 ) {
    shape.resize(1);
    shape[0] = getPntrToArgument(0)->getShape()[1];
  } else {
    shape.resize(2);
    shape[0]=getPntrToArgument(0)->getShape()[1];
    shape[1]=getPntrToArgument(0)->getShape()[0];
  }
  addValue( shape );
  setNotPeriodic();
  getPntrToComponent(0)->buildDataStore();
  if( shape.size()==2 ) {
    getPntrToComponent(0)->reshapeMatrixStore( shape[1] );
  }
}

void TransposeMatrix::prepare() {
  Value* myval = getPntrToComponent(0);
  Value* myarg = getPntrToArgument(0);
  if( myarg->getRank()==1 ) {
    if( myval->getShape()[0]!=1 || myval->getShape()[1]!=myarg->getShape()[0] ) {
      std::vector<unsigned> shape(2);
      shape[0] = 1;
      shape[1] = myarg->getShape()[0];
      myval->setShape( shape );
      myval->reshapeMatrixStore( shape[1] );
    }
  } else if( myarg->getShape()[0]==1 ) {
    if( myval->getShape()[0]!=myarg->getShape()[1] ) {
      std::vector<unsigned> shape(1);
      shape[0] = myarg->getShape()[1];
      myval->setShape( shape );
    }
  } else if( myarg->getShape()[0]!=myval->getShape()[1] || myarg->getShape()[1]!=myval->getShape()[0] ) {
    std::vector<unsigned> shape(2);
    shape[0] = myarg->getShape()[1];
    shape[1] = myarg->getShape()[0];
    myval->setShape( shape );
    myval->reshapeMatrixStore( shape[1] );
  }
}

void TransposeMatrix::calculate() {
  // Retrieve the non-zero pairs
  Value* myarg=getPntrToArgument(0);
  Value* myval=getPntrToComponent(0);
  if( myarg->getRank()<=1 || myval->getRank()==1 ) {
    if( myarg->getRank()<=1 && myval->getShape()[1]!=myarg->getShape()[0] ) {
      std::vector<unsigned> shape( 2 );
      shape[0] = 1;
      shape[1] = myarg->getShape()[0];
      myval->setShape( shape );
      myval->reshapeMatrixStore( shape[1] );
    } else if( myval->getRank()==1 && myval->getShape()[0]!=myarg->getShape()[1] ) {
      std::vector<unsigned> shape( 1 );
      shape[0] = myarg->getShape()[1];
      myval->setShape( shape );
    }
    unsigned nv=myarg->getNumberOfValues();
    for(unsigned i=0; i<nv; ++i) {
      myval->set( i, myarg->get(i) );
    }
  } else {
    if( myarg->getShape()[0]!=myval->getShape()[1] || myarg->getShape()[1]!=myval->getShape()[0] ) {
      std::vector<unsigned> shape( 2 );
      shape[0] = myarg->getShape()[1];
      shape[1] = myarg->getShape()[0];
      myval->setShape( shape );
      myval->reshapeMatrixStore( shape[1] );
    }
    std::vector<double> vals;
    std::vector<std::pair<unsigned,unsigned> > pairs;
    std::vector<unsigned> shape( myval->getShape() );
    unsigned nedge=0;
    myarg->retrieveEdgeList( nedge, pairs, vals );
    for(unsigned i=0; i<nedge; ++i) {
      myval->set( pairs[i].second*shape[1] + pairs[i].first, vals[i] );
    }
  }
}

void TransposeMatrix::apply() {
  if( doNotCalculateDerivatives() ) {
    return;
  }

  // Apply force on the matrix
  if( getPntrToComponent(0)->forcesWereAdded() ) {
    Value* myarg=getPntrToArgument(0);
    Value* myval=getPntrToComponent(0);
    if( myarg->getRank()<=1 || myval->getRank()==1 ) {
      unsigned nv=myarg->getNumberOfValues();
      for(unsigned i=0; i<nv; ++i) {
        myarg->addForce( i, myval->getForce(i) );
      }
    } else {
      MatrixOperationBase::apply();
    }
  }
}

double TransposeMatrix::getForceOnMatrixElement( const unsigned& jrow, const unsigned& kcol ) const {
  return getConstPntrToComponent(0)->getForce(kcol*getConstPntrToComponent(0)->getShape()[1]+jrow);
}


}
}
