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

This action takes a matrix in input and calculates the input matrix's [tranpose](https://en.wikipedia.org/wiki/Transpose).
The following example shows how you can use this to calculate coordination numbers of species A with species B and vice
versa.

```plumed
# Calculate the contact matrix between the two groups
c1: CONTACT_MATRIX GROUPA=1-10 GROUPB=11-30 SWITCH={RATIONAL R_0=0.1}
# Calculate the cooordination numbers for the atoms in group A by multiplying by a vector of ones
onesB: ONES SIZE=20
coordA: MATRIX_VECTOR_PRODUCT ARG=c1,onesB
# Transpose the contact matrix
c1T: TRANSPOSE ARG=c1
# Calculate the coordination number for the atoms in group B by multiplying the transpose by a vector of ones
onesA: ONES SIZE=10
coordB: MATRIX_VECTOR_PRODUCT ARG=c1T,onesA
# Output the two vectors of coordination numbers to a file
PRINT ARG=coordA,coordB FILE=colvar
```

Another useful example where the transpose can be used is shown below.  In this input the [DISTANCE](DISTANCE.md) command
is used to calculate the orientation of a collection of molecules.  We then can then use the [VSTACK](VSTACK.md), TRANSPOSE and the
[MATRIX_PRODUCT](MATRIX_PRODUCT.md) commands to calculate the dot products between all these vectors

```plumed
# Calculate the vectors connecting these three pairs of atoms
d: DISTANCE COMPONENTS ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6
# Construct a matrix that contains all the components of the vectors calculated
v: VSTACK ARG=d.x,d.y,d.z
# Transpose v
vT: TRANSPOSE ARG=v
# And now calculate the 3x3 matrix of dot products
m: MATRIX_PRODUCT ARG=v,vT
# And output the matrix product to a file
PRINT ARG=m FILE=colvar
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace matrixtools {

class TransposeMatrix : public MatrixOperationBase {
private: 
/// Holds the lengths of the rows
  std::vector<unsigned> row_lengths;
/// Holds the lengths of the columns
  std::vector<unsigned> column_lengths;
/// The array for matrix bookeeping
  std::vector<unsigned> matrix_bookeeping;
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
  double getForceOnMatrixElement( const unsigned& jrow, const unsigned& krow ) const override { plumed_error(); }
};

PLUMED_REGISTER_ACTION(TransposeMatrix,"TRANSPOSE")

void TransposeMatrix::registerKeywords( Keywords& keys ) {
  MatrixOperationBase::registerKeywords( keys );
  keys.addInputKeyword("compulsory","ARG","vector/matrix","the label of the vector or matrix that should be transposed");
  keys.setValueDescription("vector/matrix","the transpose of the input matrix");
  keys.add("hidden","MASKED_INPUT_ALLOWED","turns on that you are allowed to use masked inputs ");
}

TransposeMatrix::TransposeMatrix(const ActionOptions& ao):
  Action(ao),
  MatrixOperationBase(ao) {
  if( getPntrToArgument(0)->isSymmetric() ) {
    error("input matrix is symmetric.  Transposing will achieve nothing!");
  }
  std::vector<std::size_t> shape;
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
    row_lengths.resize( shape[1] );
    column_lengths.resize( shape[0] );
  }
  addValue( shape );
  if( getPntrToArgument(0)->isPeriodic() ) {
    std::string smin, smax;
    getPntrToArgument(0)->getDomain( smin, smax );
    setPeriodic( smin, smax );
  } else {
    setNotPeriodic();
  }
}

void TransposeMatrix::prepare() {
  Value* myval = getPntrToComponent(0);
  Value* myarg = getPntrToArgument(0);
  if( myarg->getRank()==1 ) {
    if( myval->getShape()[0]!=1 || myval->getShape()[1]!=myarg->getShape()[0] ) {
      std::vector<std::size_t> shape(2);
      shape[0] = 1;
      shape[1] = myarg->getShape()[0];
      myval->setShape( shape );
    } 
    if( myval->getNumberOfColumns()!=myarg->getShape()[0] ) {
      myval->reshapeMatrixStore( myarg->getShape()[0]  );
    }
  } else if( myarg->getShape()[0]==1 ) {
    if( myval->getShape()[0]!=myarg->getShape()[1] ) {
      std::vector<std::size_t> shape(1);
      shape[0] = myarg->getShape()[1];
      myval->setShape( shape );
    }
  } else if( myarg->getShape()[0]!=myval->getShape()[1] || myarg->getShape()[1]!=myval->getShape()[0] ) {
    std::vector<std::size_t> shape(2);
    shape[0] = myarg->getShape()[1];
    shape[1] = myarg->getShape()[0];
    myval->setShape( shape );
    row_lengths.resize( shape[1] );
    column_lengths.resize( shape[0] );
  }
}

void TransposeMatrix::calculate() {
  // Retrieve the non-zero pairs
  Value* myarg=getPntrToArgument(0);
  Value* myval=getPntrToComponent(0);
  if( myarg->getRank()<=1 || myval->getRank()==1 ) {
    unsigned nv=myarg->getNumberOfValues();
    for(unsigned i=0; i<nv; ++i) {
      myval->set( i, myarg->get(i) );
    }
  } else {
    // Find the lengths of all the columns
    std::fill( column_lengths.begin(), column_lengths.end(), 0 );
    for(unsigned i=0; i<myarg->getShape()[0]; ++i) { 
        unsigned nr = myarg->getRowLength(i);
        for(unsigned j=0; j<nr; ++j) {
            unsigned ind = myarg->getRowIndex( i, j );    
            column_lengths[ind]++;
        }
    }
    // Find the longest column
    unsigned maxcol = column_lengths[0];
    for(unsigned i=1; i<column_lengths.size(); ++i) {
        if( column_lengths[i]>maxcol ) {
            maxcol = column_lengths[i];
        }
    }
    myval->reshapeMatrixStore( maxcol );
    matrix_bookeeping.resize( myval->getShape()[0]*(1+maxcol) );
    for(unsigned i=0; i<column_lengths.size(); ++i) {
        matrix_bookeeping[i*(1+maxcol)] = 0;
    }
    // Now get the bookeeping 
    unsigned arg_ncol = myarg->getNumberOfColumns();
    for(unsigned i=0; i<myarg->getShape()[0]; ++i) {
        unsigned nr = myarg->getRowLength(i);
        for(unsigned j=0; j<nr; ++j) {
            unsigned ind = myarg->getRowIndex( i, j );
            unsigned startpos = ind*(1+maxcol); 
            matrix_bookeeping[startpos + 1 + matrix_bookeeping[startpos]] = i; 
            myval->set( ind*maxcol + matrix_bookeeping[startpos], myarg->get( i*arg_ncol+j, false)  );
            matrix_bookeeping[startpos]++;
        }
    }
    // And copy all the bookeeping to the output value
    for(unsigned i=0; i<matrix_bookeeping.size(); ++i) {
        myval->setMatrixBookeepingElement( i, matrix_bookeeping[i] );
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
      unsigned narg_cols = myarg->getNumberOfColumns();
      unsigned nval_cols = myval->getNumberOfColumns();
      std::fill( row_lengths.begin(), row_lengths.end(), 0 );
      for(unsigned i=0; i<myval->getShape()[0]; ++i) {
          unsigned nr = myval->getRowLength(i);
          for(unsigned j=0; j<nr; ++j) {
              unsigned ind = myval->getRowIndex( i, j );
              myarg->addForce( narg_cols*ind + row_lengths[ind], myval->getForce(nval_cols*i + j), false );
              row_lengths[ind]++;       
          }
      }
    }
  }
}

}
}
