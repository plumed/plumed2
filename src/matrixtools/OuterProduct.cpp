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
#include "core/ActionWithMatrix.h"
#include "core/ActionRegister.h"
#include "tools/LeptonCall.h"

//+PLUMEDOC COLVAR OUTER_PRODUCT
/*
Calculate the outer product matrix of two vectors

This action can be used to calculate the [outer product](https://en.wikipedia.org/wiki/Outer_product) of two
vectors.  As a (useless) example of what can be done with this action consider the following simple input:

```plumed
d1: DISTANCE ATOMS1=1,2 ATOMS2=3,4
d2: DISTANCE ATOMS1=5,6 ATOMS2=7,8 ATOMS3=9,10
pp: OUTER_PRODUCT ARG=d1,d2
PRINT ARG=pp FILE=colvar
```

This input outputs a $2 \times 3$ matrix. If we call the 2 dimensional vector output by the first DISTANCE action
$d$ and the 3 dimensional vector output by the second DISTANCE action $h$ then the $(i,j)$ element of the matrix
output by the action with the label `pp` is given by:

$$
p_{ij} = d_i h_j
$$

These outer product matrices are useful if you are trying to calculate an adjacency matrix that says atoms are
connected if they are within a certain distance of each other and if they satisfy a certain criterion.  For example,
consider the following input:

```plumed
# Determine if atoms are within 5.3 nm of each other
c1: CONTACT_MATRIX GROUP=1-100 SWITCH={GAUSSIAN D_0=5.29 R_0=0.01 D_MAX=5.3}
# Calculate the coordination numbers
ones: ONES SIZE=100
cc: MATRIX_VECTOR_PRODUCT ARG=c1,ones
# Now use MORE_THAN to work out which atoms have a coordination number that is bigger than six
cf: MORE_THAN ARG=cc SWITCH={RATIONAL D_0=5.5 R_0=0.5}
# Now recalculate the contact matrix above as first step towards calculating adjacency matrix that measures if
# atoms are close to each other and both have a coordination number that is bigger than six
c2: CONTACT_MATRIX GROUP=1-100 SWITCH={GAUSSIAN D_0=5.29 R_0=0.01 D_MAX=5.3}
# Now make a matrix in which element i,j is one if atom i and atom j both have a coordination number that is greater than 6
cfm: OUTER_PRODUCT ARG=cf,cf
# And multiply this by our contact matrix to determine the desired adjacency matrix
m: CUSTOM ARG=c2,cfm FUNC=x*y PERIODIC=NO
PRINT ARG=m FILE=colvar
```

This input calculates a adjacency matrix which has element $(i,j)$ equal to one if atoms $i$ and $j$ have coordination numbers
that are greater than 6 and if they are within 5.3 nm of each other.

Notice that you can specify the function of the two input vectors that is to be calculated by using the `FUNC` keyword which accepts
mathematical expressions of $x$ and $y$.  In other words, the elements of the outer product are calculated using the lepton library
that is used in the [CUSTOM](CUSTOM.md) action.  In addition, you can set `FUNC=min` or `FUNC=max` to set the elements of the outer product equal to
the minimum of the two input variables or the maximum respectively.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace matrixtools {

class OuterProduct : public ActionWithMatrix {
private:
  bool domin, domax, diagzero;
  LeptonCall function;
public:
  static void registerKeywords( Keywords& keys );
  explicit OuterProduct(const ActionOptions&);
  unsigned getNumberOfDerivatives();
  void prepare() override ;
  unsigned getNumberOfColumns() const override ;
  void setupForTask( const unsigned& task_index, std::vector<unsigned>& indices, MultiValue& myvals ) const ;
  void performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const override;
  void runEndOfRowJobs( const unsigned& ival, const std::vector<unsigned> & indices, MultiValue& myvals ) const override ;
};

PLUMED_REGISTER_ACTION(OuterProduct,"OUTER_PRODUCT")

void OuterProduct::registerKeywords( Keywords& keys ) {
  ActionWithMatrix::registerKeywords(keys);
  keys.addInputKeyword("compulsory","ARG","vector","the labels of the two vectors from which the outer product is being computed");
  keys.addInputKeyword("optional","MASK","matrix","a matrix that is used to used to determine which elements of the output matrix to compute");
  keys.add("compulsory","FUNC","x*y","the function of the input vectors that should be put in the elements of the outer product");
  keys.addFlag("ELEMENTS_ON_DIAGONAL_ARE_ZERO",false,"set all diagonal elements to zero");
  keys.setValueDescription("matrix","a matrix containing the outer product of the two input vectors that was obtained using the function that was input");
}

OuterProduct::OuterProduct(const ActionOptions&ao):
  Action(ao),
  ActionWithMatrix(ao),
  domin(false),
  domax(false) {
  unsigned nargs=getNumberOfArguments();
  if( getNumberOfMasks()>0 ) {
    nargs = nargs - getNumberOfMasks();
  }
  if( nargs!=2 ) {
    error("should be two arguments to this action, they should both be vectors");
  }
  if( getPntrToArgument(0)->getRank()!=1 || getPntrToArgument(0)->hasDerivatives() ) {
    error("first argument to this action should be a vector");
  }
  if( getPntrToArgument(1)->getRank()!=1 || getPntrToArgument(1)->hasDerivatives() ) {
    error("first argument to this action should be a vector");
  }
  if( getNumberOfMasks()==1 ) {
    if( getPntrToArgument(2)->getRank()!=2 || getPntrToArgument(2)->hasDerivatives() ) {
      error("mask argument should be a matrix");
    }
    if( getPntrToArgument(2)->getShape()[0]!=getPntrToArgument(0)->getShape()[0] ) {
      error("mask argument has wrong size");
    }
    if( getPntrToArgument(2)->getShape()[1]!=getPntrToArgument(1)->getShape()[0] ) {
      error("mask argument has wrong size");
    }
  }

  std::string func;
  parse("FUNC",func);
  if( func=="min") {
    domin=true;
    log.printf("  taking minimum of two input vectors \n");
  } else if( func=="max" ) {
    domax=true;
    log.printf("  taking maximum of two input vectors \n");
  } else {
    log.printf("  with function : %s \n", func.c_str() );
    std::vector<std::string> var(2);
    var[0]="x";
    var[1]="y";
    function.set( func, var, this );
  }
  parseFlag("ELEMENTS_ON_DIAGONAL_ARE_ZERO",diagzero);
  if( diagzero ) {
    log.printf("  setting diagonal elements equal to zero\n");
  }

  std::vector<std::size_t> shape(2);
  shape[0]=getPntrToArgument(0)->getShape()[0];
  shape[1]=getPntrToArgument(1)->getShape()[0];
  addValue( shape );
  setNotPeriodic();
  if( getPntrToArgument(0)->isDerivativeZeroWhenValueIsZero() || getPntrToArgument(1)->isDerivativeZeroWhenValueIsZero() ) {
    getPntrToComponent(0)->setDerivativeIsZeroWhenValueIsZero();
  }
}

unsigned OuterProduct::getNumberOfDerivatives() {
  return getPntrToArgument(0)->getNumberOfStoredValues() + getPntrToArgument(1)->getNumberOfStoredValues();
}

unsigned OuterProduct::getNumberOfColumns() const {
  if( getNumberOfMasks()>0 ) {
    return getPntrToArgument(2)->getNumberOfColumns();
  }
  return getConstPntrToComponent(0)->getShape()[1];
}

void OuterProduct::prepare() {
  ActionWithVector::prepare();
  Value* myval=getPntrToComponent(0);
  if( myval->getShape()[0]==getPntrToArgument(0)->getShape()[0] && myval->getShape()[1]==getPntrToArgument(1)->getShape()[0] ) {
    return;
  }
  std::vector<std::size_t> shape(2);
  shape[0] = getPntrToArgument(0)->getShape()[0];
  shape[1] = getPntrToArgument(1)->getShape()[0];
  myval->setShape( shape );
}

void OuterProduct::setupForTask( const unsigned& task_index, std::vector<unsigned>& indices, MultiValue& myvals ) const {
  unsigned start_n = getPntrToArgument(0)->getShape()[0];
  if( getNumberOfMasks()>0 ) {
    Value* maskarg = getPntrToArgument(2);
    unsigned size_v = maskarg->getRowLength(task_index);
    if( indices.size()!=size_v+1 ) {
      indices.resize( size_v+1 );
    }
    for(unsigned i=0; i<size_v; ++i) {
      indices[i+1] = start_n + maskarg->getRowIndex( task_index, i );
    }
    myvals.setSplitIndex( 1 + size_v );
    return;
  }

  unsigned size_v = getPntrToArgument(1)->getShape()[0];
  if( diagzero ) {
    if( indices.size()!=size_v ) {
      indices.resize( size_v );
    }
    unsigned k=1;
    for(unsigned i=0; i<size_v; ++i) {
      if( task_index==i ) {
        continue ;
      }
      indices[k] = size_v + i;
      k++;
    }
    myvals.setSplitIndex( size_v );
  } else {
    if( indices.size()!=size_v+1 ) {
      indices.resize( size_v+1 );
    }
    for(unsigned i=0; i<size_v; ++i) {
      indices[i+1] = start_n + i;
    }
    myvals.setSplitIndex( size_v + 1 );
  }
}

void OuterProduct::performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const {
  unsigned ind2=index2;
  if( index2>=getPntrToArgument(0)->getShape()[0] ) {
    ind2 = index2 - getPntrToArgument(0)->getShape()[0];
  }
  if( diagzero && index1==ind2 ) {
    return;
  }

  double fval;
  unsigned jarg = 0, kelem = index1;
  std::vector<double> args(2);
  args[0] = getPntrToArgument(0)->get( index1 );
  args[1] = getPntrToArgument(1)->get( ind2 );
  if( domin ) {
    fval=args[0];
    if( args[1]<args[0] ) {
      fval=args[1];
      jarg=getPntrToArgument(0)->getNumberOfStoredValues();
      kelem=ind2;
    }
  } else if( domax ) {
    fval=args[0];
    if( args[1]>args[0] ) {
      fval=args[1];
      jarg=getPntrToArgument(0)->getNumberOfStoredValues();
      kelem=ind2;
    }
  } else {
    fval=function.evaluate( args );
  }

  myvals.addValue( 0, fval );
  if( doNotCalculateDerivatives() ) {
    return ;
  }

  if( domin || domax ) {
    myvals.addDerivative( 0, jarg + kelem, 1.0 );
    myvals.updateIndex( 0, jarg + kelem );
  } else {
    myvals.addDerivative( 0, index1, function.evaluateDeriv( 0, args ) );
    myvals.updateIndex( 0, index1 );
    myvals.addDerivative( 0, getPntrToArgument(0)->getNumberOfStoredValues() + ind2, function.evaluateDeriv( 1, args ) );
    myvals.updateIndex( 0, getPntrToArgument(0)->getNumberOfStoredValues() + ind2 );
  }
  if( doNotCalculateDerivatives() ) {
    return ;
  }
  unsigned nmat_ind = myvals.getNumberOfMatrixRowDerivatives();
  myvals.getMatrixRowDerivativeIndices()[nmat_ind] = getPntrToArgument(0)->getNumberOfStoredValues() + ind2;
  myvals.setNumberOfMatrixRowDerivatives( nmat_ind+1 );
}

void OuterProduct::runEndOfRowJobs( const unsigned& ival, const std::vector<unsigned> & indices, MultiValue& myvals ) const {
  if( doNotCalculateDerivatives() ) {
    return ;
  }
  unsigned nmat_ind = myvals.getNumberOfMatrixRowDerivatives();
  myvals.getMatrixRowDerivativeIndices()[nmat_ind] = ival;
  myvals.setNumberOfMatrixRowDerivatives( nmat_ind+1 );
}

}
}
