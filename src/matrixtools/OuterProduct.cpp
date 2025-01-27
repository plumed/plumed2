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

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace matrixtools {

class OuterProduct : public ActionWithMatrix {
private:
  bool domin, domax, diagzero;
  LeptonCall function;
  unsigned nderivatives;
  bool stored_vector1, stored_vector2;
public:
  static void registerKeywords( Keywords& keys );
  explicit OuterProduct(const ActionOptions&);
  unsigned getNumberOfDerivatives();
  void prepare() override ;
  unsigned getNumberOfColumns() const override {
    return getConstPntrToComponent(0)->getShape()[1];
  }
  void setupForTask( const unsigned& task_index, std::vector<unsigned>& indices, MultiValue& myvals ) const ;
  void performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const override;
  void runEndOfRowJobs( const unsigned& ival, const std::vector<unsigned> & indices, MultiValue& myvals ) const override ;
};

PLUMED_REGISTER_ACTION(OuterProduct,"OUTER_PRODUCT")

void OuterProduct::registerKeywords( Keywords& keys ) {
  ActionWithMatrix::registerKeywords(keys);
  keys.addInputKeyword("compulsory","ARG","vector","the labels of the two vectors from which the outer product is being computed");
  keys.add("compulsory","FUNC","x*y","the function of the input vectors that should be put in the elements of the outer product");
  keys.addFlag("ELEMENTS_ON_DIAGONAL_ARE_ZERO",false,"set all diagonal elements to zero");
  keys.setValueDescription("matrix","a matrix containing the outer product of the two input vectors that was obtained using the function that was input");
}

OuterProduct::OuterProduct(const ActionOptions&ao):
  Action(ao),
  ActionWithMatrix(ao),
  domin(false),
  domax(false) {
  if( getNumberOfArguments()!=2 ) {
    error("should be two arguments to this action, a matrix and a vector");
  }
  if( getPntrToArgument(0)->getRank()!=1 || getPntrToArgument(0)->hasDerivatives() ) {
    error("first argument to this action should be a vector");
  }
  if( getPntrToArgument(1)->getRank()!=1 || getPntrToArgument(1)->hasDerivatives() ) {
    error("first argument to this action should be a vector");
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

  std::vector<unsigned> shape(2);
  shape[0]=getPntrToArgument(0)->getShape()[0];
  shape[1]=getPntrToArgument(1)->getShape()[0];
  addValue( shape );
  setNotPeriodic();
  nderivatives = buildArgumentStore(0);
  std::string headstr=getFirstActionInChain()->getLabel();
  stored_vector1 = getPntrToArgument(0)->ignoreStoredValue( headstr );
  stored_vector2 = getPntrToArgument(1)->ignoreStoredValue( headstr );
  if( getPntrToArgument(0)->isDerivativeZeroWhenValueIsZero() || getPntrToArgument(1)->isDerivativeZeroWhenValueIsZero() ) {
    getPntrToComponent(0)->setDerivativeIsZeroWhenValueIsZero();
  }
}

unsigned OuterProduct::getNumberOfDerivatives() {
  return nderivatives;
}

void OuterProduct::prepare() {
  ActionWithVector::prepare();
  Value* myval=getPntrToComponent(0);
  if( myval->getShape()[0]==getPntrToArgument(0)->getShape()[0] && myval->getShape()[1]==getPntrToArgument(1)->getShape()[0] ) {
    return;
  }
  std::vector<unsigned> shape(2);
  shape[0] = getPntrToArgument(0)->getShape()[0];
  shape[1] = getPntrToArgument(1)->getShape()[0];
  myval->setShape( shape );
}

void OuterProduct::setupForTask( const unsigned& task_index, std::vector<unsigned>& indices, MultiValue& myvals ) const {
  unsigned start_n = getPntrToArgument(0)->getShape()[0], size_v = getPntrToArgument(1)->getShape()[0];
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
  unsigned ostrn = getConstPntrToComponent(0)->getPositionInStream(), ind2=index2;
  if( index2>=getPntrToArgument(0)->getShape()[0] ) {
    ind2 = index2 - getPntrToArgument(0)->getShape()[0];
  }
  if( diagzero && index1==ind2 ) {
    return;
  }

  double fval;
  unsigned jarg = 0, kelem = index1;
  bool jstore=stored_vector1;
  std::vector<double> args(2);
  args[0] = getArgumentElement( 0, index1, myvals );
  args[1] = getArgumentElement( 1, ind2, myvals );
  if( domin ) {
    fval=args[0];
    if( args[1]<args[0] ) {
      fval=args[1];
      jarg=1;
      kelem=ind2;
      jstore=stored_vector2;
    }
  } else if( domax ) {
    fval=args[0];
    if( args[1]>args[0] ) {
      fval=args[1];
      jarg=1;
      kelem=ind2;
      jstore=stored_vector2;
    }
  } else {
    fval=function.evaluate( args );
  }

  myvals.addValue( ostrn, fval );
  if( doNotCalculateDerivatives() ) {
    return ;
  }

  if( domin || domax ) {
    addDerivativeOnVectorArgument( jstore, 0, jarg, kelem, 1.0, myvals );
  } else {
    addDerivativeOnVectorArgument( stored_vector1, 0, 0, index1, function.evaluateDeriv( 0, args ), myvals );
    addDerivativeOnVectorArgument( stored_vector2, 0, 1, ind2, function.evaluateDeriv( 1, args ), myvals );
  }
  if( doNotCalculateDerivatives() || !matrixChainContinues() ) {
    return ;
  }
  unsigned nmat = getConstPntrToComponent(0)->getPositionInMatrixStash(), nmat_ind = myvals.getNumberOfMatrixRowDerivatives( nmat );
  myvals.getMatrixRowDerivativeIndices( nmat )[nmat_ind] = arg_deriv_starts[1] + ind2;
  myvals.setNumberOfMatrixRowDerivatives( nmat, nmat_ind+1 );
}

void OuterProduct::runEndOfRowJobs( const unsigned& ival, const std::vector<unsigned> & indices, MultiValue& myvals ) const {
  if( doNotCalculateDerivatives() || !matrixChainContinues() ) {
    return ;
  }
  unsigned nmat = getConstPntrToComponent(0)->getPositionInMatrixStash(), nmat_ind = myvals.getNumberOfMatrixRowDerivatives( nmat );
  myvals.getMatrixRowDerivativeIndices( nmat )[nmat_ind] = ival;
  myvals.setNumberOfMatrixRowDerivatives( nmat, nmat_ind+1 );
}

}
}
