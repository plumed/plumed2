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
#include "core/ActionWithMatrix.h"
#include "core/ActionRegister.h"

//+PLUMEDOC MCOLVAR VSTACK
/*
Create a matrix by stacking vectors together

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace valtools {

class VStack : public ActionWithMatrix {
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit VStack(const ActionOptions&);
/// Get the number of derivatives
  unsigned getNumberOfDerivatives() override {
    return 0;
  }
///
  void prepare() override ;
///
  unsigned getNumberOfColumns() const override {
    return getNumberOfArguments();
  }
///
  void setupForTask( const unsigned& task_index, std::vector<unsigned>& indices, MultiValue& myvals ) const override ;
///
  int checkTaskIsActive( const unsigned& itask ) const override ;
///
  void performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const override ;
///
  void runEndOfRowJobs( const unsigned& ival, const std::vector<unsigned> & indices, MultiValue& myvals ) const override ;
///
  void getMatrixColumnTitles( std::vector<std::string>& argnames ) const override ;
};

PLUMED_REGISTER_ACTION(VStack,"VSTACK")

void VStack::registerKeywords( Keywords& keys ) {
  ActionWithMatrix::registerKeywords( keys );
  keys.addInputKeyword("compulsory","ARG","scalar/vector","the values that you would like to stack together to construct the output matrix");
  keys.setValueDescription("matrix","a matrix that contains the input vectors in its columns");
}

VStack::VStack(const ActionOptions& ao):
  Action(ao),
  ActionWithMatrix(ao) {
  if( getNumberOfArguments()==0 ) {
    error("no arguments were specificed");
  }
  if( getPntrToArgument(0)->getRank()>1 ) {
    error("all arguments should be vectors");
  }
  unsigned nvals=1;
  bool periodic=false;
  std::string smin, smax;
  if( getPntrToArgument(0)->getRank()==1 ) {
    nvals = getPntrToArgument(0)->getShape()[0];
  }
  if( getPntrToArgument(0)->isPeriodic() ) {
    periodic=true;
    getPntrToArgument(0)->getDomain( smin, smax );
  }

  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->getRank()>1 || (getPntrToArgument(i)->getRank()==1 && getPntrToArgument(i)->hasDerivatives()) ) {
      error("all arguments should be vectors");
    }
    if( getPntrToArgument(i)->getRank()==0 ) {
      if( nvals!=1 ) {
        error("all input vector should have same number of elements");
      }
    } else if( getPntrToArgument(i)->getShape()[0]!=nvals ) {
      error("all input vector should have same number of elements");
    }
    if( periodic ) {
      if( !getPntrToArgument(i)->isPeriodic() ) {
        error("one argument is periodic but " + getPntrToArgument(i)->getName() + " is not periodic");
      }
      std::string tmin, tmax;
      getPntrToArgument(i)->getDomain( tmin, tmax );
      if( tmin!=smin || tmax!=smax ) {
        error("domain of argument " + getPntrToArgument(i)->getName() + " is different from domain for all other arguments");
      }
    } else if( getPntrToArgument(i)->isPeriodic() ) {
      error("one argument is not periodic but " + getPntrToArgument(i)->getName() + " is periodic");
    }
  }
  // And create a value to hold the matrix
  std::vector<unsigned> shape(2);
  shape[0]=nvals;
  shape[1]=getNumberOfArguments();
  addValue( shape );
  if( periodic ) {
    setPeriodic( smin, smax );
  } else {
    setNotPeriodic();
  }
  // And store this value
  getPntrToComponent(0)->reshapeMatrixStore( shape[1] );
}

void VStack::getMatrixColumnTitles( std::vector<std::string>& argnames ) const {
  for(unsigned j=0; j<getNumberOfArguments(); ++j) {
    if( (getPntrToArgument(j)->getPntrToAction())->getName()=="COLLECT" ) {
      ActionWithArguments* aa = dynamic_cast<ActionWithArguments*>( getPntrToArgument(j)->getPntrToAction() );
      plumed_assert( aa && aa->getNumberOfArguments()==1 );
      argnames.push_back( (aa->getPntrToArgument(0))->getName() );
    } else {
      argnames.push_back( getPntrToArgument(j)->getName() );
    }
  }
}

void VStack::prepare() {
  ActionWithVector::prepare();
  if( getPntrToArgument(0)->getRank()==0 || getPntrToArgument(0)->getShape()[0]==getPntrToComponent(0)->getShape()[0] ) {
    return ;
  }
  std::vector<unsigned> shape(2);
  shape[0] = getPntrToArgument(0)->getShape()[0];
  shape[1] = getNumberOfArguments();
  getPntrToComponent(0)->setShape(shape);
  getPntrToComponent(0)->reshapeMatrixStore( shape[1] );
}

void VStack::setupForTask( const unsigned& task_index, std::vector<unsigned>& indices, MultiValue& myvals ) const {
  unsigned nargs = getNumberOfArguments();
  unsigned nvals = getConstPntrToComponent(0)->getShape()[0];
  if( indices.size()!=nargs+1 ) {
    indices.resize( nargs+1 );
  }
  for(unsigned i=0; i<nargs; ++i) {
    indices[i+1] = nvals + i;
  }
  myvals.setSplitIndex( nargs + 1 );
}

int VStack::checkTaskIsActive( const unsigned& itask ) const {
  return 1;
}

void VStack::performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const {
  unsigned ind2 = index2;
  if( index2>=getConstPntrToComponent(0)->getShape()[0] ) {
    ind2 = index2 - getConstPntrToComponent(0)->getShape()[0];
  }
  myvals.addValue( 0, getPntrToArgument(ind2)->get( index1 ) );

  if( doNotCalculateDerivatives() ) {
    return;
  }
  unsigned vstart=0;
  for(unsigned i=0; i<ind2; ++i) {
    vstart += getPntrToArgument(i)->getNumberOfStoredValues();
  }
  myvals.addDerivative( 0, vstart + index1, 1.0 );
  myvals.updateIndex( 0, vstart + index1 );
}

void VStack::runEndOfRowJobs( const unsigned& ival, const std::vector<unsigned> & indices, MultiValue& myvals ) const {
  if( doNotCalculateDerivatives() ) {
    return ;
  }

  unsigned nmat_ind = myvals.getNumberOfMatrixRowDerivatives();
  std::vector<unsigned>& matrix_indices( myvals.getMatrixRowDerivativeIndices() );
  plumed_assert( nmat_ind<matrix_indices.size() );
  unsigned ncols = getConstPntrToComponent(0)->getShape()[0];
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    matrix_indices[i] = i*ncols + ival;
  }
  myvals.setNumberOfMatrixRowDerivatives( getNumberOfArguments() );
}

}
}
