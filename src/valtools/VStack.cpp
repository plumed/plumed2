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
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "core/ActionRegister.h"

//+PLUMEDOC MCOLVAR VSTACK
/*
Create a matrix by stacking vectors together

This action takes 2 or more vectors with the same number of elements as input and outputs a matrix.
This output matrix is constructed stacking the input vectors vertically.  The first row of the output
matrix thus contains the first elements of all the input vector, the second row contains the second elements
and so on.  In other words, the first input vector is in the first column of the output matrix, the second
input vector is in the second column and so.

The following shows an example of how this action operates in practise. The [DISTANCE](DISTANCE.md) command below calculates
the vectors containing four pairs of atoms.  The VSTACK command is then used to construct a $4 \times 3$ matrix
that contains all the vectors. The 1,1, 1,2 and 1,3 components in this matrix contain the
$x$, $y$ and $z$ components of the vector connecting atoms 1 and 2. The 2,1, 2,2 and 2,3 contain the
$x$, $y$ and $z$ components of the vector connecting atoms 3 and 4 and so on.

```plumed
d1: DISTANCE ...
   COMPONENTS
   ATOMS1=1,2 ATOMS2=3,4
   ATOMS3=5,6 ATOMS4=7,8
...
m: VSTACK ARG=d1.x,d1.y,d1.z

PRINT ARG=m FILE=matrix.dat
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace valtools {

class VStack :
  public ActionWithValue,
  public ActionWithArguments {
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
  void calculate() override ;
///
  void apply() override ;
///
  void getMatrixColumnTitles( std::vector<std::string>& argnames ) const override ;
};

PLUMED_REGISTER_ACTION(VStack,"VSTACK")

void VStack::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  keys.remove("NUMERICAL_DERIVATIVES");
  keys.addInputKeyword("compulsory","ARG","scalar/vector","the values that you would like to stack together to construct the output matrix");
  keys.setValueDescription("matrix","a matrix that contains the input vectors in its columns");
}

VStack::VStack(const ActionOptions& ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao) {
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

  bool derivbool=true;
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
    if( !getPntrToArgument(i)->isDerivativeZeroWhenValueIsZero() ) {
      derivbool=false;
    }
  }
  // And create a value to hold the matrix
  std::vector<std::size_t> shape(2);
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
  if( derivbool ) {
    getPntrToComponent(0)->setDerivativeIsZeroWhenValueIsZero();
  }
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
  if( getPntrToArgument(0)->getRank()==0 || getPntrToArgument(0)->getShape()[0]==getPntrToComponent(0)->getShape()[0] ) {
    return ;
  }
  std::vector<std::size_t> shape(2);
  shape[0] = getPntrToArgument(0)->getShape()[0];
  shape[1] = getNumberOfArguments();
  getPntrToComponent(0)->setShape(shape);
  getPntrToComponent(0)->reshapeMatrixStore( shape[1] );
}

void VStack::calculate() {
  unsigned nvals=1;
  if( getPntrToArgument(0)->getRank()==1 ) {
    nvals = getPntrToArgument(0)->getShape()[0];
  }

  Value* valout = getPntrToComponent(0);
  unsigned nargs = getNumberOfArguments();
  for(unsigned i=0; i<nargs; ++i) {
    unsigned ipos = i;
    Value* myarg = getPntrToArgument(i);
    for(unsigned j=0; j<nvals; ++j) {
      valout->set( ipos, myarg->get(j) );
      ipos += nargs;
    }
  }
}

void VStack::apply() {
  if( !getPntrToComponent(0)->forcesWereAdded() ) {
    return;
  }

  unsigned nvals=1;
  if( getPntrToArgument(0)->getRank()==1 ) {
    nvals = getPntrToArgument(0)->getShape()[0];
  }
  Value* valout = getPntrToComponent(0);
  unsigned nargs = getNumberOfArguments();
  for(unsigned i=0; i<nargs; ++i) {
    unsigned ipos = i;
    Value* myarg = getPntrToArgument(i);
    for(unsigned j=0; j<nvals; ++j) {
      myarg->addForce( j, valout->getForce( ipos ) );
      ipos += nargs;
    }
  }
}

}
}
