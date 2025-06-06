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

//+PLUMEDOC FUNCTION FLATTEN
/*
Convert a matrix into a vector

This action can be used to convert an input matrix to a vector. The input matrix is flattened
in row-major order so if the input matrix is $N \times M$ the output vector has $N\times M$ elements.
The first $M$ of these elements contain the first row of the input matrix, the second $M$ contain the
second row of the input matrix and so on.

The following example illustrates how we can convert a $5\times 5$ contact matrix into a vector with
25 elements:

```plumed
c1: CONTACT_MATRIX GROUP=1,2,3,4,5 SWITCH={RATIONAL R_0=0.1}
f: FLATTEN ARG=c1
PRINT ARG=f FILE=colvar
```


*/
//+ENDPLUMEDOC

namespace PLMD {
namespace valtools {

class Flatten :
  public ActionWithValue,
  public ActionWithArguments {
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit Flatten(const ActionOptions&);
/// Get the number of derivatives
  unsigned getNumberOfDerivatives() override {
    return 0;
  }
/// Do the calculation
  void calculate() override;
///
  void apply() override;
};

PLUMED_REGISTER_ACTION(Flatten,"FLATTEN")

void Flatten::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  keys.addInputKeyword("compulsory","ARG","matrix","the label for the matrix that you would like to flatten to a vector");
  keys.add("hidden","MASKED_INPUT_ALLOWED","turns on that you are allowed to use masked inputs ");
  keys.setValueDescription("vector","a vector containing all the elements of the input matrix");
  keys.remove("NUMERICAL_DERIVATIVES");
}

Flatten::Flatten(const ActionOptions& ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao) {
  if( getNumberOfArguments()!=1 ) {
    error("should only be one argument for this action");
  }
  if( getPntrToArgument(0)->getRank()!=2 || getPntrToArgument(0)->hasDerivatives() ) {
    error("input to this action should be a matrix");
  }
  std::vector<std::size_t> inshape( getPntrToArgument(0)->getShape() );
  std::vector<std::size_t> shape( 1 );
  shape[0]=inshape[0]*inshape[1];
  addValue( shape );
  setNotPeriodic();
}

void Flatten::calculate() {
  Value* myval = getPntrToComponent(0);
  unsigned ss=getPntrToArgument(0)->getShape()[1];
  std::vector<double> vals;
  std::vector<std::pair<unsigned,unsigned> > pairs;
  bool symmetric=getPntrToArgument(0)->isSymmetric();
  unsigned nedge=0;
  getPntrToArgument(0)->retrieveEdgeList( nedge, pairs, vals );
  for(unsigned l=0; l<nedge; ++l ) {
    unsigned i=pairs[l].first, j=pairs[l].second;
    myval->set( i*ss + j, vals[l] );
    if( symmetric ) {
      myval->set( j*ss + i, vals[l] );
    }
  }
}

void Flatten::apply() {
  if( doNotCalculateDerivatives() || !getPntrToComponent(0)->forcesWereAdded() ) {
    return;
  }

  Value* myval=getPntrToComponent(0);
  Value* myarg=getPntrToArgument(0);
  unsigned nvals=myval->getNumberOfValues();
  for(unsigned j=0; j<nvals; ++j) {
    myarg->addForce( j, myval->getForce(j) );
  }
}

}
}
