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

\par Examples


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
  unsigned getNumberOfDerivatives() override { return 0; }
/// Do the calculation
  void calculate() override;
///
  void apply() override;
};

PLUMED_REGISTER_ACTION(Flatten,"FLATTEN")

void Flatten::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys ); keys.use("ARG");
  keys.setValueDescription("a vector containing all the elements of the input matrix");
}

Flatten::Flatten(const ActionOptions& ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao)
{
  if( getNumberOfArguments()!=1 ) error("should only be one argument for this action");
  if( getPntrToArgument(0)->getRank()!=2 || getPntrToArgument(0)->hasDerivatives() ) error("input to this action should be a matrix");
  getPntrToArgument(0)->buildDataStore(true);
  std::vector<unsigned> inshape( getPntrToArgument(0)->getShape() );
  std::vector<unsigned> shape( 1 ); shape[0]=inshape[0]*inshape[1];
  addValue( shape ); setNotPeriodic(); getPntrToComponent(0)->buildDataStore();
}

void Flatten::calculate() {
  Value* myval = getPntrToComponent(0); unsigned ss=getPntrToArgument(0)->getShape()[1];
  std::vector<double> vals; std::vector<std::pair<unsigned,unsigned> > pairs;
  bool symmetric=getPntrToArgument(0)->isSymmetric();
  unsigned nedge=0; getPntrToArgument(0)->retrieveEdgeList( nedge, pairs, vals );
  for(unsigned l=0; l<nedge; ++l ) {
    unsigned i=pairs[l].first, j=pairs[l].second;
    myval->set( i*ss + j, vals[l] );
    if( symmetric ) myval->set( j*ss + i, vals[l] );
  }
}

void Flatten::apply() {
  if( doNotCalculateDerivatives() || !getPntrToComponent(0)->forcesWereAdded() ) return;

  Value* myval=getPntrToComponent(0); Value* myarg=getPntrToArgument(0);
  unsigned nvals=myval->getNumberOfValues(); for(unsigned j=0; j<nvals; ++j) myarg->addForce( j, myval->getForce(j) );
}

}
}
