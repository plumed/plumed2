/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#include "ActionRegister.h"
#include "Function.h"

#include <cmath>
#include <algorithm>
#include <utility>

using namespace std;

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION LOWEST
/*
This function can be used to find the lowest colvar by magnitude in a set.

\par Examples

*/
//+ENDPLUMEDOC


class Lowest : public Function {
private:
  MultiValue tvals;
public:
  explicit Lowest(const ActionOptions&);
  void calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const ;
  static void registerKeywords(Keywords& keys);
  void transformFinalValueAndDerivatives( const std::vector<double>& buf );
};


PLUMED_REGISTER_ACTION(Lowest,"LOWEST")

void Lowest::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys); keys.use("ARG");
}

Lowest::Lowest(const ActionOptions&ao):
  Action(ao),
  Function(ao),
  tvals(0,0)
{
  if( !numberedkeys ) {
    if( getPntrToArgument(0)->getRank()>0 ) {
      if( getNumberOfArguments()>1 ) error("should only use one non-scalar argument in input for ARG keyword");
    }
    for(unsigned i=0; i<getNumberOfArguments(); ++i) getPntrToArgument(i)->buildDataStore( getLabel() );
  }
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    string s; Tools::convert(i+1,s);
    if(getPntrToArgument(i)->isPeriodic()) error("Cannot sort periodic values (check argument "+s+")");
  }
  addValueWithDerivatives(); checkRead();
}

void Lowest::calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const {
  if( args.size()>1 ) {
    double lowest = args[0]; unsigned lowind = 0;
    for(unsigned i=1; i<args.size(); ++i) {
      if( args[i]<lowest ) { lowest = args[i]; lowind = i; }
    }
    addValue( 0, lowest, myvals ); addDerivative( 0, lowind, 1.0, myvals );
  }
}

void Lowest::transformFinalValueAndDerivatives( const std::vector<double>& buf ) {
  if( !actionInChain() || getNumberOfArguments()>1 ) return;

  unsigned lind = 0, pves = 0; unsigned aind=0; double lowest = getPntrToArgument(0)->get(0);
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    Value* myarg = getPntrToArgument(i);
    for(unsigned j=0; j<myarg->getNumberOfValues(); ++j ) {
      if( myarg->get(j)<lowest ) { aind=i; lowest=myarg->get(j); lind = pves + j; }
    }
    pves += myarg->getNumberOfValues();
  }
  Value* val0 = getPntrToComponent(0); val0->set( lowest );
  if( !doNotCalculateDerivatives() ) {
    unsigned nn=0, nm=0;
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      nn += getPntrToArgument(i)->getNumberOfValues();
      if( lind<nn ) { break; }
      nm += getPntrToArgument(i)->getNumberOfValues();
    }
    tvals.clearAll(); (getPntrToArgument(aind)->getPntrToAction())->rerunTask( lind - nm, tvals );
    unsigned istrn = getPntrToArgument(aind)->getPositionInStream();
    for(unsigned i=0; i<tvals.getNumberActive(istrn); ++i) {
      unsigned ider = tvals.getActiveIndex(istrn,i); val0->addDerivative( ider, tvals.getDerivative( istrn, ider ) );
    }
  }
}

}
}


