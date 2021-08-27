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

//+PLUMEDOC FUNCTION ORDINAL_WEIGHTS 
/*
This function can be used to find the highest colvar by magnitude in a set.

\par Examples

*/
//+ENDPLUMEDOC


class OrdinalWeights : public Function {
private:
  unsigned nlowest;
  unsigned nhighest;
public:
  explicit OrdinalWeights(const ActionOptions&);
  void calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const ;
  static void registerKeywords(Keywords& keys);
  bool valuesComputedInChain() const { return false; }
  void transformFinalValueAndDerivatives( const std::vector<double>& buf );
};


PLUMED_REGISTER_ACTION(OrdinalWeights,"ORDINAL_WEIGHTS")

void OrdinalWeights::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys); keys.use("ARG"); keys.remove("PERIODIC");
  keys.add("optional","NLOWEST","give weights of one to the indices of the nlowest input cvs");
  keys.add("optional","NHIGHEST","give weights of one to the indices of the nhighest input cvs");
}

OrdinalWeights::OrdinalWeights(const ActionOptions&ao):
  Action(ao),
  Function(ao),
  nlowest(0),
  nhighest(0)
{
  if( numberedkeys ) error("numbered keys cannot be used with this action");

  if( getPntrToArgument(0)->getRank()>0 ) {
    if( getNumberOfArguments()>1 ) error("should only use one non-scalar argument in input for ARG keyword");
  }
  for(unsigned i=0; i<getNumberOfArguments(); ++i) getPntrToArgument(i)->buildDataStore( getLabel() );

  unsigned nargs=0;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    string s; Tools::convert(i+1,s); nargs += getPntrToArgument(i)->getNumberOfValues();
    if(getPntrToArgument(i)->isPeriodic()) error("Cannot sort periodic values (check argument "+s+")");
  }
  parse("NLOWEST",nlowest); parse("NHIGHEST",nhighest);
  if( nlowest>0 ) log.printf("  giving weights of 1 to %u lowest cvs \n", nlowest );
  if( nhighest>0 ) log.printf("  giving weights of 1 to %u highest cvs \n", nhighest );
  std::vector<unsigned> shape(1); shape[0] = nargs;
  ActionWithValue::addValue( shape ); setNotPeriodic(); checkRead();
}

void OrdinalWeights::calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const { }

void OrdinalWeights::transformFinalValueAndDerivatives( const std::vector<double>& buf ) {
  if( !actionInChain() || getNumberOfArguments()>1 ) return;

  // Count number of values
  unsigned nargs = 0;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      nargs += getPntrToArgument(i)->getNumberOfValues();
  }
  // Retrieve all the values
  std::vector<std::pair<double,int> > argvals( nargs ); unsigned pves = 0;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    Value* myarg = getPntrToArgument(i);
    for(unsigned j=0; j<myarg->getNumberOfValues(); ++j ) { 
        argvals[pves+j].first = myarg->get(j); argvals[pves+j].second = pves+j;
    } 
    pves += myarg->getNumberOfValues();
  }
  // Sort the values
  std::sort( argvals.begin(), argvals.end() ); Value* val0 = getPntrToComponent(0);
  // And set the weights
  for(unsigned i=0;i<nlowest;++i) val0->set( argvals[i].second, 1.0 );
  for(unsigned i=0;i<nhighest;++i) val0->set( argvals[nargs-1-i].second, 1.0 );
}

}
}


