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
public:
  explicit Lowest(const ActionOptions&);
  void calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const ;
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(Lowest,"LOWEST")

void Lowest::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys); keys.use("ARG");
}

Lowest::Lowest(const ActionOptions&ao):
  Action(ao),
  Function(ao)
{
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    string s; Tools::convert(i+1,s);
    if(getPntrToArgument(i)->isPeriodic()) error("Cannot sort periodic values (check argument "+s+")");
  }
  addValueWithDerivatives(); checkRead();
}

void Lowest::calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const {
  double lowest = args[0]; unsigned lowind = 0;
  for(unsigned i=1; i<args.size(); ++i) {
    if( args[i]<lowest ){ lowest = args[i]; lowind = 0; } 
  }
  setValue( 0, lowest, myvals ); addDerivative( 0, lowind, 1.0, myvals );
}

}
}


