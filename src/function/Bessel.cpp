/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#include "Function.h"
#include "ActionRegister.h"
#include "tools/SwitchingFunction.h"

#include <cmath>

using namespace std;

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION BESSEL
/*
Calculate the value of a Bessel function.

\par Examples

*/
//+ENDPLUMEDOC


class Bessel :
  public Function
{
  unsigned order;
public:
  explicit Bessel(const ActionOptions&);
  void calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const ;
  static void registerKeywords(Keywords& keys);
  void turnOnDerivatives() { error("derivatives for BESSEL have not been implemented"); }
};


PLUMED_REGISTER_ACTION(Bessel,"BESSEL")

void Bessel::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys); keys.use("ARG");
  keys.add("compulsory","ORDER","0","the order of Bessel function to use.  Can only be zero at the moment.");
}

Bessel::Bessel(const ActionOptions&ao):
  Action(ao),
  Function(ao)
{
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->isPeriodic() ) error("cannot use this function on periodic functions");
  }
  parse("ORDER",order); log.printf("  computing %dth order bessel function \n", order );
  if( order!=0 ) error("only zero order bessel function is implemented");
  addValueWithDerivatives();
  checkRead();
}

void Bessel::calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const {
  plumed_dbg_assert( args.size()==1 ); 
  if( order==0 ) {
      if (fabs(args[0])<3.75) {
        double y = Tools::fastpow( args[0]/3.75, 2 );
        addValue( 0, 1 + y*(3.5156229 +y*(3.0899424 + y*(1.2067492+y*(0.2659732+y*(0.0360768+y*0.0045813))))), myvals );
        return;
      }
      double ax=fabs(args[0]), y=3.75/ax, bx=std::exp(ax)/sqrt(ax);
      ax=0.39894228+y*(0.01328592+y*(0.00225319+y*(-0.00157565+y*(0.00916281+y*(-0.02057706+y*(0.02635537+y*(-0.01647633+y*0.00392377)))))));
      addValue( 0, ax*bx, myvals ); 
  } else plumed_error();
}

}
}


