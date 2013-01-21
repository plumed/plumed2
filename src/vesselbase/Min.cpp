/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#include "VesselRegister.h"
#include "FunctionVessel.h"
#include "ActionWithVessel.h"

namespace PLMD {
namespace vesselbase{

class Min : public FunctionVessel {
private:
  std::vector<double> df;
  double beta;
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  Min( const VesselOptions& da );
  unsigned getNumberOfTerms(){ return 1; }
  std::string function_description();
  bool calculate();
  void finish();
};

PLUMED_REGISTER_VESSEL(Min,"MIN")

void Min::registerKeywords( Keywords& keys ){
  FunctionVessel::registerKeywords( keys );
  keys.add("compulsory","BETA","the value of beta for the equation in the manual");
}

void Min::reserveKeyword( Keywords& keys ){
  keys.reserve("optional","MIN","calculate the minimum value. "
                                "To make this quantity continuous the minimum is calculated using "
                                "\\f$ \\textrm{min} = \\frac{\\beta}{ \\log \\sum_i \\exp\\left( \\frac{\\beta}{s_i} \\right) } \\f$ "
                                "The value of \\f$\\beta\\f$ in this function is specified using (BETA=\\f$\\beta\\f$)",true);
}

Min::Min( const VesselOptions& da ) :
FunctionVessel(da),
df(1)
{
  if( getAction()->isPeriodic() ) error("min is not a meaningful option for periodic variables");
  parse("BETA",beta);
}

std::string Min::function_description(){
  std::string str_beta; Tools::convert( beta, str_beta );
  return "the minimum value. Beta is equal to " + str_beta;
}

bool Min::calculate(){
  double val=getAction()->getElementValue(0);
  double dval, f = exp(beta/val); dval=f/(val*val);
  bool addval=addValue(0,f);
  if(addval) getAction()->chainRuleForElementDerivatives( 0, 0, dval, this );
  return addval;
}

void Min::finish(){
  double valin=getFinalValue(0); double dist=beta/std::log( valin );
  setOutputValue( dist ); df[0]=dist*dist/valin;  
  mergeFinalDerivatives( df );
}

}
}
