/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "VesselRegister.h"
#include "FunctionVessel.h"
#include "ActionWithVessel.h"

namespace PLMD {
namespace vesselbase{

class Max : public FunctionVessel {
private:
  std::vector<double> df;
  double beta;
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  Max( const VesselOptions& da );
  std::string function_description();
  bool calculate();
  void finish();
};

PLUMED_REGISTER_VESSEL(Max,"MAX")

void Max::registerKeywords( Keywords& keys ){
  FunctionVessel::registerKeywords( keys );
  keys.add("compulsory","BETA","the value of beta for the equation in the manual");
}

void Max::reserveKeyword( Keywords& keys ){
  keys.reserve("optional","MAX","calculate the maximum value. "
                                "To make this quantity continuous the maximum is calculated using "
                                "\\f$ \\textrm{max} = \\beta \\log \\sum_i \\exp\\left( \\frac{s_i}{\\beta}\\right) } \\f$ "
                                "The value of \\f$\\beta\\f$ in this function is specified using (BETA=\\f$\\beta\\f$)",true);
}

Max::Max( const VesselOptions& da ) :
FunctionVessel(da),
df(2)
{
  if( getAction()->isPeriodic() ) error("max is not a meaningful option for periodic variables");
  parse("BETA",beta);

  if( diffweight ) error("can't calculate max if weight is differentiable");
}

std::string Max::function_description(){
  std::string str_beta; Tools::convert( beta, str_beta );
  return "the maximum value. Beta is equal to " + str_beta;
}

bool Max::calculate(){
  double val=getAction()->getElementValue(0);
  double dval, f = exp(val/beta); dval=f/beta;
  addValueIgnoringTolerance(0,f);
  getAction()->chainRuleForElementDerivatives( 0, 0, dval, this );
  return true;
}

void Max::finish(){
  double valin=getFinalValue(0); double dist=beta*std::log( valin );
  setOutputValue( dist ); df[0]=beta/valin; df[1]=0.0;  
  mergeFinalDerivatives( df );
}

}
}
