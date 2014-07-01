/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014 The plumed team
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
#include "tools/SwitchingFunction.h"
#include "ActionWithVessel.h"

namespace PLMD {
namespace vesselbase{

class LessThan : public FunctionVessel {
private:
  SwitchingFunction sf;
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys ); 
  LessThan( const VesselOptions& da );
  std::string function_description();
  bool calculate();
  void finish();
};

PLUMED_REGISTER_VESSEL(LessThan,"LESS_THAN")

void LessThan::registerKeywords( Keywords& keys ){
  FunctionVessel::registerKeywords( keys );
  SwitchingFunction::registerKeywords( keys );
}

void LessThan::reserveKeyword( Keywords& keys ){
  keys.reserve("numbered","LESS_THAN","calculate the number of variables less than a certain target value. "
                                 "This quantity is calculated using \\f$\\sum_i \\sigma(s_i)\\f$, where \\f$\\sigma(s)\\f$ "
                                 "is a \\ref switchingfunction.",true);
  keys.addOutputComponent("less_than","LESS_THAN","the number of values less than a target value. This is calculated using one of the " 
                                                  "formula described in the description of the keyword so as to make it continuous. "
                                                  "You can calculate this quantity multiple times using different parameters."); 
}

LessThan::LessThan( const VesselOptions& da ) :
FunctionVessel(da)
{
  if( getAction()->isPeriodic() ) error("LESS_THAN is not a meaningful option for periodic variables");
  std::string errormsg; sf.set( getAllInput(), errormsg ); 
  if( errormsg.size()!=0 ) error( errormsg ); 
}

std::string LessThan::function_description(){
  return "the number of values less than " + sf.description();
}

bool LessThan::calculate(){
  double weight=getAction()->getElementValue(1);
  plumed_dbg_assert( weight>=getTolerance() );

  double val=getAction()->getElementValue(0);
  double dval, f = sf.calculate(val, dval); dval*=val;
  double contr=weight*f;
  bool addval=addValueUsingTolerance(0,contr);
  if(addval){ 
    getAction()->chainRuleForElementDerivatives(0, 0, weight*dval, this);
    if(diffweight) getAction()->chainRuleForElementDerivatives(0, 1, f, this);
  }
  return ( contr>getNLTolerance() ); 
}

void LessThan::finish(){
  setOutputValue( getFinalValue(0) ); 
  std::vector<double> df(2); df[0]=1.0; df[1]=0.0;
  mergeFinalDerivatives( df );
}

}
}
