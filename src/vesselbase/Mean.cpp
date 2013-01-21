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
namespace vesselbase {

class Mean : public FunctionVessel {
private:
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  Mean( const vesselbase::VesselOptions& da );
  unsigned getNumberOfTerms(){ return 2; }
  std::string function_description();
  bool calculate();
  void finish();
};

PLUMED_REGISTER_VESSEL(Mean,"MEAN")

void Mean::registerKeywords( Keywords& keys ){
  FunctionVessel::registerKeywords(keys);
}

void Mean::reserveKeyword( Keywords& keys ){
  keys.reserveFlag("MEAN",false,"take the mean of these variables.",true);
}

Mean::Mean( const vesselbase::VesselOptions& da ) :
FunctionVessel(da)
{
  if( getAction()->isPeriodic() ) error("MEAN cannot be used with periodic variables");
}

std::string Mean::function_description(){
  return "the mean value";
}

bool Mean::calculate(){
  double val=getAction()->getElementValue(0);
  bool ignore=addValue(1,1.0);
  bool addval=addValue( 0, getAction()->getElementValue(0) );
  getAction()->chainRuleForElementDerivatives( 0, 0, 1.0, this );
  return true;
}

void Mean::finish(){
  setOutputValue( getFinalValue(0) / getFinalValue(1) ); 
  std::vector<double> df(2); df[0]=1.0 / getFinalValue(1); df[1]=0.0; 
  mergeFinalDerivatives( df );
}

}
}
