/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2014 The plumed team
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
namespace vesselbase {

class Mean : public FunctionVessel {
private:
  unsigned wnum;
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  Mean( const vesselbase::VesselOptions& da );
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
  keys.addOutputComponent("mean","MEAN","the mean value. The output component can be refererred to elsewhere in the input "
                                        "file by using the label.mean");
}

Mean::Mean( const vesselbase::VesselOptions& da ) :
FunctionVessel(da)
{
  wnum=getAction()->getIndexOfWeight();
  if( getAction()->isPeriodic() ) error("MEAN cannot be used with periodic variables");
}

std::string Mean::function_description(){
  return "the mean value";
}

bool Mean::calculate(){
  double weight=getAction()->getElementValue(wnum);
  plumed_dbg_assert( weight>=getTolerance() );
  addValueIgnoringTolerance( 1, weight );
  double colvar=getAction()->getElementValue(0);
  addValueIgnoringTolerance( 0, weight*colvar  );
  getAction()->chainRuleForElementDerivatives( 0, 0, weight, this );
  if(diffweight){
     getAction()->chainRuleForElementDerivatives( 0, wnum, colvar, this );
     getAction()->chainRuleForElementDerivatives( 1, wnum, 1.0, this );
  }
  return true;
}

void Mean::finish(){
  setOutputValue( getFinalValue(0) / getFinalValue(1) ); 
  double denom=getFinalValue(1);
  std::vector<double> df(2); 
  df[0] = 1.0 / denom; 
  if(diffweight) df[1] = -getFinalValue(0) / (denom*denom); 
  else df[1]=0.0;
  mergeFinalDerivatives( df );
}

}
}
