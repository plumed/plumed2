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
#include "FunctionVessel.h"
#include "VesselRegister.h"

namespace PLMD {
namespace vesselbase{

class Sum : public FunctionVessel {
private:
  std::vector<double> df;
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  Sum( const VesselOptions& da );
  std::string function_description();
  bool calculate();
  void finish();
};

PLUMED_REGISTER_VESSEL(Sum,"SUM")

void Sum::registerKeywords( Keywords& keys ){
  FunctionVessel::registerKeywords( keys );
}

void Sum::reserveKeyword( Keywords& keys ){
  keys.reserveFlag("SUM",false,"calculate the sum of all the quantities.",true);
}

Sum::Sum( const VesselOptions& da ) :
FunctionVessel(da),
df(2)
{
}

std::string Sum::function_description(){
  return "the sum of all the values"; 
}

bool Sum::calculate(){
  double weight=getAction()->getElementValue(1);
  if( weight<getTolerance() ) return false;

  double val=getAction()->getElementValue(0);
  bool addval=addValue(0,weight*val);
  if(addval) getAction()->chainRuleForElementDerivatives( 0, 0, weight, this );
  if(diffweight) getAction()->chainRuleForElementDerivatives(0, 1, val, this);
  return addval;
}

void Sum::finish(){
  setOutputValue( getFinalValue(0) ); 
  df[0]=1.0; df[1]=0.0;
  mergeFinalDerivatives( df );
}



}
}
