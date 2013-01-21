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
#include "tools/HistogramBead.h"
#include "ActionWithVessel.h"

namespace PLMD {
namespace vesselbase {

class Between : public FunctionVessel {
private:
  bool norm;
  std::vector<double> df;
  HistogramBead hist;
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  Between( const VesselOptions& da );
  unsigned getNumberOfTerms();
  std::string function_description();
  bool calculate();
  void finish();
};

PLUMED_REGISTER_VESSEL(Between,"BETWEEN")

void Between::registerKeywords( Keywords& keys ){
  FunctionVessel::registerKeywords( keys );
  HistogramBead::registerKeywords( keys );
  keys.addFlag("NORM",false,"calculate the fraction of values rather than the number");
}

void Between::reserveKeyword( Keywords& keys ){
  keys.reserve("numbered","BETWEEN","calculate the number of values that are within a certain range. "
                                    "These quantities are calculated using kernel density estimation as described on "
                                    "\\ref histogrambead.",true); 
}

Between::Between( const VesselOptions& da ) :
FunctionVessel(da)
{ 

  bool isPeriodic=getAction()->isPeriodic();
  double min, max; std::string str_min, str_max;
  if( isPeriodic ){
      getAction()->retrieveDomain( str_min, str_max );
      Tools::convert(str_min,min); Tools::convert(str_max,max);
  }

  parseFlag("NORM",norm); std::string errormsg;
  if(norm){ df.resize(2); df[1]=0.0; }
  else { df.resize(1); df[0]=1.0; }

  hist.set( getAllInput(),"",errormsg );
  if( !isPeriodic ) hist.isNotPeriodic();
  else hist.isPeriodic( min, max ); 
  if( errormsg.size()!=0 ) error( errormsg );
}

std::string Between::function_description(){
  if(norm) return "the fraction of values " + hist.description();
  return "the number of values " + hist.description();
}

unsigned Between::getNumberOfTerms(){
  if(norm) return 2; 
  return 1;
}

bool Between::calculate(){
  double val=getAction()->getElementValue(0);
  double dval, f = hist.calculate(val, dval);
  bool addval=addValue(0,f);
  getAction()->chainRuleForElementDerivatives( 0, 0, dval, this ); 
  if(norm){
     bool ignore=addValue(1,1.0);
     return true;
  } 
  return addval;
}

void Between::finish(){
  if(norm){
     setOutputValue( getFinalValue(0) / getFinalValue(1) ); 
     df[0]=1.0 / getFinalValue(1); 
     mergeFinalDerivatives( df );
  } else {
     setOutputValue( getFinalValue(0) );
     mergeFinalDerivatives( df );
  } 
}

}
}
