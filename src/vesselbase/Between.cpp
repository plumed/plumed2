/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013,2014 The plumed team
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
#include "Between.h"

namespace PLMD {
namespace vesselbase {

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
  keys.addOutputComponent("between","BETWEEN","the number/fraction of values within a certain range. This is calculated using one of the "
                                              "formula described in the description of the keyword so as to make it continuous. "
                                              "You can calculate this quantity multiple times using different parameters."); 
}

Between::Between( const VesselOptions& da ) :
FunctionVessel(da)
{ 
  wnum=getAction()->getIndexOfWeight();
  bool isPeriodic=getAction()->isPeriodic();
  double min, max; std::string str_min, str_max;
  if( isPeriodic ){
      getAction()->retrieveDomain( str_min, str_max );
      Tools::convert(str_min,min); Tools::convert(str_max,max);
  }

  parseFlag("NORM",norm); std::string errormsg; df.resize(2); 

  hist.set( getAllInput(),"",errormsg );
  if( !isPeriodic ) hist.isNotPeriodic();
  else hist.isPeriodic( min, max ); 
  if( errormsg.size()!=0 ) error( errormsg );
}

std::string Between::function_description(){
  if(norm) return "the fraction of values " + hist.description();
  return "the number of values " + hist.description();
}

bool Between::calculate(){
  double weight=getAction()->getElementValue(wnum);
  plumed_dbg_assert( weight>=getTolerance() );
  double val=getAction()->getElementValue(0);
  double dval, f = hist.calculate(val, dval);

  bool bigw=addValueUsingTolerance(1,weight);
  if( !bigw ) return false;

  double contr=weight*f;
  bool addval=addValueUsingTolerance(0,contr);
  if( addval ){
     getAction()->chainRuleForElementDerivatives( 0, 0, weight*dval, this );
     if(diffweight){
        getAction()->chainRuleForElementDerivatives( 0, wnum, f, this ); 
        if(norm) getAction()->chainRuleForElementDerivatives( 1, wnum, 1.0, this );
     }
  }
  return ( contr>getNLTolerance() );
}

void Between::finish(){
  double denom=getFinalValue(1);
  if( norm && diffweight ){ 
     df[0] = 1.0 / denom;
     setOutputValue( getFinalValue(0) / denom ); 
     df[1] = -getFinalValue(0) / ( denom*denom );
     mergeFinalDerivatives( df );
  } else if (norm) {
     df[0] = 1.0 / denom; df[1]=0.0;
     setOutputValue( getFinalValue(0) / denom );
     mergeFinalDerivatives( df );
  } else {
     setOutputValue( getFinalValue(0) );
     df[0] = 1.0; df[1]=0.0;
     mergeFinalDerivatives( df );
  }
}

double Between::getCutoff(){
  return std::numeric_limits<double>::max();
} 

}
}
