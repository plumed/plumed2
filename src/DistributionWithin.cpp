/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The PLUMED team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of PLUMED, version 2.0.

   PLUMED is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   PLUMED is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with PLUMED.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "DistributionFunctions.h"

namespace PLMD {

std::string within::documentation(){
  std::ostringstream ostr;
  ostr<<"To make this quantity continuous it is calculated using "<<HistogramBead::documentation(false)<<". ";
  ostr<<"Adding the NORM flag allows you to calculate the fraction of colvars in the particular range rather than the total number ";
  ostr<<"(N.B. this option should probably not used if you are using neighbor lists and the derivatives of the WITHIN value).";
  return ostr.str();
}

within::within( const std::string& parameters ) :
DistributionFunction(parameters),
usenorm(false)
{ 
  std::string errormsg;
  std::vector<std::string> data=Tools::getWords(parameters);
  Tools::parseFlag(data,"NORM",usenorm);
  hist.set( parameters, "", errormsg );
  if( errormsg.size()!=0 ) error( errormsg );
  addAccumulator( true );
  addAccumulator( false );
}

std::string within::message(){
  std::ostringstream ostr;
  if(usenorm) ostr<<"fraction of values "<<hist.description();
  else ostr<<"number of values "<<hist.description(); 
  return ostr.str();
}

void within::printKeywords( Log& log ){
  hist.printKeywords( log );
}

std::string within::getLabel(){
  std::string lb, ub;
  Tools::convert( hist.getlowb(), lb );
  Tools::convert( hist.getbigb(), ub );
  return "between" + lb + "&" + ub;
}

void within::calculate( Value* value_in, std::vector<Value>& aux ){
  copyValue( 0, value_in ); 
  double df, f; f=hist.calculate( value_in->get() , df );
  chainRule(0, df); setValue(0, f);
  setValue( 1, 1.0 );
}

void within::finish( Value* value_out ){
  extractDerivatives( 0, value_out );
  if( usenorm ) {
    double nval=getPntrToAccumulator(1)->get();
    value_out->chainRule(1.0/nval); value_out->set(getPntrToAccumulator(0)->get()/nval);
  } else {
    value_out->set( getPntrToAccumulator(0)->get() );
  }
}

}
