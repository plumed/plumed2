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

mean::mean( const std::string& parameters ) :
DistributionFunction(parameters)
{
  Tools::convert(parameters,nval); 
  addAccumulator( true );
  addAccumulator( false );
}

std::string mean::message(){
  std::ostringstream ostr;
  ostr<<"the average value"; 
  return ostr.str();
}

void mean::printKeywords( Log& log ){
  plumed_massert( 0, "it should be impossible to get here");
}

std::string mean::getLabel(){
  return "average";
}

void mean::calculate( Value* value_in, std::vector<Value>& aux ){
  copyValue( 0, value_in ); 
  setValue( 1, 1.0 );
}

void mean::finish( Value* value_out ){
  if ( getPntrToAccumulator(1)->get()!=nval ) printf("WARNING: A neighbor list is causing discontinuities in an average");
  extractDerivatives( 0, value_out );
  value_out->chainRule(1.0/nval); value_out->set(getPntrToAccumulator(0)->get()/nval); 
}

}
