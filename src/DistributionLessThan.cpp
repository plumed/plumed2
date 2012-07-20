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
#include "DistributionFunctions.h"

namespace PLMD {

std::string less_than::documentation(){
  std::ostringstream ostr;
  ostr<<SwitchingFunction::documentation();
  return ostr.str();
}

less_than::less_than( const std::string& parameters ) :
DistributionFunction(parameters)
{
  std::string errormsg;
  sf.set( parameters, errormsg ); 
  if( errormsg.size()!=0 ) error( errormsg ); 
  addAccumulator( true );
}

std::string less_than::message(){
  std::ostringstream ostr;
  ostr<<"number of values less than "<<sf.description();
  return ostr.str();
}

void less_than::printKeywords( Log& log ){
  sf.printKeywords( log );
}

std::string less_than::getLabel(){
  std::string vv;
  Tools::convert( sf.get_r0(), vv );
  return "lt" + vv;
}

void less_than::calculate( Value* value_in, std::vector<Value>& aux ){
  copyValue( 0, value_in );
  double p, df, f; p=value_in->get(); 
  f=sf.calculate(p, df); df*=p;
  chainRule( 0, df ); setValue( 0, f );
}

void less_than::finish( Value* value_out ){
  extractDerivatives( 0, value_out );
  value_out->set( getPntrToAccumulator(0)->get() );
}

}
