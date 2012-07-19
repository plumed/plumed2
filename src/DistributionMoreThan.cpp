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

std::string more_than::documentation(){
  std::ostringstream ostr;
  ostr<<"This is calculated using \\f$1.0 - s(r)\\f$, where \\f$s(r)\\f$ is a switching function. ";
  ostr<<SwitchingFunction::documentation();
  return ostr.str();
}

more_than::more_than( const std::string& parameters ) :
DistributionFunction(parameters)
{
  std::string errormsg;
  sf.set( parameters, errormsg );
  if( errormsg.size()!=0 ) error( errormsg );
  addAccumulator( true );
}

std::string more_than::message(){
  std::ostringstream ostr;
  ostr<<"number of values more than "<<sf.description();
  return ostr.str();
}

void more_than::printKeywords( Log& log ){
  sf.printKeywords( log );
}

std::string more_than::getLabel(){
  std::string vv;
  Tools::convert( sf.get_r0(), vv );
  return "gt" + vv;
}

void more_than::calculate( Value* value_in, std::vector<Value>& aux ){
  copyValue( 0, value_in );
  double p, df, f; p=value_in->get(); 
  f=1.0 - sf.calculate(p, df); df*=-p;
  chainRule( 0, df ); setValue( 0, f );
}

void more_than::finish( Value* value_out ){
  extractDerivatives( 0, value_out );
  value_out->set( getPntrToAccumulator(0)->get() );
}

}
