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
#include "SwitchingFunction.h"
#include "ActionWithDistribution.h"

namespace PLMD {

class VesselLessThan : public SumVessel {
private:
  SwitchingFunction sf;
public:
  static void reserveKeyword( Keywords& keys ); 
  VesselLessThan( const VesselOptions& da );
  double compute( const unsigned& i, const double& val, double& df ); 
  void printKeywords( Log& log );
};

PLUMED_REGISTER_VESSEL(VesselLessThan,"LESS_THAN")

void VesselLessThan::reserveKeyword( Keywords& keys ){
  keys.reserve("optional","LESS_THAN","calculate the number of variables less than a certain target value. "
                                      "This quantity is calculated using \\f$\\sum_i \\sigma(s_i)\\f$, where \\f$\\sigma(s)\\f$ "
                                      "is a \\ref switchingfunction. The final value can be referenced using "
                                      "\\e label.lt\\f$r_0\\f$.");  
}

VesselLessThan::VesselLessThan( const VesselOptions& da ) :
SumVessel(da)
{
  if( getAction()->isPeriodic() ) error("less than is not a meaningful option for periodic variables");
 
  std::string errormsg; sf.set( da.parameters, errormsg ); 
  if( errormsg.size()!=0 ) error( errormsg ); 
  std::string vv; Tools::convert( sf.get_r0(), vv );
  addOutput("lt" + vv);
  log.printf("  value %s.lt%s contains number of values less than %s\n",(getAction()->getLabel()).c_str(),vv.c_str(),(sf.description()).c_str() );
}

void VesselLessThan::printKeywords( Log& log ){
  sf.printKeywords( log );
}

double VesselLessThan::compute( const unsigned& i, const double& val, double& df ){
  plumed_assert( i==0 );
  double f; f = sf.calculate(val, df); df*=val;
  return f;
}

}
