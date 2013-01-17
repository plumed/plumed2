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
#include "SumVessel.h"
#include "tools/SwitchingFunction.h"
#include "ActionWithVessel.h"

namespace PLMD {
namespace vesselbase{

class VesselMoreThan : public SumVessel {
private:
  SwitchingFunction sf;
public:
  static void reserveKeyword( Keywords& keys );
  VesselMoreThan( const VesselOptions& da );
  double compute( const unsigned& i, const double& val, double& df ); 
  void printKeywords();
};

PLUMED_REGISTER_VESSEL(VesselMoreThan,"MORE_THAN")

void VesselMoreThan::reserveKeyword( Keywords& keys ){
  keys.reserve("numbered","MORE_THAN","calculate the number of variables more than a certain target value. "
                                      "This quantity is calculated using \\f$\\sum_i 1.0 - \\sigma(s_i)\\f$, where \\f$\\sigma(s)\\f$ "
                                      "is a \\ref switchingfunction. The final value can be referenced using "
                                      "\\e label.gt\\f$(d_0 + r_0)\\f$.");
}

VesselMoreThan::VesselMoreThan( const VesselOptions& da ) :
SumVessel(da)
{
  if( getAction()->isPeriodic() ) error("more than is not a meaningful option for periodic variables");

  std::string errormsg;
  sf.set( da.parameters, errormsg );
  if( errormsg.size()!=0 ) error( errormsg );
  std::string vv; Tools::convert( sf.get_d0() + sf.get_r0(), vv );
  ActionWithValue* aval=dynamic_cast<ActionWithValue*>( getAction() ); plumed_assert( aval );
  if( aval->exists("gt"+vv) ) error("Either d_0 or r_0 must be different if there are multiple instances of the MORE_THAN keyword on one line");
  addOutput("gt" + vv);
  log.printf("  value %s.gt%s contains the number of values more than %s\n",(getAction()->getLabel()).c_str(),vv.c_str(),(sf.description()).c_str());
}

void VesselMoreThan::printKeywords(){
  sf.printKeywords( log );
}

double VesselMoreThan::compute( const unsigned& i, const double& val, double& df ){
  plumed_dbg_assert( i==0 );
  double f; f=1.0 - sf.calculate(val, df); df*=-val;
  return f;
}

}
}
