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

class more_than : public SumVessel {
private:
  SwitchingFunction sf;
public:
  static void reserveKeyword( Keywords& keys );
  more_than( const VesselOptions& da );
  double compute( const unsigned& i, const double& val, double& df ); 
  void printKeywords( Log& log );
};

PLUMED_REGISTER_VESSEL(more_than,"MORE_THAN")

void more_than::reserveKeyword( Keywords& keys ){
  keys.reserve("numbered","MORE_THAN","take the number of variables more than the specified target and store it in a value called gt<target>. "
                          "This is calculated using \\f$1.0 - s(r)\\f$, where \\f$s(r)\\f$ is a switching function. " +
                           SwitchingFunction::documentation() );
}

more_than::more_than( const VesselOptions& da ) :
SumVessel(da)
{
  std::string errormsg;
  sf.set( da.parameters, errormsg );
  if( errormsg.size()!=0 ) error( errormsg );
  std::string vv; Tools::convert( sf.get_r0(), vv );
  addOutput("gt" + vv);
  log.printf("  value %s.gt%s contains the number of values more than %s\n",(getAction()->getLabel()).c_str(),vv.c_str(),(sf.description()).c_str());
}

void more_than::printKeywords( Log& log ){
  sf.printKeywords( log );
}

double more_than::compute( const unsigned& i, const double& val, double& df ){
  plumed_assert( i==0 );
  double f; f=1.0 - sf.calculate(val, df); df*=-val;
  return f;
}

}
