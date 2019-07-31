/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2019 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

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
#include "FunctionVessel.h"
#include "tools/SwitchingFunction.h"
#include "ActionWithVessel.h"

namespace PLMD {
namespace vesselbase {

class MoreThan : public FunctionVessel {
private:
  SwitchingFunction sf;
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  explicit MoreThan( const VesselOptions& da );
  std::string value_descriptor() override;
  double calcTransform( const double& val, double& dv ) const override;
};

PLUMED_REGISTER_VESSEL(MoreThan,"MORE_THAN")

void MoreThan::registerKeywords( Keywords& keys ) {
  FunctionVessel::registerKeywords( keys );
  SwitchingFunction::registerKeywords( keys );
}

void MoreThan::reserveKeyword( Keywords& keys ) {
  keys.reserve("vessel","MORE_THAN","calculate the number of variables more than a certain target value. "
               "This quantity is calculated using \\f$\\sum_i 1.0 - \\sigma(s_i)\\f$, where \\f$\\sigma(s)\\f$ "
               "is a \\ref switchingfunction.");
  keys.addOutputComponent("morethan","MORE_THAN","the number of values more than a target value. This is calculated using one of the "
                          "formula described in the description of the keyword so as to make it continuous. "
                          "You can calculate this quantity multiple times using different parameters.");
}

MoreThan::MoreThan( const VesselOptions& da ) :
  FunctionVessel(da)
{
  usetol=true;
  if( getAction()->isPeriodic() ) error("more than is not a meaningful option for periodic variables");
  std::string errormsg; sf.set( getAllInput(), errormsg );
  if( errormsg.size()!=0 ) error( errormsg );
}

std::string MoreThan::value_descriptor() {
  return "the number of values more than " + sf.description();
}

double MoreThan::calcTransform( const double& val, double& dv ) const {
  double f = 1.0 - sf.calculate(val, dv); dv*=-val; return f;
}

}
}
