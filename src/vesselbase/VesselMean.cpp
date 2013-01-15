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
#include "WeightedSumVessel.h"
#include "ActionWithVessel.h"

namespace PLMD {
namespace vesselbase {

class VesselMean : public vesselbase::WeightedSumVessel {
private:
public:
  static void reserveKeyword( Keywords& keys );
  VesselMean( const vesselbase::VesselOptions& da );
  double getWeight( const unsigned& i, bool& hasDerivatives );
  double compute( const unsigned& i, const unsigned& j, const double& val, double& df );
};

PLUMED_REGISTER_VESSEL(VesselMean,"AVERAGE")

void VesselMean::reserveKeyword( Keywords& keys ){
  keys.reserveFlag("AVERAGE",false,"take the average value of these variables and store it in value called average.");
}

VesselMean::VesselMean( const vesselbase::VesselOptions& da ) :
WeightedSumVessel(da)
{
  if( getAction()->isPeriodic() ) error("MEAN cannot be used with periodic variables");

  useNorm();
  addOutput("average");
  log.printf("  value %s.average contains the average value\n",(getAction()->getLabel()).c_str());
}

double VesselMean::compute( const unsigned& i, const unsigned& j, const double& val, double& df ){
  plumed_dbg_assert( j==0 );
  df=1.0; return val;
}

double VesselMean::getWeight( const unsigned& i, bool& hasDerivatives ){
  plumed_dbg_assert( !getAction()->isPossibleToSkip() ); 
  hasDerivatives=false;
  return 1.0;
}

}
}
