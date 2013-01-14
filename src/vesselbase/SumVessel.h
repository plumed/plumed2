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
#ifndef __PLUMED_vesselbase_SumVessel_h
#define __PLUMED_vesselbase_SumVessel_h

#include <string>
#include <cstring>
#include <vector>
#include "VesselAccumulator.h"

namespace PLMD {
namespace vesselbase{

class SumVessel : public VesselAccumulator {
public:
  SumVessel( const VesselOptions& );
/// This retrieves data from action and calculates
  bool calculate( const unsigned& , const double& );
/// Compute the ith component and the derivatives
  virtual double compute( const unsigned& , const double& , double& )=0;
/// This does the final step of the calculation
  void finish( const double& tolerance );
/// Do any final compuations
  virtual double final_computations( const unsigned& , const double& , double& );
};

}
}
#endif
