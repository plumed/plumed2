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
#include "FunctionVessel.h"
#include "VesselRegister.h"

namespace PLMD {
namespace vesselbase {

class Sum : public FunctionVessel {
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  explicit Sum( const VesselOptions& da );
  std::string value_descriptor();
  double calcTransform( const double& val, double& dv ) const ;
};

PLUMED_REGISTER_VESSEL(Sum,"SUM")

void Sum::registerKeywords( Keywords& keys ) {
  FunctionVessel::registerKeywords( keys );
}

void Sum::reserveKeyword( Keywords& keys ) {
  keys.reserve("vessel","SUM","calculate the sum of all the quantities.");
  keys.addOutputComponent("sum","SUM","the sum of values");
}

Sum::Sum( const VesselOptions& da ) :
  FunctionVessel(da)
{
}

std::string Sum::value_descriptor() {
  return "the sum of all the values";
}

double Sum::calcTransform( const double& val, double& dv ) const {
  dv=1.0; return val;
}

}
}
