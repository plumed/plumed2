/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#include "OrderingVessel.h"
#include "VesselRegister.h"

namespace PLMD {
namespace vesselbase {

class Highest : public OrderingVessel {
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  explicit Highest( const VesselOptions& da );
  std::string value_descriptor() override;
  bool compare( const double&, const double& ) override;
};

PLUMED_REGISTER_VESSEL(Highest,"HIGHEST")

void Highest::registerKeywords( Keywords& keys ) {
  OrderingVessel::registerKeywords( keys );
}

void Highest::reserveKeyword( Keywords& keys ) {
  keys.reserve("vessel","HIGHEST","this flag allows you to recover the highest of these variables.");
  keys.addOutputComponent("highest","HIGHEST","the highest of the quantities calculated by this action");
}

Highest::Highest( const VesselOptions& da ) :
  OrderingVessel(da)
{
}

std::string Highest::value_descriptor() {
  return "the highest of the individual colvar values";
}

bool Highest::compare( const double& val1, const double& val2 ) {
  return val1>val2;
}

}
}
