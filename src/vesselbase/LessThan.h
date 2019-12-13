/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2019 The plumed team
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
#ifndef __PLUMED_vesselbase_LessThan_h
#define __PLUMED_vesselbase_LessThan_h

#include "FunctionVessel.h"
#include "tools/SwitchingFunction.h"
#include "ActionWithVessel.h"

namespace PLMD {
namespace vesselbase {

class LessThan : public FunctionVessel {
private:
  SwitchingFunction sf;
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  explicit LessThan( const VesselOptions& da );
  std::string value_descriptor() override;
  double calcTransform( const double& val, double& dv ) const override;
  double getCutoff();
};

}
}
#endif
