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
#include "ActionWithVessel.h"

namespace PLMD {
namespace vesselbase {

class Mean : public FunctionVessel {
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  explicit Mean( const vesselbase::VesselOptions& da );
  std::string value_descriptor();
  double calcTransform( const double& val, double& dv ) const ;
};

PLUMED_REGISTER_VESSEL(Mean,"MEAN")

void Mean::registerKeywords( Keywords& keys ) {
  FunctionVessel::registerKeywords(keys);
}

void Mean::reserveKeyword( Keywords& keys ) {
  keys.reserve("vessel","MEAN","take the mean of these variables.");
  keys.addOutputComponent("mean","MEAN","the mean value. The output component can be referred to elsewhere in the input "
                          "file by using the label.mean");
}

Mean::Mean( const vesselbase::VesselOptions& da ) :
  FunctionVessel(da)
{
  if( getAction()->isPeriodic() ) error("MEAN cannot be used with periodic variables");
  norm=true;   // Makes sure we calculate the average
}

std::string Mean::value_descriptor() {
  return "the mean value";
}

double Mean::calcTransform( const double& val, double& dv ) const {
  dv=1.0; return val;
}

}
}
