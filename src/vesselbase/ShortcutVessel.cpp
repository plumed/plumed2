/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#include "ShortcutVessel.h"
#include "ActionWithVessel.h"

namespace PLMD {
namespace vesselbase {

void ShortcutVessel::registerKeywords( Keywords& keys ) {
  Vessel::registerKeywords( keys ); keys.remove("LABEL");
  plumed_assert( keys.size()==0 );
}

ShortcutVessel::ShortcutVessel( const VesselOptions& da):
  Vessel(da)
{
}

void ShortcutVessel::addVessel( const std::string& name, const std::string& input ) {
  unsigned numlab=1;
  for(unsigned i=0; i<(getAction()->functions).size(); ++i) {
    if( (getAction()->functions[i])->getName()==name ) numlab++;
  }
  getAction()->addVessel( name, input, numlab );
}

}
}
