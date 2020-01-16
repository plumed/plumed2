/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "HistogramBase.h"
#include "SphericalKDE.h"

namespace PLMD {
namespace gridtools {

class SphericalKDEShortcut : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit SphericalKDEShortcut(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(SphericalKDEShortcut,"SPHERICAL_KDE")

void SphericalKDEShortcut::registerKeywords( Keywords& keys ) {
  SphericalKDE::registerKeywords( keys ); 
}

SphericalKDEShortcut::SphericalKDEShortcut(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  std::string height, height_str; parse("HEIGHTS",height); if( height.length()>0 ) height_str = " HEIGHTS=" + height;
  HistogramBase::createKDEObject( getShortcutLabel(), "SPHERICAL_KDE", height, height_str, this );
}

}
}
