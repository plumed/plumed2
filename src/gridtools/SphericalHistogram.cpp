/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2020 The plumed team
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
#include "SphericalKDE.h"
#include "core/ActionRegister.h"
#include "core/ActionShortcut.h"

namespace PLMD {
namespace gridtools {

class SphericalHistogram : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit SphericalHistogram(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(SphericalHistogram,"SPHERICAL_HISTOGRAM")

void SphericalHistogram::registerKeywords( Keywords& keys ) {
  SphericalKDE::registerKeywords( keys );
  HistogramBase::histogramKeywords( keys );
  keys.use("UPDATE_FROM"); keys.use("UPDATE_UNTIL");
}

SphericalHistogram::SphericalHistogram(const ActionOptions&ao):
Action(ao),
ActionShortcut(ao)
{
    // Parse the keymap
    std::map<std::string,std::string> keymap; HistogramBase::readHistogramKeywords( keymap, this ); 
    // Make the kde object
    readInputLine( getShortcutLabel() + "_kde: SPHERICAL_KDE " + convertInputLineToString(), true ); 
    // And the averaging object
    HistogramBase::createAveragingObject( getShortcutLabel() + "_kde", getShortcutLabel(), keymap, this );
}

}
}
