/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
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
#include "multicolvar/MultiColvarBase.h"

namespace PLMD {
namespace adjmat {

class ClusterProperties : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit ClusterProperties(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(ClusterProperties,"CLUSTER_PROPERTIES")

void ClusterProperties::registerKeywords(Keywords& keys){
  ActionShortcut::registerKeywords( keys );
  keys.add("optional","ARG","calculate the sum of the arguments calculated by this action for the cluster");
  multicolvar::MultiColvarBase::shortcutKeywords( keys );
}

ClusterProperties::ClusterProperties(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  // Read the property we are interested in
  std::string argstr; parse("ARG",argstr);
  // Read in the shortcut keywords
  std::map<std::string,std::string> keymap; multicolvar::MultiColvarBase::readShortcutKeywords( keymap, this );
  // Create a cluster weights object
  readInputLine( getShortcutLabel() + ": CLUSTER_WEIGHTS FROM_PROPERTIES=true " + convertInputLineToString() );
  // Now do the multicolvar bit
  multicolvar::MultiColvarBase::expandFunctions( getShortcutLabel(), argstr, getShortcutLabel(), keymap, this );
}

}
}
