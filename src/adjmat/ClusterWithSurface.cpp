/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2019 The plumed team
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

namespace PLMD {
namespace adjmat {

class ClusterWithSurface : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit ClusterWithSurface(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(ClusterWithSurface,"CLUSTER_WITHSURFACE")

void ClusterWithSurface::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords(keys);
  keys.add("optional","RCUT_SURF",""); 
}

ClusterWithSurface::ClusterWithSurface(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  // Read rcut input
  std::string rcut_surf_str; parse("RCUT_SURF",rcut_surf_str);
  // Create a cluster weights object
  readInputLine( getShortcutLabel() + "_wnosurf: CLUSTER_WEIGHTS " + convertInputLineToString() );
  // Now create a contact matrix
  readInputLine( getShortcutLabel() + "_cmat: CONTACT_MATRIX GROUP=" + getShortcutLabel() + "_wnosurf SWITCH={" + rcut_surf_str +"}" );
  // Now create a custom matrix
  readInputLine( getShortcutLabel() + "_cwmat: CUSTOM_MATRIX GROUP1=" + getShortcutLabel() + "_wnosurf FUNC=max");
  // Product of matrices
  readInputLine( getShortcutLabel() + "_pmat: MATHEVAL ARG1=" + getShortcutLabel() + "_cmat.w ARG2=" + getShortcutLabel() + "_cwmat FUNC=x*y PERIODIC=NO");
  // DFS clustering
  readInputLine( getShortcutLabel() + "_clust: DFSCLUSTERING MATRIX=" + getShortcutLabel() + "_pmat"); 
  // And final cluster weights
  readInputLine( getShortcutLabel() + ": CLUSTER_WEIGHTS CLUSTERS=" + getShortcutLabel() + "_clust CLUSTER=1");

}

}
}
