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
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"

//+PLUMEDOC CONCOMP CLUSTER_WITHSURFACE
/*
Determine the atoms that are within a certain cutoff of the atoms in a cluster.

This shortcut action can be used to identify the atoms within the largest connected cluster in a system as well
as the atoms around the cluster as shown in the example input below:

```plumed
# Calculate a matrix with elements that are non zero if the corresponding atoms are within 0.38 of each other
cmat: CONTACT_MATRIX GROUP=1-1996 SWITCH={CUBIC D_0=0.34 D_MAX=0.38}
# Calculate the coordination numbers of all the atoms
ones: ONES SIZE=1996
c1: MATRIX_VECTOR_PRODUCT ARG=cmat,ones
# Now identify those atoms that have a coordination number that is less than 13.5
cf: LESS_THAN ARG=c1 SWITCH={CUBIC D_0=13 D_MAX=13.5}

# Now construct an adjacency matrix with elements that are equal to one if atoms are within 0.38 of each other
# and if they both have a coordination number that is less than 13.5
cmat2: CONTACT_MATRIX GROUP=1-1996 SWITCH={CUBIC D_0=0.34 D_MAX=0.38}
dmat: OUTER_PRODUCT ARG=cf,cf
mat: CUSTOM ARG=cmat2,dmat FUNC=x*y PERIODIC=NO
# Find the connected componets of this matrix
dfs: DFSCLUSTERING ARG=mat

# Now identify the atoms that are within the largest cluster and the atoms that are within 0.38 nm of the
# atoms that are within this cluster.  Determine how many atoms in total satisfy this condition and the
# distance betwen the two atoms in this set that are farthest appart.
clust2a: CLUSTER_WITHSURFACE ATOMS=1-1996 CLUSTERS=dfs CLUSTER=1 RCUT_SURF={CUBIC D_0=0.34 D_MAX=0.38}
size2a: SUM ARG=clust2a PERIODIC=NO
dia2a: CLUSTER_DIAMETER ATOMS=1-1996 ARG=clust2a
PRINT ARG=size2a,dia2a
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace clusters {

class ClusterWithSurface : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit ClusterWithSurface(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(ClusterWithSurface,"CLUSTER_WITHSURFACE")

void ClusterWithSurface::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords(keys);
  keys.add("optional","RCUT_SURF","the cutoff to use for determining which atoms are connected to the surface of the cluster");
  keys.add("compulsory","ATOMS","the atoms that were used to calculate the matrix that was clustered");
  keys.add("compulsory","CLUSTERS","the label of the action that does the clustering");
  keys.add("compulsory","CLUSTER","1","which cluster would you like to look at 1 is the largest cluster, 2 is the second largest, 3 is the the third largest and so on.");
  keys.setValueDescription("vector","a vector that is one for those atoms that are within the cluster or that are within a cetain cutoff of one of the atoms in the cluster and zero otherwise");
  keys.needsAction("CLUSTER_WEIGHTS");
  keys.needsAction("CONTACT_MATRIX");
  keys.needsAction("OUTER_PRODUCT");
  keys.needsAction("CUSTOM");
  keys.needsAction("DFSCLUSTERING");
  keys.addDOI("https://doi.org/10.1021/acs.jctc.6b01073");
}

ClusterWithSurface::ClusterWithSurface(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // Read atoms for contact matrix
  std::string atdata;
  parse("ATOMS",atdata);
  // Read rcut input
  std::string rcut_surf_str;
  parse("RCUT_SURF",rcut_surf_str);
  // Create a cluster weights object
  readInputLine( getShortcutLabel() + "_wnosurf: CLUSTER_WEIGHTS " + convertInputLineToString() );
  // Now create a contact matrix
  readInputLine( getShortcutLabel() + "_cmat: CONTACT_MATRIX GROUP=" + atdata + " SWITCH={" + rcut_surf_str +"}" );
  // Now create a custom matrix
  readInputLine( getShortcutLabel() + "_cwmat: OUTER_PRODUCT ARG=" + getShortcutLabel() + "_wnosurf," + getShortcutLabel() + "_wnosurf FUNC=max MASK=" + getShortcutLabel() + "_cmat");
  // Product of matrices
  readInputLine( getShortcutLabel() + "_pmat: CUSTOM ARG=" + getShortcutLabel() + "_cmat," + getShortcutLabel() + "_cwmat FUNC=x*y PERIODIC=NO");
  // DFS clustering
  readInputLine( getShortcutLabel() + "_clust: DFSCLUSTERING ARG=" + getShortcutLabel() + "_pmat");
  // And final cluster weights
  readInputLine( getShortcutLabel() + ": CLUSTER_WEIGHTS CLUSTERS=" + getShortcutLabel() + "_clust CLUSTER=1");
}

}
}
