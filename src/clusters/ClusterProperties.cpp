/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2020 The plumed team
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
#include "multicolvar/MultiColvarShortcuts.h"

//+PLUMEDOC CONCOMP CLUSTER_PROPERTIES
/*
Calculate properties of the distribution of some quantities that are part of a connected component

__This shortcut action is present to ensure that inputs for older PLUMED versions remain compatible.  We STRONGLY
encourage you to use the newer (and simpler) syntax for clustering calculations__  Look at the documentation for
CLUSTER_WEIGHTS and the expanded version of the input below to see how this new input syntax operates.

This collective variable was developed for looking at nucleation phenomena, where you are
interested in using studying the behavior of atoms in small aggregates or nuclei.  In these sorts of
problems you might be interested in the degree the atoms in a nucleus have adopted their crystalline
structure or (in the case of heterogeneous nucleation of a solute from a solvent) you might be
interested in how many atoms are present in the largest cluster.

The input below calculates the coordination numbers of atoms 1-100 and then computes the an adjacency
matrix whose elements measures whether atoms $i$ and $j$ are within 0.55 nm of each other.  The action
labelled dfs then treats the elements of this matrix as zero or ones and thus thinks of the matrix as defining
a graph.  This dfs action then finds the largest connected component in this graph.  The sum of the coordination
numbers for the atoms in this largest connected component are then computed and this quantity is output to a colvar
file.

```plumed
lq: COORDINATIONNUMBER SPECIES=1-100 SWITCH={CUBIC D_0=0.45  D_MAX=0.55}
cm: CONTACT_MATRIX GROUP=lq  SWITCH={CUBIC D_0=0.45  D_MAX=0.55}
dfs: DFSCLUSTERING ARG=cm
clust1: CLUSTER_PROPERTIES ARG=lq CLUSTERS=dfs CLUSTER=1 SUM
PRINT ARG=clust1.sum FILE=colvar
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace clusters {

class ClusterProperties : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit ClusterProperties(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(ClusterProperties,"CLUSTER_PROPERTIES")

void ClusterProperties::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.setDeprecated("CLUSTER_WEIGHTS");
  keys.add("compulsory","ARG","calculate the sum of the arguments calculated by this action for the cluster");
  keys.add("compulsory","CLUSTERS","the label of the action that does the clustering");
  keys.add("compulsory","CLUSTER","1","which cluster would you like to look at 1 is the largest cluster, 2 is the second largest, 3 is the the third largest and so on.");
  keys.setValueDescription("vector","a vector that is one if the atom is part of the cluster or interest and zero otherwise");
  multicolvar::MultiColvarShortcuts::shortcutKeywords( keys );
  keys.needsAction("CLUSTER_WEIGHTS");
}

ClusterProperties::ClusterProperties(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // Read the property we are interested in
  std::string argstr;
  parse("ARG",argstr);
  // Read in the shortcut keywords
  std::map<std::string,std::string> keymap;
  multicolvar::MultiColvarShortcuts::readShortcutKeywords( keymap, this );
  // Create a cluster weights object
  readInputLine( getShortcutLabel() + ": CLUSTER_WEIGHTS " + convertInputLineToString() );
  // Now do the multicolvar bit
  multicolvar::MultiColvarShortcuts::expandFunctions( getShortcutLabel(), argstr, getShortcutLabel(), keymap, this );
}

}
}
