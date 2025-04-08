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

//+PLUMEDOC CONCOMP CLUSTER_NATOMS
/*
Calculate the number of atoms in the cluster of interest

An example input that determines the number of atoms in the second largest cluster that is identified by
analysing the connected components of a [CONTACT_MATRIX](CONTACT_MATRIX.md) using [DFSCLUSTERING](DFSCLUSTERING.md) is shown below:

```plumed
# Calculate a contact matrix between the first 100 atoms in the configuration
cm: CONTACT_MATRIX GROUP=1-100 SWITCH={CUBIC D_0=0.45  D_MAX=0.55}
# Find the connected components from the contact matrix
dfs: DFSCLUSTERING ARG=cm
# And determine the size of the second largest cluster that was identified
c1: CLUSTER_NATOMS CLUSTERS=dfs CLUSTER=2
PRINT ARG=c1 FILE=colvar
```

__The output from this action is NOT differentiable__

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace clusters {

class ClusterNatoms : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit ClusterNatoms(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(ClusterNatoms,"CLUSTER_NATOMS")

void ClusterNatoms::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","CLUSTERS","the label of the action that does the clustering");
  keys.add("compulsory","CLUSTER","1","which cluster would you like to look at 1 is the largest cluster, 2 is the second largest, 3 is the the third largest and so on.");
  keys.setValueDescription("scalar","the number of atoms in the cluster");
  keys.needsAction("CLUSTER_WEIGHTS");
  keys.needsAction("SUM");
}

ClusterNatoms::ClusterNatoms(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // Create a cluster weights object
  readInputLine( getShortcutLabel() + "_weights: CLUSTER_WEIGHTS " + convertInputLineToString() );
  // Add all the weights together (weights are 1 or 0)
  readInputLine( getShortcutLabel() + ": SUM ARG=" + getShortcutLabel() + "_weights PERIODIC=NO");
}

}
}
