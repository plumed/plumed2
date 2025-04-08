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

//+PLUMEDOC CONCOMP CLUSTER_DIAMETER
/*
Print out the diameter of one of the connected components.

An example input that determines the diameter of the second largest cluster that is identified by
analysing the connected components of a [CONTACT_MATRIX](CONTACT_MATRIX.md) using [DFSCLUSTERING](DFSCLUSTERING.md) is shown below:

```plumed
# Calculate a contact matrix between the first 100 atoms in the configuration
cm: CONTACT_MATRIX GROUP=1-100 SWITCH={CUBIC D_0=0.45  D_MAX=0.55}
# Find the connected components from the contact matrix
dfs: DFSCLUSTERING ARG=cm
# Returns a 100-dimensional vector that is 1 if the correpsonding atom index is in the largest cluster and is zero otherwise
clust: CLUSTER_WEIGHTS CLUSTERS=dfs CLUSTER=1
# And determine the size of the largest cluster that was identified
c1: CLUSTER_DIAMETER ARG=clust ATOMS=1-100
PRINT ARG=c1 FILE=colvar
```

The diameter here is defined as the largest of all the distances between the pairs of atom in the cluster

__The output from this action is NOT differentiable__

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace clusters {

class ClusterDiameter : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit ClusterDiameter(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(ClusterDiameter,"CLUSTER_DIAMETER")

void ClusterDiameter::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("optional","ARG","calculate ths radius of the cluster that are in this particular cluster");
  keys.add("compulsory","ATOMS","the atoms that were used to calculate the matrix that was clustered");
  keys.setValueDescription("scalar","the largest of all the distances between the pairs of atom in the cluster");
  keys.needsAction("DISTANCE_MATRIX");
  keys.needsAction("OUTER_PRODUCT");
  keys.needsAction("CUSTOM");
  keys.needsAction("FLATTEN");
  keys.needsAction("HIGHEST");
  keys.addDOI("https://doi.org/10.1021/acs.jctc.6b01073");
}

ClusterDiameter::ClusterDiameter(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // Read in the argument
  std::string arg_str, atdata;
  parse("ARG",arg_str);
  parse("ATOMS",atdata);
  // Distance matrix
  readInputLine( getShortcutLabel() + "_dmat: DISTANCE_MATRIX GROUP=" + atdata );
  // Matrix of bonds in cluster
  readInputLine( getShortcutLabel() + "_bmat: OUTER_PRODUCT FUNC=x*y ARG=" + arg_str + "," + arg_str + " MASK=" + getShortcutLabel() + "_dmat" );
  // Product of matrices
  readInputLine( getShortcutLabel() + "_dcls: CUSTOM ARG=" + getShortcutLabel() + "_dmat," + getShortcutLabel() + "_bmat FUNC=x*y PERIODIC=NO");
  // Convert matrix to a vector to get highest
  readInputLine( getShortcutLabel() + "_vdcls: FLATTEN ARG=" + getShortcutLabel() + "_dcls" );
  // And take the highest value
  readInputLine( getShortcutLabel() + ": HIGHEST ARG=" + getShortcutLabel() + "_vdcls");
}

}
}
