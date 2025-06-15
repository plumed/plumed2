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

//+PLUMEDOC MCOLVAR BRIDGE
/*
Calculate the number of briding atoms between two groups

This quantity calculates:

$$
f(x) = \sum_{ijk} s_A(r_{ij})s_B(r_{ik})
$$

where the sum over $i$ is over all the "bridging atoms" and
$s_A$ and $s_B$ are switching functions.

The following example instructs plumed to calculate the number of water molecules
that are bridging between atoms 1-10 and atoms 11-20 and to print the value
to a file

```plumed
w1: BRIDGE BRIDGING_ATOMS=100-200 GROUPA=1-10 GROUPB=11-20 SWITCH={RATIONAL R_0=0.2}
PRINT ARG=w1 FILE=colvar
```

You can use the following input instead to calculate the number of atoms that bridge between the
45 pairs of distinct pairs of atoms that can be selected from the set of atoms with indexes between
1 and 10.

```plumed
w1: BRIDGE BRIDGING_ATOMS=100-200 GROUP=1-10 SWITCH={RATIONAL R_0=0.2}
PRINT ARG=w1 FILE=colvar
```

In the above inputs the switching functions $s_A$ and $s_B$ are identical.  If you want to use two different
switching functions you would use an input like the one shown below:

```plumed
w1: BRIDGE ...
  BRIDGING_ATOMS=100-200 GROUPA=1-10 GROUPB=11-20
  SWITCHA={RATIONAL R_0=0.2}
  SWITCHB={RATIONAL R_0=0.4}
...
PRINT ARG=w1 FILE=colvar
```

With this input an atom is calculated as bridging if it is (approximately) within 0.2 nm of one of the atoms in
GROUPA and within 0.4 nm of one of the atoms in GROUPB.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class Bridge : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit Bridge(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(Bridge,"BRIDGE")

void Bridge::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords(keys);
  keys.add("atoms","GROUP","the atoms for which you would like to calculate the adjacency matrix");
  keys.add("atoms","GROUPA","when you are calculating the adjacency matrix between two sets of atoms this keyword is used to specify the atoms along with the keyword GROUPB");
  keys.add("atoms","GROUPB","when you are calculating the adjacency matrix between two sets of atoms this keyword is used to specify the atoms along with the keyword GROUPA");
  keys.add("atoms","BRIDGING_ATOMS","The list of atoms that can form the bridge between the two interesting parts "
           "of the structure.");
  keys.add("optional","SWITCH","The parameters of the two switching functions in the above formula");
  keys.linkActionInDocs("SWITCH","LESS_THAN");
  keys.add("optional","SWITCHA","The switching function on the distance between bridging atoms and the atoms in group A");
  keys.linkActionInDocs("SWITCHA","LESS_THAN");
  keys.add("optional","SWITCHB","The switching function on the distance between the bridging atoms and the atoms in group B");
  keys.linkActionInDocs("SWITCHB","LESS_THAN");
  keys.needsAction("BRIDGE_MATRIX");
  keys.needsAction("SUM");
  keys.setValueDescription("scalar","the number of bridging atoms between the two groups");
}

Bridge::Bridge(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // Need to read in switch
  std::string s_inp, sfinput;
  parse("SWITCH",sfinput);
  if( sfinput.length()>0 ) {
    s_inp += "SWITCH={" + sfinput +"} ";
  }
  std::string sfinputa;
  parse("SWITCHA",sfinputa);
  if( sfinputa.length()>0 ) {
    s_inp += "SWITCHA={" + sfinputa +"} ";
  }
  std::string sfinputb;
  parse("SWITCHB",sfinputb);
  if( sfinputb.length()>0 ) {
    s_inp += "SWITCHB={" + sfinputb +"} ";
  }
  // Create the matrix object
  readInputLine( getShortcutLabel() + "_mat: BRIDGE_MATRIX " + s_inp + convertInputLineToString() );
  // Add all the elements of the matrix together
  readInputLine( getShortcutLabel() + ": SUM ARG=" + getShortcutLabel() + "_mat PERIODIC=NO");
}

}
}
