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
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace adjmat {

//+PLUMEDOC MCOLVAR BRIDGE_MATRIX
/*
Calculate the number of atoms that bridge two parts of a structure

This adjacency matrix is used to implement the [BRIDGE](BRIDGE.md) shortcut. The action outputs a adjacency matrix
in the same way as [CONTACT_MATRIX](CONTACT_MATRIX.md).  However, the  $j,k$ element of the adjacency matrix is calculated
using:

$$
M_{jk} = \sum_i s_A(r_{ij})s_B(r_{ik})
$$

In this expression, the sum runs over all the atoms that were specified using the `BRIDGING_ATOMS` keyword, $s_A$ and
$s_B$ are switching functions, and $r_{ij}$ and $r_{ik}$ are the distances between atom $i$ and $j$ and between atoms
$i$ and $k$.  Less formally, this formula ensures that $j,k$ element of the output matrix is one if there is a bridging
atom between atom $j$ and $k$.

Notice that the BRIDGE_MATRIX action is also a shortcut. If we have a single group of non-bridging atoms then we can compute
the contact matrix $\mathbf{C}$ between the non-bridging and bridging atoms. The BRIDGE_MATRIX, $\mathbf{M}$, can by calculated
from $\mathbf{C}$ by multiplying $\mathbf{C}$ by its transpose.  Similarly, if we have two groups of non-bridging atoms
we can calculate the BRIDGE_MATRIX by calculating a [CONTACT_MATRIX](CONTACT_MATRIX.md) between the first set of non-bridging atoms and
a second [CONTACT_MATRIX](CONTACT_MATRIX.md) between the bridging atoms and the bridging atoms and the second set of non-bridging atoms.
The bridge matrix is the product of these two matrices.

In the following example input atoms 100-200 can serve as bridging atoms between the atoms in GROUPA and GROUPB and the
two switching functions $s_A$ and $s_B$ in the formula above are identical.

```plumed
w1: BRIDGE_MATRIX ...
   BRIDGING_ATOMS=100-200
   GROUPA=1-10 GROUPB=11-20
   SWITCH={RATIONAL R_0=0.2}
...
```

If you use a single GROUP keyword as in the input below as a single SWITCH keyword the output matrix is symmetric.

```plumed
w2: BRIDGE_MATRIX ...
   BRIDGING_ATOMS=100-200 GROUP=1-10
   SWITCH={RATIONAL R_0=0.2}
...
```

However, if the two switching functions are not identical, as in the following example, then the output matrix is __not__ symmetric
even if GROUP is used rather than GROUPA/GROUPB.

```plumed
w2: BRIDGE_MATRIX ...
   BRIDGING_ATOMS=100-200 GROUP=1-10
   SWITCHA={RATIONAL R_0=0.2}
   SWITCHB={RATIONAL R_0=0.4}
...
```

Notice that in all the inputs above the $r_{ij}$ and $r_{ik}$ values that enter the formula above are calculated in a way that takes the
periodic boundary conditions into account.  If you want to ignore the periodic boundary conditions you can use the NOPBC flag as shown below.

```plumed
w2: BRIDGE_MATRIX ...
   BRIDGING_ATOMS=100-200 GROUP=1-10
   SWITCH={RATIONAL R_0=0.2}
   NOPBC
...
```

*/
//+ENDPLUMEDOC

class BridgeMatrix : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  BridgeMatrix(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(BridgeMatrix,"BRIDGE_MATRIX")

void BridgeMatrix::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("atoms","GROUP","the atoms for which you would like to calculate the adjacency matrix");
  keys.add("atoms","GROUPA","when you are calculating the adjacency matrix between two sets of atoms this keyword is used to specify the atoms along with the keyword GROUPB");
  keys.add("atoms","GROUPB","when you are calculating the adjacency matrix between two sets of atoms this keyword is used to specify the atoms along with the keyword GROUPA");
  keys.add("atoms","BRIDGING_ATOMS","The list of atoms that can form the bridge between the two interesting parts "
           "of the structure.");
  keys.add("optional","SWITCH","The parameters of the two switching functions in the above formula");
  keys.linkActionInDocs("SWITCH","LESS_THAN");
  keys.add("optional","SWITCHA","The switching function on the distance between bridging atoms and the atoms in "
           "group A");
  keys.linkActionInDocs("SWITCHA","LESS_THAN");
  keys.setValueDescription("matrix","a matrix containing the weights for the bonds between each pair of atoms");
  keys.add("optional","SWITCHB","The switching function on the distance between the bridging atoms and the atoms in "
           "group B");
  keys.addFlag("NOPBC", false, "don't use periodic boundary conditions");
  keys.linkActionInDocs("SWITCHB","LESS_THAN");
  keys.needsAction("CONTACT_MATRIX");
  keys.needsAction("TRANSPOSE");
  keys.needsAction("MATRIX_PRODUCT");
}

BridgeMatrix::BridgeMatrix(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  std::string b_atoms;
  parse("BRIDGING_ATOMS",b_atoms);
  std::string swstr, swA_str, swB_str;
  parse("SWITCH",swstr);
  if( swstr.length()==0 ) {
    parse("SWITCHA",swA_str);
    if( swA_str.length()==0 ) {
      error("must set either SWITCH or SWITCHA+SWITCHB");
    }
    parse("SWITCHB",swB_str);
    if( swB_str.length()==0 ) {
      error("found SWITCHA but no SWITCHB");
    }
  }
  bool nopbc;
  std::string pbcstr="";
  parseFlag("NOPBC",nopbc);
  if( nopbc ) {
    pbcstr = " NOPBC";
  }

  std::string grp_str;
  if( grp_str.length()>0 ) {
    if( swstr.length()>0 ) {
      readInputLine( getShortcutLabel() + "_cmat: CONTACT_MATRIX GROUPA=" + grp_str + " GROUPB=" + b_atoms + " SWITCH={" + swstr + "}" + pbcstr );
      readInputLine( getShortcutLabel() + "_cmatT: TRANSPOSE ARG=" + getShortcutLabel() + "_cmat" );
    } else {
      readInputLine( getShortcutLabel() + "_cmat: CONTACT_MATRIX GROUPA=" + grp_str + " GROUPB=" + b_atoms + " SWITCH={" + swA_str + "}" + pbcstr );
      readInputLine( getShortcutLabel() + "_cmatT: CONTACT_MATRIX GROUPA=" + b_atoms + " GROUPB=" + grp_str + " SWITCH={" + swB_str + "}" + pbcstr );
    }
    readInputLine( getShortcutLabel() + ": MATRIX_PRODUCT ARG=" + getShortcutLabel() + "_cmat," + getShortcutLabel() + "_cmatT");
  } else {
    if( swA_str.length()==0 ) {
      swA_str = swB_str = swstr;
    }
    std::string grpA_str, grpB_str;
    parse("GROUPA",grpA_str);
    parse("GROUPB",grpB_str);
    readInputLine( getShortcutLabel() + "_cmat: CONTACT_MATRIX GROUPA=" + grpA_str + " GROUPB=" + b_atoms + " SWITCH={" + swA_str + "}" + pbcstr );
    readInputLine( getShortcutLabel() + "_cmatT: CONTACT_MATRIX GROUPA=" + b_atoms + " GROUPB=" + grpB_str + " SWITCH={" + swB_str + "}" + pbcstr );
    readInputLine( getShortcutLabel() + ": MATRIX_PRODUCT ARG=" + getShortcutLabel() + "_cmat," + getShortcutLabel() + "_cmatT");
  }
}

}
}
