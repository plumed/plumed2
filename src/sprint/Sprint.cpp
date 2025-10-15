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
#include "core/ActionWithValue.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

//+PLUMEDOC MATRIXF SPRINT
/*
Calculate SPRINT topological variables from an adjacency matrix.

The SPRINT topological variables are calculated from the largest eigenvalue, $\lambda$ of
an $n\times n$ adjacency matrix and its corresponding eigenvector, $V$, using:

$$
s_i = \sqrt{n} \lambda v_i
$$

The example input below calculates the 7 SPRINT coordinates for a 7 atom cluster of Lennard-Jones
atoms and prints their values to a file.

```plumed
ss: SPRINT GROUP=1-7 SWITCH={RATIONAL R_0=0.1}
PRINT ARG=ss.* FILE=colvar
```

This example input calculates the 14 SPRINT coordinates for a molecule composed of 7 hydrogen and
7 carbon atoms.

```plumed
ss: SPRINT ...
  GROUP1=1-7 GROUP2=8-14
  SWITCH11={RATIONAL R_0=2.6 NN=6 MM=12}
  SWITCH12={RATIONAL R_0=2.2 NN=6 MM=12}
  SWITCH22={RATIONAL R_0=2.2 NN=6 MM=12}
...

PRINT ARG=ss.* FILE=colvar
```

If you explore the inputs above you can see that when PLUMED reads them it creates a more complicated input file
for calculating the SPRINT CVs. You can get a sense of how these CVs are calculated by looking at the
expanded versions of the shortcuts in the inputs above. The insight into these methods that you can obtain by looking
at these expanded input should hopefully give you ideas for developing new versions of these methods that use the same
body of theory.  For example, if you look at the inputs above you can see that one or more [CONTACT_MATRIX](CONTACT_MATRIX.md) actions are
used to calculate sprint.  These CONTACT_MATRIX determine whether atoms are adjacent or not.  However, you can
use different quantities to measure whether or not two given atoms/molecules are
adjacent or not and compute a different type of adjacency matrix. For example you can say that two molecules are
adjacent if they are within a certain distance of each other and if they have similar orientations or you can argue that
two molecules are adjacent if there is a hydrogen bond between them.

In the example input below we measure adjacency between atoms using a [BRIDGE_MATRIX](BRIDGE_MATRIX.md) and then compute the SPRINT CVs
by diagonalizing this matrix.

```plumed
b: BRIDGE_MATRIX GROUP=1-7 BRIDGING_ATOMS=8-100 SWITCH={RATIONAL R_0=0.2}
s: SPRINT MATRIX=b
PRINT ARG=s.* FILE=colvar
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace sprint {

class Sprint : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit Sprint(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(Sprint,"SPRINT")

void Sprint::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("optional","MATRIX","the matrix that you would like to perform SPRINT on");
  keys.add("numbered","GROUP","specifies the list of atoms that should be assumed indistinguishable");
  keys.add("numbered","SWITCH","specify the switching function to use between two sets of indistinguishable atoms");
  keys.needsAction("CONTACT_MATRIX");
  keys.needsAction("DIAGONALIZE");
  keys.needsAction("CUSTOM");
  keys.needsAction("SELECT_COMPONENTS");
  keys.needsAction("SORT");
  keys.needsAction("COMBINE");
  keys.addOutputComponent("coord","default","scalar","the sprint coordinates");
  keys.addDOI("10.1103/PhysRevLett.107.085504");
}

Sprint::Sprint(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  std::string matinp;
  parse("MATRIX",matinp);
  if( matinp.length()==0 ) {
    readInputLine( getShortcutLabel() + "_jmat: CONTACT_MATRIX " + convertInputLineToString() );
    matinp = getShortcutLabel() + "_jmat";
  }
  std::vector<unsigned> nin_group;
  unsigned ntot_atoms=0;
  for(unsigned i=1;; ++i) {
    std::string inum;
    Tools::convert( i, inum );
    ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( matinp + inum + inum );
    if( !av ) {
      break ;
    }
    unsigned natoms = (av->copyOutput(0))->getShape()[0];
    nin_group.push_back( natoms );
    ntot_atoms += natoms;
  }
  if( nin_group.size()==0 ) {
    ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( matinp );
    unsigned natoms = (av->copyOutput(0))->getShape()[0];
    nin_group.push_back( natoms );
    ntot_atoms = natoms;
  }

  // Diagonalization
  readInputLine( getShortcutLabel() + "_diag: DIAGONALIZE ARG=" + matinp + " VECTORS=1");
  // Compute sprint coordinates as product of eigenvalue and eigenvector times square root of number of atoms in all groups
  std::string str_natoms;
  Tools::convert( ntot_atoms, str_natoms );
  readInputLine( getShortcutLabel() + "_sp: CUSTOM ARG=" + getShortcutLabel() + "_diag.vecs-1," + getShortcutLabel() +
                 "_diag.vals-1 FUNC=sqrt(" + str_natoms + ")*x*y PERIODIC=NO");
  // Sort sprint coordinates for each group of atoms
  unsigned k=0, kk=0;
  for(unsigned j=0; j<nin_group.size(); ++j) {
    std::string jnum;
    std::string knum;
    std::string nnum;
    Tools::convert( j+1, jnum );
    Tools::convert(k+1, knum);
    k++;
    std::string sort_act = getShortcutLabel() + "_selection" + jnum + ": SELECT_COMPONENTS ARG=" + getShortcutLabel() + "_sp COMPONENTS=" + knum;
    for(unsigned n=1; n<nin_group[j]; ++n) {
      Tools::convert( k+1, knum );
      sort_act += ","+ knum;
      k++;
    }
    readInputLine( sort_act );
    readInputLine( getShortcutLabel() + jnum + ": SORT ARG=" + getShortcutLabel() + "_selection" + jnum );
    for(unsigned n=0; n<nin_group[j]; ++n) {
      Tools::convert( kk, knum );
      Tools::convert( n+1, nnum );
      kk++;
      readInputLine( getShortcutLabel() + "_coord-" + knum + ": COMBINE ARG=" + getShortcutLabel() + jnum + "." + nnum + " PERIODIC=NO" );
    }
  }
}

}
}
