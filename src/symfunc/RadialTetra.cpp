/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The plumed team
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
#include "core/ActionRegister.h"
#include "core/ActionShortcut.h"
#include "core/ActionWithValue.h"
#include "CoordinationNumbers.h"
#include "multicolvar/MultiColvarShortcuts.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

//+PLUMEDOC MCOLVAR TETRA_RADIAL
/*
Calculate the radial tetra CV

This shortcut calculates a [symmetry function](https://www.plumed-tutorials.org/lessons/23/001/data/SymmetryFunction.html). The particular function that is being
evaluated for the coordination sphere here is as follows:

$$
s_i = 1 - \frac{\sum_{j=1}^4 r_{ij}^2 - z_i\sum_{j=1}^4 r_{ij}}{12 z_i^2} \qquad \textrm{where} \qquad z_i = \frac{1}{4} \sum_{j=1}^4 r_{ij}
$$

In this expression the 4 atoms in the sums over $j$ are the four atoms that are nearest to atom $i$ and $r_{ij}$ is the distance between atoms $i$ and $j$.
The CV is large if the four atoms nearest atom $i$ are arranged on the vertices of a regular tetrahedron
and small otherwise.  The following example shows how you can use this action to measure the degree of tetrahedral order in a system.

```plumed
# Calculate a vector that contains 64 values for the symmetry function.
acv: TETRA_RADIAL SPECIES=1-64
# Sum the elements of the vector and calculate the mean value on the atoms from this sum.
acv_sum: SUM ARG=acv PERIODIC=NO
acv_mean: CUSTOM ARG=acv FUNC=x/64 PERIODIC=NO
# Print out the positions of the 64 atoms for which the symmetry function was calculated
# to an xyz file along with the values of the symmetry function
DUMPATOMS ATOMS=1-64 ARG=acv FILE=mcolv.xyz
# Print out the average value of the symmetry function and the sum of all the symmetry functions
PRINT ARG=acv_sum,acv_mean FILE=colvar
```

In the input above we have only one type of atom.  If you want to measure the degree of tetrahedral order amongst the bonds from atoms of SPECIESA to atoms of SPECIESB
you use an input like the one shown below:

```plumed
acv: TETRA_RADIAL SPECIESA=1-64 SPECIESB=65-128
acv_sum: SUM ARG=acv PERIODIC=NO
acv_mean: CUSTOM ARG=acv FUNC=x/64 PERIODIC=NO
DUMPATOMS ATOMS=1-64 ARG=acv FILE=mcolv.xyz
PRINT ARG=acv_sum,acv_mean FILE=colvar
```

By default the vectors connecting atoms that appear in the expression above are calculated in a way that takes periodic boundary conditions into account.  If you want
not to take the periodic boundary conditions into account for any reason you use the NOPBC flag as shown below:

```plumed
acv: TETRA_RADIAL SPECIES=1-64 NOPBC
acv_sum: SUM ARG=acv PERIODIC=NO
acv_mean: CUSTOM ARG=acv FUNC=x/64 PERIODIC=NO
DUMPATOMS ATOMS=1-64 ARG=acv FILE=mcolv.xyz
PRINT ARG=acv_sum,acv_mean FILE=colvar
```

## Deprecated syntax

More information on the deprecated keywords that are given below is available in the documentation for the [DISTANCES](DISTANCES.md) command.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace symfunc {

class RadialTetra : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit RadialTetra(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(RadialTetra,"TETRA_RADIAL")

void RadialTetra::registerKeywords( Keywords& keys ) {
  CoordinationNumbers::shortcutKeywords( keys );
  keys.remove("MASK");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.add("compulsory","CUTOFF","-1","ignore distances that have a value larger than this cutoff");
  keys.setValueDescription("vector","the value of the radial tetrahedrality parameter for each of the input atoms");
  keys.remove("NN");
  keys.remove("MM");
  keys.remove("D_0");
  keys.remove("R_0");
  keys.remove("SWITCH");
  keys.needsAction("DISTANCE_MATRIX");
  keys.needsAction("NEIGHBORS");
  keys.needsAction("CUSTOM");
  keys.needsAction("ONES");
  keys.needsAction("MATRIX_VECTOR_PRODUCT");
}

RadialTetra::RadialTetra( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // Read species input and create the matrix
  bool nopbc;
  parseFlag("NOPBC",nopbc);
  std::string pbcstr="";
  if( nopbc ) {
    pbcstr = " NOPBC";
  }
  std::string sp_str, rcut;
  parse("SPECIES",sp_str);
  parse("CUTOFF",rcut);
  if( sp_str.length()>0 ) {
    readInputLine( getShortcutLabel() + "_mat: DISTANCE_MATRIX GROUP=" + sp_str + " CUTOFF=" + rcut + pbcstr );
  } else {
    std::string specA, specB;
    parse("SPECIESA",specA);
    parse("SPECIESB",specB);
    if( specA.length()==0 ) {
      error("missing input atoms");
    }
    if( specB.length()==0 ) {
      error("missing SPECIESB keyword");
    }
    readInputLine( getShortcutLabel() + "_mat: DISTANCE_MATRIX GROUPA=" + specA + " GROUPB=" + specB + " CUTOFF=" + rcut + pbcstr);
  }
  // Get the neighbors matrix
  readInputLine( getShortcutLabel() + "_neigh: NEIGHBORS ARG=" + getShortcutLabel() + "_mat NLOWEST=4");
  // Now get distance matrix that just contains four nearest distances
  readInputLine( getShortcutLabel() + "_near4: CUSTOM ARG=" + getShortcutLabel() + "_mat," + getShortcutLabel() + "_neigh MASK=" + getShortcutLabel() + "_neigh FUNC=x*y PERIODIC=NO");
  //Now compute sum of four nearest distances
  ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_mat");
  plumed_assert( av && av->getNumberOfComponents()>0 && (av->copyOutput(0))->getRank()==2 );
  std::string size;
  Tools::convert( (av->copyOutput(0))->getShape()[1], size );
  readInputLine( getShortcutLabel() + "_ones: ONES SIZE=" + size );
  readInputLine( getShortcutLabel() + "_sum4: MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_near4," + getShortcutLabel() + "_ones");
  // Now compute squares of four nearest distance
  readInputLine( getShortcutLabel() + "_near4_2: CUSTOM ARG=" + getShortcutLabel() + "_near4 FUNC=x*x PERIODIC=NO");
  // Now compute sum of the squares of the four nearest distances
  readInputLine( getShortcutLabel() + "_sum4_2: MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_near4_2," + getShortcutLabel() + "_ones");
  // Evaluate the average distance to the four nearest neighbors
  readInputLine( getShortcutLabel() + "_meanr: CUSTOM ARG=" + getShortcutLabel() + "_sum4 FUNC=0.25*x PERIODIC=NO");
  // Now evaluate the actual per atom CV
  readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_sum4," + getShortcutLabel() + "_sum4_2," + getShortcutLabel() + "_meanr " +
                 "FUNC=(1-(y-x*z)/(12*z*z)) PERIODIC=NO");
  // And get the things to do with the quantities we have computed
  std::map<std::string,std::string> keymap;
  multicolvar::MultiColvarShortcuts::readShortcutKeywords( keymap, this );
  multicolvar::MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel(), "", keymap, this );
}

}
}
