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
#include "CoordinationNumbers.h"
#include "multicolvar/MultiColvarShortcuts.h"

//+PLUMEDOC MCOLVAR TETRA_ANGULAR
/*
Calculate the angular tetra CV

This shortcut calculates a [symmetry function](https://www.plumed-tutorials.org/lessons/23/001/data/SymmetryFunction.html). The particular function that is being
evaluated for the coordination sphere here is as follows:

$$
s_i = 1 - \frac{3}{8}\sum_{j=2}^4 \sum_{k=1}^{j-1} \left( \cos(\theta_{jik} + \frac{1}{2} \right)^2
$$

In this expression the 4 atoms in the sums over $j$ and $k$ are the four atoms that are nearest to atom $i$.  $\theta_{jik}$ is the angle between the vectors connecting
atoms $i$ and $j$ and the vector connecting atoms $i$ and $k$.  The CV is large if the four atoms nearest atom $i$ are arranged on the vertices of a regular tetrahedron
and small otherwise.  The following example shows how you can use this action to measure the degree of tetrahedral order in a system.

```plumed
# Calculate a vector that contains 64 values for the symmetry function.
acv: TETRA_ANGULAR SPECIES=1-64
# Sum the elements of the vector and calculate the mean value on the atoms from this sum.
acv_sum: SUM ARG=acv PERIODIC=NO
acv_mean: MEAN ARG=acv PERIODIC=NO
# Print out the positions of the 64 atoms for which the symmetry function was calculated
# to an xyz file along with the values of the symmetry function
DUMPATOMS ATOMS=1-64 ARG=acv FILE=mcolv.xyz
# Print out the average value of the symmetry function and the sum of all the symmetry functions
PRINT ARG=acv_sum,acv_mean FILE=colvar
```

In the input above we have only one type of atom.  If you want to measure the degree of tetrahedral order amongst the bonds from atoms of SPECIESA to atoms of SPECIESB
you use an input like the one shown below:

```plumed
acv: TETRA_ANGULAR SPECIESA=1-64 SPECIESB=65-128
acv_sum: SUM ARG=acv PERIODIC=NO
acv_mean: MEAN ARG=acv PERIODIC=NO
DUMPATOMS ATOMS=1-64 ARG=acv FILE=mcolv.xyz
PRINT ARG=acv_sum,acv_mean FILE=colvar
```

By default the vectors connecting atoms that appear in the expression above are calculated in a way that takes periodic boundary conditions into account.  If you want
not to take the periodic boundary conditions into account for any reason you use the NOPBC flag as shown below:

```plumed
acv: TETRA_ANGULAR SPECIES=1-64 NOPBC
acv_sum: SUM ARG=acv PERIODIC=NO
acv_mean: MEAN ARG=acv PERIODIC=NO
DUMPATOMS ATOMS=1-64 ARG=acv FILE=mcolv.xyz
PRINT ARG=acv_sum,acv_mean FILE=colvar
```

## Deprecated syntax

More information on the deprecated keywords that are given below is available in the documentation for the [DISTANCES](DISTANCES.md) command.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace symfunc {

class AngularTetra : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit AngularTetra(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(AngularTetra,"TETRA_ANGULAR")

void AngularTetra::registerKeywords( Keywords& keys ) {
  CoordinationNumbers::shortcutKeywords( keys );
  keys.remove("MASK");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.add("compulsory","CUTOFF","-1","ignore distances that have a value larger than this cutoff");
  keys.setValueDescription("vector","the value of the angular tetehedrality parameter for each of the input atoms");
  keys.remove("NN");
  keys.remove("MM");
  keys.remove("D_0");
  keys.remove("R_0");
  keys.remove("SWITCH");
  keys.needsAction("DISTANCE_MATRIX");
  keys.needsAction("NEIGHBORS");
  keys.needsAction("GSYMFUNC_THREEBODY");
  keys.needsAction("CUSTOM");
}

AngularTetra::AngularTetra( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  bool nopbc;
  parseFlag("NOPBC",nopbc);
  std::string pbcstr="";
  if( nopbc ) {
    pbcstr = " NOPBC";
  }
  // Read species input and create the matrix
  std::string sp_str, rcut;
  parse("SPECIES",sp_str);
  parse("CUTOFF",rcut);
  if( sp_str.length()>0 ) {
    readInputLine( getShortcutLabel() + "_mat: DISTANCE_MATRIX COMPONENTS GROUP=" + sp_str + " CUTOFF=" + rcut + pbcstr );
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
    readInputLine( getShortcutLabel() + "_mat: DISTANCE_MATRIX COMPONENTS GROUPA=" + specA + " GROUPB=" + specB + " CUTOFF=" + rcut + pbcstr );
  }
  // Get the neighbors matrix
  readInputLine( getShortcutLabel() + "_neigh: NEIGHBORS ARG=" + getShortcutLabel() + "_mat.w NLOWEST=4");
  // Now construct the symmetry function (sum of cos(a) + 1/3)
  readInputLine( getShortcutLabel() + "_g8: GSYMFUNC_THREEBODY WEIGHT=" + getShortcutLabel() + "_neigh " +
                 "ARG=" + getShortcutLabel() + "_mat.x," + getShortcutLabel() + "_mat.y," + getShortcutLabel() + "_mat.z FUNCTION1={FUNC=(cos(ajik)+1/3)^2 LABEL=g8}");
  // Now evaluate the actual per atom CV
  readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_g8.g8 FUNC=(1-(3*x/8)) PERIODIC=NO");
  // And get the things to do with the quantities we have computed
  std::map<std::string,std::string> keymap;
  multicolvar::MultiColvarShortcuts::readShortcutKeywords( keymap, this );
  multicolvar::MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel(), "", keymap, this );
}

}
}
