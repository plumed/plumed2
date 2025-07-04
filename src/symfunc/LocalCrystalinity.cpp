/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionWithValue.h"
#include "CoordinationNumbers.h"
#include "multicolvar/MultiColvarShortcuts.h"

#include <complex>

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVAR LOCAL_CRYSTALINITY
/*
Calculate the local crystalinity symmetry function

This shortcut provides an implementation of the local crystalinity order parameter that is described in the paper in the bibliography.
To use this CV you define a series of unit vectors, $g_k$, using multiple instances of the GVECTOR keyword.  This allows you to define a
[symmetry function](https://www.plumed-tutorials.org/lessons/23/001/data/SymmetryFunction.html) for the $i$th atom as:

$$
s_i = \sum_k \left| \frac{\sum_j \sigma(|r_{ij}|) e^{ig_k r_{ij}}}{\sum_j \sigma(|r_{ij}|)} \right|^2
$$

In this expression $r_{ij}$ is the vector connecting atom $i$ to atom $j$ and $\sigma$ is a switching function that acts upon the magnidue of this vector, $|r_{ij}|$.
The following input is an example that shows how this symmetry function can be used in practice.

```plumed
lc: LOCAL_CRYSTALINITY ...
   SPECIES=1-64 D_0=3.0 R_0=1.5
   GVECTOR1=1,1,1 GVECTOR2=1,0.5,0.5 GVECTOR3=0.5,1.0,1.0 SUM
...
PRINT ARG=lc_sum FILE=colvar
```

As you can see if you expand the shortcut in this input, the sum over $k$ in the above expression has three terms in this input as 3 GVECTORS are specified.
Sixty four values for the expression above are computed.  These sixty four numbers are then added together in order to give a global mesuare of the crystallinity
for the simulated system.

In the input above $\sigma(|r_{ij}|)$ is a rational [switching function](LESS_THAN.md)
with the parameters above. We would recommend using SWITCH syntax
rather than the syntax above when giving the parameters for the switching function as you can then use any of the switching functions described
in the documentation for [LESS_THAN](LESS_THAN.md).  More importantly, however, using this syntax allows you to set the D_MAX parameter for the
switching function as demonstrated below:

```plumed
lc: LOCAL_CRYSTALINITY ...
   SPECIES=1-64 SWITCH={RATIONAL D_0=3.0 R_0=1.5 D_MAX=3.0}
   GVECTOR1=1,1,1 GVECTOR2=1,0.5,0.5 GVECTOR3=0.5,1.0,1.0 SUM
...
PRINT ARG=lc_sum FILE=colvar
```

## Working with two types of atom

If you would like to calculate the function above for the bonds connecting the atoms in GROUPA to the atoms in GROUPB you can use an input like the one
shown below:

```plumed
lc: LOCAL_CRYSTALINITY ...
   SPECIESA=1-5 SPECIESB=1-64 SWITCH={RATIONAL D_0=3.0 R_0=1.5 D_MAX=3.0}
   GVECTOR1=1,1,1 GVECTOR2=1,0.5,0.5 GVECTOR3=0.5,1.0,1.0 SUM
...
PRINT ARG=lc_sum FILE=colvar
```

## The MASK keyword

You can use the MASK keyword with this action in the same way that it is used with [COORDINATIONNUMBER](COORDINATIONNUMBER.md).  This keyword thus expects a vector in
input, which tells PLUMED the atoms for which you do not need to calculate the function.  As illustrated below, this is useful if you are using functionality
from the [volumes module](module_volumes.md) to calculate the average value of the TETRAHEDRAL parameter for only those atoms that lie in a certain part of the simulation box.

```plumed
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=1-64 CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# Calculate the local crystallinity parameters
lc: LOCAL_CRYSTALINITY ...
   SPECIES=1-64 SWITCH={RATIONAL D_0=3.0 R_0=1.5 D_MAX=3.0}
   GVECTOR1=1,1,1 GVECTOR2=1,0.5,0.5 GVECTOR3=0.5,1.0,1.0 MASK=sphere
...
# Multiply local crystallinity parameters numbers by sphere vector
prod: CUSTOM ARG=lc,sphere FUNC=x*y PERIODIC=NO
# Sum of coordination numbers for atoms that are in the sphere of interest
numer: SUM ARG=prod PERIODIC=NO
# Number of atoms that are in sphere of interest
denom: SUM ARG=sphere PERIODIC=NO
# Average coordination number for atoms in sphere of interest
av: CUSTOM ARG=numer,denom FUNC=x/y PERIODIC=NO
# And print out final CV to a file
PRINT ARG=av FILE=colvar STRIDE=1
```

This input calculate the average local crystallinity parameter for only those atoms that are within a spherical region that is centered on the point
$(2.5,2.5,2.5)$.

*/
//+ENDPLUMEDOC


class LocalCrystallinity : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit LocalCrystallinity(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(LocalCrystallinity,"LOCAL_CRYSTALINITY")

void LocalCrystallinity::registerKeywords( Keywords& keys ) {
  CoordinationNumbers::shortcutKeywords( keys );
  keys.add("numbered","GVECTOR","the coefficients of the linear combinations to compute for the CV");
  keys.setValueDescription("vector","the value of the local crystalinity for each of the input atoms");
  keys.needsAction("ONES");
  keys.needsAction("MATRIX_VECTOR_PRODUCT");
  keys.needsAction("COMBINE");
  keys.needsAction("CUSTOM");
  keys.addDOI("10.1063/1.4822877");
}

LocalCrystallinity::LocalCrystallinity( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // This builds an adjacency matrix
  std::string sp_str, specA, specB;
  parse("SPECIES",sp_str);
  parse("SPECIESA",specA);
  parse("SPECIESB",specB);
  CoordinationNumbers::expandMatrix( true, getShortcutLabel(), sp_str, specA, specB, this );
  // Input for denominator (coord)
  ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_mat");
  plumed_assert( av && av->getNumberOfComponents()>0 && (av->copyOutput(0))->getRank()==2 );
  std::string size;
  Tools::convert( (av->copyOutput(0))->getShape()[1], size );
  readInputLine( getShortcutLabel() + "_ones: ONES SIZE=" + size );
  readInputLine( getShortcutLabel() + "_denom: MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_mat.w," + getShortcutLabel() + "_ones");
  // Input for numerator
  std::string finput = "";
  for(unsigned i=1;; ++i) {
    std::vector<std::string> gvec;
    std::string istr;
    Tools::convert( i, istr );
    if( !parseNumberedVector("GVECTOR",i,gvec) ) {
      break;
    }
    if( gvec.size()!=3 ) {
      error("gvectors should have size 3");
    }
    // This is the dot product between the input gvector and the bond vector
    readInputLine( getShortcutLabel() + "_dot" + istr + ": COMBINE ARG=" + getShortcutLabel() + "_mat.x," + getShortcutLabel() + "_mat.y," + getShortcutLabel() + "_mat.z "
                   "PERIODIC=NO COEFFICIENTS=" + gvec[0] + "," + gvec[1] + "," + gvec[2] );
    // Now calculate the sine and cosine of the dot product
    readInputLine( getShortcutLabel() + "_cos" + istr + ": CUSTOM ARG=" + getShortcutLabel() +"_mat.w," + getShortcutLabel() + "_dot" + istr + " FUNC=x*cos(y) PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_sin" + istr + ": CUSTOM ARG=" + getShortcutLabel() +"_mat.w," + getShortcutLabel() + "_dot" + istr + " FUNC=x*sin(y) PERIODIC=NO");
    // And sum up the sine and cosine over the coordination spheres
    readInputLine( getShortcutLabel() + "_cossum" + istr + ": MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_cos" + istr + "," + getShortcutLabel() + "_ones");
    readInputLine( getShortcutLabel() + "_sinsum" + istr + ": MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_sin" + istr + "," + getShortcutLabel() + "_ones");
    // And average the sine and cosine over the number of bonds
    readInputLine( getShortcutLabel() + "_cosmean" + istr + ": CUSTOM FUNC=x/y PERIODIC=NO ARG=" + getShortcutLabel() + "_cossum" + istr + "," + getShortcutLabel() + "_denom");
    readInputLine( getShortcutLabel() + "_sinmean" + istr + ": CUSTOM FUNC=x/y PERIODIC=NO ARG=" + getShortcutLabel() + "_sinsum" + istr + "," + getShortcutLabel() + "_denom");
    // And work out the square modulus of this complex number
    readInputLine( getShortcutLabel() + "_" + istr + ": CUSTOM FUNC=x*x+y*y PERIODIC=NO ARG=" + getShortcutLabel() + "_cosmean" + istr + "," + getShortcutLabel() + "_sinmean" + istr);
    // These are all the kvectors that we are adding together in the final combine for the final CV
    if( i>1 ) {
      finput += ",";
    }
    finput += getShortcutLabel() + "_" + istr;
  }
  // This computes the final CV
  readInputLine( getShortcutLabel() + ": COMBINE NORMALIZE PERIODIC=NO ARG=" + finput );
  // Now calculate the total length of the vector
  multicolvar::MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel(), "", this );
}

}
}

