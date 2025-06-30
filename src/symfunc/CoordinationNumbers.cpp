/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include "CoordinationNumbers.h"
#include "multicolvar/MultiColvarShortcuts.h"
#include "core/ActionWithValue.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

#include <string>
#include <cmath>

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVAR COORDINATIONNUMBER
/*
Calculate the coordination numbers of atoms so that you can then calculate functions of the distribution of
 coordination numbers such as the minimum, the number less than a certain quantity and so on.

The coordination number of a atom $i$ is the number of atoms that are within a certain cutoff distance of it.
This quantity can be calculated as follows:

$$
s_i = \sum_j \sigma(r_{ij})
$$

where $r_{ij}$ is the distance between atoms $i$ and $j$ and $\sigma$ is a switching function. The typical switching
function that is used in metadynamics is this one:

$$
s(r) = \frac{ 1 - \left(\frac{r-d_0}{r_0}\right)^n } { 1 - \left(\frac{r-d_0}{r_0}\right)^m }
$$

The following example shows how you can use this shortcut action to calculate and print the coordination numbers of
one hundred atoms with themselves:

```plumed
c: COORDINATIONNUMBER SPECIES=1-100 R_0=1.0
DUMPATOMS ATOMS=c ARG=c FILE=coords.xyz
```

This input will produce an output file called coords that contains the coordination numbers of the 100 input atoms.  The cutoff
that is used to calculate the coordination number in this case is 1.0.  In the input above we use a rational [switching function](LESS_THAN.md)
with the parameters above. We would recommend using SWITCH syntax
rather than the syntax above when giving the parameters for the switching function as you can then use any of the switching functions described
in the documentation for [LESS_THAN](LESS_THAN.md).  More importantly, however, using this syntax allows you to set the D_MAX parameter for the
switching function as demonstrated below:

```plumed
c: COORDINATIONNUMBER SPECIES=1-100 SWITCH={RATIONAL R_0=1.0}
DUMPATOMS ATOMS=c ARG=c FILE=coords.xyz
```

Setting the `D_MAX` can substantially improve PLUMED performance as it turns on the linked list algorithm that is discussed in the optimisation details part
of the documentation for [CONTACT_MATRIX](CONTACT_MATRIX.md).

The vectors that are output by the COORDINATIONNUMBER shortcut can be used in the input for many other functions that are within
PLUMED. In addition, in order to ensure compatibility with older versions of PLUMED you can add additional keywords on the input
line for COORDINATIONNUMBER in order to calculate various functions of the input.  For example, the following input tells plumed
ato calculate the coordination numbers of atoms 1-100 with themselves.  The minimum coordination number is then calculated.

```plumed
c: COORDINATIONNUMBER SPECIES=1-100 R_0=1.0 MIN={BETA=0.1}
```

By constrast, this input tells plumed to calculate how many atoms from 1-100 are within 3.0 of each of the atoms
from 101-110.  In the first 101 is the central atom, in the second 102 is the central atom and so on.  The
number of coordination numbers that are more than 6 is then computed.

```plumed
c: COORDINATIONNUMBER SPECIESA=101-110 SPECIESB=1-100 R_0=3.0
mt: MORE_THAN ARG=c SWITCH={RATIONAL R_0=6.0 NN=6 MM=12 D_0=0}
s: SUM ARG=mt PERIODIC=NO
```

## The MASK keyword

Supppose that you want to calculate the average coordination number for the atoms that are within a sphere in the center of your simulation box. You can do so by exploiting an input similar to the one shown
below:

```plumed
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=1-400 CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# Calculate the coordination numbers
cc: COORDINATIONNUMBER SPECIES=1-400 SWITCH={RATIONAL D_0=3.0 R_0=1.5 D_MAX=6.0}
# Multiply fccubic parameters numbers by sphere vector
prod: CUSTOM ARG=cc,sphere FUNC=x*y PERIODIC=NO
# Sum of coordination numbers for atoms that are in the sphere of interest
numer: SUM ARG=prod PERIODIC=NO
# Number of atoms that are in sphere of interest
denom: SUM ARG=sphere PERIODIC=NO
# Average coordination number for atoms in sphere of interest
av: CUSTOM ARG=numer,denom FUNC=x/y PERIODIC=NO
# And print out final CV to a file
PRINT ARG=av FILE=colvar STRIDE=1
```

This calculation is slow because you have to calculate the coordination numbers of all the atoms even though only a small subset of these quanitties are required to compute the average coordination number in the
sphere.  To avoid all these unecessary calculations you use the `MASK` keyword as shown below:

```plumed
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=1-400 CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# Calculate the coordination numbers
cc: COORDINATIONNUMBER SPECIES=1-400 MASK=sphere SWITCH={RATIONAL D_0=3.0 R_0=1.5 D_MAX=6.0}
# Multiply fccubic parameters numbers by sphere vector
prod: CUSTOM ARG=cc,sphere FUNC=x*y PERIODIC=NO
# Sum of coordination numbers for atoms that are in the sphere of interest
numer: SUM ARG=prod PERIODIC=NO
# Number of atoms that are in sphere of interest
denom: SUM ARG=sphere PERIODIC=NO
# Average coordination number for atoms in sphere of interest
av: CUSTOM ARG=numer,denom FUNC=x/y PERIODIC=NO
# And print out final CV to a file
PRINT ARG=av FILE=colvar STRIDE=1
```

Adding the instruction `MASK=sphere` to the CONTACT_MATRIX line in this input tells PLUMED to only calculate the $i$th row in the adjacency matrix if the $i$th element of the vector `sphere` is non-zero.
In other words, by adding this command we have ensured that we are not calculating coordination numbers for atoms that are not in the sphere that is of interest.  In this way we can thus reduce the computational
expense of the calculation enormously.

Notice, that there are other places where we can use this same trick.  Some further examples are given in the documentation for [CONTACT_MATRIX](CONTACT_MATRIX.md).

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR COORDINATION_MOMENTS
/*
Calculate moments of the distribution of distances in the first coordination sphere

This is the CV that was developed by White and Voth and is described in the paper in the bibliograhy below. This action provides a way of indirectly biasing radial distribution functions and computes the following function

$$
s_i = \sum_j \sigma(r_{ij})r_{ij}^k
$$

where $k$ is the value that is input using the R_POWER keyword, $r_{ij}$ is the distance between atoms $i$ and $j$ and $\sigma$ is a switching function.

The following example shows how this action can be used.

```plumed
cn1: COORDINATION_MOMENTS SPECIES=1-10 R_0=1.0 R_POWER=1
cn1_mean: MEAN ARG=cn1 PERIODIC=NO
PRINT ARG=cn1_mean FILE=colvar
```

As you can see, the action works similarlly to [COORDINATIONNUMBER](COORDINATIONNUMBER.md).

In the input above we use a rational [switching function](LESS_THAN.md)
with the parameters above. We would recommend using SWITCH syntax
rather than the syntax above when giving the parameters for the switching function as you can then use any of the switching functions described
in the documentation for [LESS_THAN](LESS_THAN.md).  More importantly, however, using this syntax allows you to set the D_MAX parameter for the
switching function as demonstrated below:


```plumed
cn0: COORDINATIONNUMBER SPECIES=1-10 SWITCH={RATIONAL R_0=1.0 D_MAX=8}
cn0_mean: MEAN ARG=cn0 PERIODIC=NO
cn1: COORDINATION_MOMENTS SPECIES=1-10 SWITCH={RATIONAL R_0=1.0 D_MAX=8} R_POWER=1
cn1_mean: MEAN ARG=cn1 PERIODIC=NO
cn2: COORDINATION_MOMENTS SPECIES=1-10 SWITCH={RATIONAL R_0=1.0 D_MAX=8} R_POWER=2
cn2_mean: MEAN ARG=cn2 PERIODIC=NO
PRINT ARG=cn0_mean,cn1_mean,cn2_mean STRIDE=1 FILE=cn_out
```

Setting the `D_MAX` can substantially improve PLUMED performance as it turns on the linked list algorithm that is discussed in the optimisation details part
of the documentation for [CONTACT_MATRIX](CONTACT_MATRIX.md).

## Working with two types of atom

If you would like a way of indirectly biasing the radial distribution function that describes how the atoms in GROUPB are arranged around the atoms in GROUPA you use an input like the one
shown below:

```plumed
d: COORDINATION_MOMENTS ...
   SPECIESA=1-64 SPECIESB=65-200 R_POWER=1
   SWITCH={RATIONAL D_0=3.0 R_0=1.5 D_MAX=6.0}
...
s: MEAN ARG=d PERIODIC=NO
PRINT ARG=s FILE=colv
```

## The MASK keyword

You can use the MASK keyword with this action in the same way that it is used with [COORDINATIONNUMBER](COORDINATIONNUMBER.md).  This keyword thus expects a vector in
input, which tells COORDINATION_MOMENTS that it is safe not to calculate the COORDINATION_MOMENTS parameter for some of the atoms.  As illustrated below, this is useful if you are using functionality
from the [volumes module](module_volumes.md) to calculate the average value of the COORDINATION_MOMENTS parameter for only those atoms that lie in a certain part of the simulation box.

```plumed
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=1-400 CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# Calculate the coordination moments of the atoms
cc: COORDINATION_MOMENTS ...
  SPECIES=1-400 MASK=sphere R_POWER=1
  SWITCH={RATIONAL D_0=3.0 R_0=1.5 D_MAX=6.0}
...
# Multiply coordination moments by sphere vector
prod: CUSTOM ARG=cc,sphere FUNC=x*y PERIODIC=NO
# Sum of coordination numbers for atoms that are in the sphere of interest
numer: SUM ARG=prod PERIODIC=NO
# Number of atoms that are in sphere of interest
denom: SUM ARG=sphere PERIODIC=NO
# Average coordination number for atoms in sphere of interest
av: CUSTOM ARG=numer,denom FUNC=x/y PERIODIC=NO
# And print out final CV to a file
PRINT ARG=av FILE=colvar STRIDE=1
```

This input calculate the average value of the COORDINATION_MOMENTS parameter for only those atoms that are within a spherical region that is centered on the point
$(2.5,2.5,2.5)$.

## Deprecated syntax

More information on the deprecated keywords that are given below is available in the documentation for the [DISTANCES](DISTANCES.md) command.


*/
//+ENDPLUMEDOC


PLUMED_REGISTER_ACTION(CoordinationNumbers,"COORDINATIONNUMBER")
PLUMED_REGISTER_ACTION(CoordinationNumbers,"COORDINATION_MOMENTS")

void CoordinationNumbers::shortcutKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("atoms-3","SPECIES","the list of atoms for which the symmetry function is being calculated and the atoms that can be in the environments");
  keys.add("atoms-4","SPECIESA","the list of atoms for which the symmetry function is being calculated.  This keyword must be used in conjunction with SPECIESB, which specifies the atoms that are in the environment.");
  keys.add("atoms-4","SPECIESB","the list of atoms that can be in the environments of each of the atoms for which the symmetry function is being calculated.  This keyword must be used in conjunction with SPECIESA, which specifies the atoms for which the symmetry function is being calculated.");
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","the switching function that it used in the construction of the contact matrix");
  keys.add("optional","MASK","the label for a vector that is used to determine which rows of the matrix are computed");
  keys.linkActionInDocs("SWITCH","LESS_THAN");
  multicolvar::MultiColvarShortcuts::shortcutKeywords( keys );
  keys.needsAction("CONTACT_MATRIX");
  keys.needsAction("GROUP");
  keys.addDOI("10.1021/ct500320c");
}

void CoordinationNumbers::expandMatrix( const bool& components, const std::string& lab, const std::string& sp_str,
                                        const std::string& spa_str, const std::string& spb_str, ActionShortcut* action ) {
  if( sp_str.length()==0 && spa_str.length()==0 ) {
    return;
  }

  std::string matinp = lab  + "_mat: CONTACT_MATRIX";
  if( sp_str.length()>0 ) {
    matinp += " GROUP=" + sp_str;
    action->readInputLine( lab + "_grp: GROUP ATOMS=" + sp_str );
  } else if( spa_str.length()>0 ) {
    matinp += " GROUPA=" + spa_str + " GROUPB=" + spb_str;
    action->readInputLine( lab + "_grp: GROUP ATOMS=" + spa_str );
  }

  std::string sw_str;
  action->parse("SWITCH",sw_str);
  if( sw_str.length()>0 ) {
    matinp += " SWITCH={" + sw_str + "}";
  } else {
    std::string r0;
    action->parse("R_0",r0);
    std::string d0;
    action->parse("D_0",d0);
    if( r0.length()==0 ) {
      action->error("missing switching function parameters use SWITCH/R_0");
    }
    std::string nn;
    action->parse("NN",nn);
    std::string mm;
    action->parse("MM",mm);
    matinp += " R_0=" + r0 + " D_0=" + d0 + " NN=" + nn + " MM=" + mm;
  }
  if( components ) {
    matinp += " COMPONENTS";
  }
  std::string maskstr;
  action->parse("MASK",maskstr);
  if( maskstr.length()>0 ) {
    matinp += " MASK=" + maskstr;
  }
  action->readInputLine( matinp );
}

void CoordinationNumbers::registerKeywords( Keywords& keys ) {
  shortcutKeywords( keys );
  if( keys.getDisplayName()=="COORDINATION_MOMENTS" ) {
    keys.add("compulsory","R_POWER","the power to which you want to raise the distance");
  }
  keys.addDeprecatedFlag("LOWMEM","");
  keys.add("optional","MOMENTS","the list of moments that you would like to calculate");
  keys.addOutputComponent("moment","MOMENTS","scalar","the moments of the distribution");
  keys.reset_style("MOMENTS","deprecated");
  keys.needsAction("MATRIX_VECTOR_PRODUCT");
  keys.needsAction("ONES");
  keys.needsAction("MOMENTS");
  keys.setValueDescription("vector","the coordination numbers of the specified atoms");
}

CoordinationNumbers::CoordinationNumbers(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  bool lowmem;
  parseFlag("LOWMEM",lowmem);
  if( lowmem ) {
    warning("LOWMEM flag is deprecated and is no longer required for this action");
  }
  // Setup the contract matrix if that is what is needed
  std::string matlab, sp_str, specA, specB;
  parse("SPECIES",sp_str);
  parse("SPECIESA",specA);
  parse("SPECIESB",specB);
  if( sp_str.length()>0 || specA.length()>0 ) {
    matlab = getShortcutLabel() + "_mat";
    bool comp=false;
    if( getName()=="COORDINATION_MOMENTS" ) {
      comp=true;
      matlab = getShortcutLabel() + "_mat";
    }
    expandMatrix( comp, getShortcutLabel(), sp_str, specA, specB, this );
  } else {
    error("missing atoms input use SPECIES or SPECIESA/SPECIESB");
  }
  ActionWithValue* mb=plumed.getActionSet().selectWithLabel<ActionWithValue*>( matlab );
  if( !mb ) {
    error("could not find action with name " + matlab );
  }
  Value*  arg=mb->copyOutput(0);
  if( arg->getRank()!=2 || arg->hasDerivatives() ) {
    error("the input to this action should be a matrix or scalar");
  }
  // Create vector of ones to multiply input matrix by
  std::string nones;
  Tools::convert( arg->getShape()[1], nones );
  readInputLine( getShortcutLabel() + "_ones: ONES SIZE=" + nones );
  if( getName()=="COORDINATION_MOMENTS" ) {
    // Calculate the lengths of the vectors
    std::string r_power;
    parse("R_POWER",r_power);
    readInputLine( getShortcutLabel() + "_pow: CUSTOM ARG=" + matlab + ".x," + matlab + ".y," + matlab + ".z," + matlab + ".w VAR=x,y,z,w "
                   + "PERIODIC=NO FUNC=w*(sqrt(x*x+y*y+z*z)^" + r_power +")");
    matlab = getShortcutLabel() + "_pow";
  }
  // Calcualte coordination numbers as matrix vector times vector of ones
  readInputLine( getShortcutLabel() + ": MATRIX_VECTOR_PRODUCT  ARG=" + matlab + "," + getShortcutLabel() + "_ones");
  std::vector<std::string> moments;
  parseVector("MOMENTS",moments);
  Tools::interpretRanges( moments );
  if( moments.size()>0 ) {
    readInputLine( getShortcutLabel() + "_caverage: MEAN ARG=" + getShortcutLabel() + " PERIODIC=NO");
    for(unsigned i=0; i<moments.size(); ++i) {
      readInputLine( getShortcutLabel() + "_diffpow-" + moments[i] + ": CUSTOM ARG=" + getShortcutLabel() + "," + getShortcutLabel() + "_caverage PERIODIC=NO FUNC=(x-y)^" + moments[i] );
      readInputLine( getShortcutLabel() + "_moment-" + moments[i] + ": MEAN ARG=" + getShortcutLabel() + "_diffpow-" + moments[i] + " PERIODIC=NO");
    }
  }
  // Read in all the shortcut stuff
  std::map<std::string,std::string> keymap;
  multicolvar::MultiColvarShortcuts::readShortcutKeywords( keymap, this );
  multicolvar::MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel(), "", keymap, this );
}


}
}
