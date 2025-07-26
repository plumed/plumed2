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
#include "core/ActionWithValue.h"
#include "multicolvar/MultiColvarShortcuts.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "CoordinationNumbers.h"
#include <string>
#include <cmath>

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVARF LOCAL_AVERAGE
/*
Calculate averages over spherical regions centered on atoms

As is explained in <a href="http://www.youtube.com/watch?v=iDvZmbWE5ps"> this video </a> certain PLUMED actions
calculate one scalar quantity or one vector for each of the atoms in the system.  For example
[COORDINATIONNUMBER](COORDINATIONNUMBER.md) measures the coordination number of each of the atoms in the system and [Q4](Q4.md) measures
the 4th order Steinhardt parameter for each of the atoms in the system.  These quantities provide tell us something about
the disposition of the atoms in the first coordination sphere of each of the atoms of interest.  In the paper in the bibliography Lechner and Dellago
have suggested that one can probe local order in a system by taking the average value of such symmetry functions over
the atoms within a spherical cutoff of each of these atoms in the systems.  When this is done with Steinhardt parameters
they claim this gives a coordinate that is better able to distinguish solid and liquid configurations of Lennard-Jones atoms.

You can calculate such locally averaged quantities within plumed by using the LOCAL_AVERAGE command.  This command calculates
the following atom-centered quantities:

$$
s_i = \frac{ c_i + \sum_j \sigma(r_{ij})c_j }{ 1 + \sum_j \sigma(r_{ij}) }
$$

where the $c_i$ and $c_j$ values can be any vector of [symmetry functions](https://www.plumed-tutorials.org/lessons/23/001/data/SymmetryFunction.html)
 that can be calculated using plumed multicolvars.  The function $\sigma( r_{ij} )$ is a switching function that acts on the distance between
atoms $i$ and $j$.  Lechner and Dellago suggest that the parameters of this function should be set so that it the function is equal to one
when atom $j$ is in the first coordination sphere of atom $i$ and is zero otherwise.

To see how this works in practice consider the following example input.

```plumed
d1: COORDINATIONNUMBER SPECIES=1-64 SWITCH={RATIONAL D_0=1.3 R_0=0.2 D_MAX=3.0}
la: LOCAL_AVERAGE SPECIES=d1 D_0=1.3 R_0=0.2 NN=6 MM=12
mtf: MORE_THAN ARG=la SWITCH={RATIONAL R_0=4}
la_morethan: SUM ARG=mtf PERIODIC=NO
PRINT ARG=la_morethan FILE=colvar
```

This input calculates the coordination numbers for all the atoms in the system.  These coordination numbers are then averaged over
spherical regions.  The number of averaged coordination numbers that are greater than 4 is then output to a file.  Furthermore, if you
expand the input above you can see how the LOCAL_AVERAGE command is a shortcut action that expands to a longer input that you should be able to
interpret.

In the input above we use a rational [switching function](LESS_THAN.md) with the parameters above. We would recommend using SWITCH syntax
rather than the syntax above when giving the parameters for the switching function as you can then use any of the switching functions described
in the documentation for [LESS_THAN](LESS_THAN.md).  More importantly, however, using this syntax allows you to set the D_MAX parameter for the
switching function as demonstrated below:

```plumed
d1: COORDINATIONNUMBER SPECIES=1-64 SWITCH={RATIONAL D_0=1.3 R_0=0.2 D_MAX=3.0}
la: LOCAL_AVERAGE SPECIES=d1 SWITCH={RATIONAL D_0=1.3 R_0=0.2 D_MAX=3.0}
mtf: MORE_THAN ARG=la SWITCH={RATIONAL R_0=4}
la_morethan: SUM ARG=mtf PERIODIC=NO
PRINT ARG=la_morethan FILE=colvar
```

Setting the `D_MAX` can substantially improve PLUMED performance as it turns on the linked list algorithm that is discussed in the optimisation details part
of the documentation for [CONTACT_MATRIX](CONTACT_MATRIX.md).

## Locally averaged Steinhardt Parameters

What Lechner and Dellago did in their paper was a little more complicated than the example in the previous section. To reproduce what they did you would use
an input something like this:

```plumed
q4: Q4 SPECIES=1-64 SWITCH={RATIONAL D_0=1.3 R_0=0.2}
la: LOCAL_AVERAGE SPECIES=q4 SWITCH={RATIONAL D_0=1.3 R_0=0.2}
la_mean: MEAN ARG=la PERIODIC=NO
PRINT ARG=la_mean FILE=colvar
```

This example input calculates the [Q4](Q4.md) vectors for each of the atoms in the system.  These vectors are then averaged
component by component over a spherical region.  The average value for this quantity is then outputeed to a file.  If you want
to understand more about the shortcut that is used here you can read [this page](https://www.plumed-tutorials.org/lessons/23/001/data/Steinhardt.html).

## Using SPECIESA/SPECIESB

If you only want to calculate the local averages for a subset of the atoms you might use an input like the one shown below:

```plumed
# Calculate the coordination numbers for the two types of atom
d1a: COORDINATIONNUMBER SPECIES=1-5 SWITCH={RATIONAL D_0=1.3 R_0=0.2 D_MAX=3.0}
d1b: COORDINATIONNUMBER SPECIES=1-64 SWITCH={RATIONAL D_0=1.3 R_0=0.2 D_MAX=3.0}
# Now compute the local average
la: LOCAL_AVERAGE SPECIESA=d1a SPECIESB=d1b SWITCH={RATIONAL D_0=1.3 R_0=0.2 D_MAX=3.0}
DUMPATOMS ARG=la ATOMS=1-5 FILE=atoms.xyz
```

In the input above the local averages are computed for atoms 1-5 only. However, to calculate these local averages we still need to calculate the coordination numbers of all
the atoms in the system.

## The MASK keyword

You can use the MASK keyword with this action in the same way that it is used with [COORDINATIONNUMBER](COORDINATIONNUMBER.md).  This keyword thus expects a vector in
input, which tells PLUMED the atoms for which you do not need to calculate the function.  As illustrated below, this is useful if you are using functionality
from the [volumes module](module_volumes.md) to calculate the average value of the local averages for only those atoms that
lie in a certain part of the simulation box.

```plumed
# Calculate the coordination number for all the atoms
cc: COORDINATIONNUMBER SPECIES=1-400 SWITCH={RATIONAL D_0=1.3 R_0=0.2 D_MAX=3.0}
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=1-400 CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# Calculate the local average
lcc: LOCAL_AVERAGE SPECIES=cc SWITCH={RATIONAL D_0=1.3 R_0=0.2 D_MAX=3.0} MASK=sphere
# Multiply local averages parameters numbers by sphere vector
prod: CUSTOM ARG=lcc,sphere FUNC=x*y PERIODIC=NO
# Sum of coordination numbers for atoms that are in the sphere of interest
numer: SUM ARG=prod PERIODIC=NO
# Number of atoms that are in sphere of interest
denom: SUM ARG=sphere PERIODIC=NO
# Average coordination number for atoms in sphere of interest
av: CUSTOM ARG=numer,denom FUNC=x/y PERIODIC=NO
# And print out final CV to a file
PRINT ARG=av FILE=colvar STRIDE=1
```

This input calculates the average value of the local average for those atoms that are within a spherical region that is centered on the point $(2.5,2.5,2.5)$.
By including the MASK keyword in the LOCAL_AVERAGE line we reduce the number of local averages that are computed. These quantities are only calculated for those
atoms that are within the spherical region of interest.  However, we are still asking PLUMED to calculate the coordination numbers for many atoms that
will not contribute to the final averaged quantity in the COORDINATIONNUMBER commmand with label `cc`.  If you use an input like the one shown below you can ensure
that only those coordiantion numbers that compute to the final average are computed.

```plumed
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=1-400 CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
ones: ONES SIZE=400
# Calculate a elements of adjacency matrix for those elements that are inside the sphere
vol_cmap: CONTACT_MATRIX GROUP=1-400 SWITCH={RATIONAL D_0=1.3 R_0=0.2 D_MAX=3.0} MASK=sphere
# Create a matrix that for the sphere stuff
sphere_mat: OUTER_PRODUCT ARG=sphere,ones MASK=vol_cmap
# Multiply this contact matrix by the map
bonds_mat: CUSTOM ARG=vol_cmap,sphere_mat FUNC=x*y PERIODIC=NO
# Transpose the above matrix
bonds_matT: TRANSPOSE ARG=bonds_mat
# And multiply by the volume
bonds: MATRIX_VECTOR_PRODUCT ARG=bonds_matT,sphere
#  And make sure we calculate the coordination numbers for the atoms in the volume too
mask: CUSTOM ARG=bonds,sphere FUNC=x+y PERIODIC=NO
# Calculate the coordination number for all the atoms
cc: COORDINATIONNUMBER SPECIES=1-400 SWITCH={RATIONAL D_0=1.3 R_0=0.2 D_MAX=3.0} MASK=mask
# Calculate the local average
lcc: LOCAL_AVERAGE SPECIES=cc SWITCH={RATIONAL D_0=1.3 R_0=0.2 D_MAX=3.0} MASK=sphere
# Multiply local averages parameters numbers by sphere vector
prod: CUSTOM ARG=lcc,sphere FUNC=x*y PERIODIC=NO
# Sum of coordination numbers for atoms that are in the sphere of interest
numer: SUM ARG=prod PERIODIC=NO
# Number of atoms that are in sphere of interest
denom: SUM ARG=sphere PERIODIC=NO
# Average coordination number for atoms in sphere of interest
av: CUSTOM ARG=numer,denom FUNC=x/y PERIODIC=NO
# And print out final CV to a file
PRINT ARG=av FILE=colvar STRIDE=1
```

In the above input the only elements of the [CONTACT_MATRIX](CONTACT_MATRIX.md) with label `vol_cmap` that are non-zero are those that are bonded to one of the atoms
in the sphere of interest.  By transposing this matrix and multiplying by a vector of ones we thus get a vector with a length equal to total number of coordination
numbers in the system. The only elements in this vector that are non-zero are those that are bonded to atoms that are within the sphere of interest.  Consequently,
if we add the vector `sphere` to this we can a vector in which element $i$ is only non-zero if we have to calculate that coordination number to evaluate the final
quantity of interest.

*/
//+ENDPLUMEDOC

class LocalAverage : public ActionShortcut {
private:
  std::string getMomentumSymbol( const int& m ) const ;
public:
  static void registerKeywords( Keywords& keys );
  explicit LocalAverage(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(LocalAverage,"LOCAL_AVERAGE")

void LocalAverage::registerKeywords( Keywords& keys ) {
  CoordinationNumbers::shortcutKeywords( keys );
  keys.needsAction("ONES");
  keys.needsAction("MATRIX_VECTOR_PRODUCT");
  keys.needsAction("VSTACK");
  keys.needsAction("CUSTOM");
  keys.needsAction("OUTER_PRODUCT");
  keys.setValueDescription("vector","the values of the local averages");
  keys.addDOI("10.1063/1.2977970");
  keys.addFlag("USEGPU",false,"run part of this calculation on the GPU");
}

LocalAverage::LocalAverage(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  std::string sp_str, specA, specB;
  parse("SPECIES",sp_str);
  parse("SPECIESA",specA);
  parse("SPECIESB",specB);

  bool usegpu;
  parseFlag("USEGPU",usegpu);
  const std::string doUSEGPU = usegpu?" USEGPU":"";

  CoordinationNumbers::expandMatrix( false,
                                     getShortcutLabel(),
                                     sp_str,
                                     specA,
                                     specB,
                                     this );
  std::map<std::string,std::string> keymap;
  multicolvar::MultiColvarShortcuts::readShortcutKeywords( keymap, this );
  if( sp_str.length()>0 ) {
    specA=specB=sp_str;
  }
  const std::string matLab=getShortcutLabel() + "_mat";
  const std::string onesLab=getShortcutLabel() + "_ones";
  const std::string coordLab=getShortcutLabel() + "_coord";
  // Calculate the coordination numbers
  ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>(matLab);
  plumed_assert( av && av->getNumberOfComponents()>0 && (av->copyOutput(0))->getRank()==2 );
  std::string size;
  Tools::convert( (av->copyOutput(0))->getShape()[1], size );
  readInputLine(onesLab + ": ONES SIZE=" + size );
  readInputLine(coordLab + ": MATRIX_VECTOR_PRODUCT ARG=" + matLab + "," + onesLab + doUSEGPU);

  int l=-1;
  std::vector<ActionShortcut*> shortcuts=plumed.getActionSet().select<ActionShortcut*>();
  for(unsigned i=0; i<shortcuts.size(); ++i) {
    if( specA==shortcuts[i]->getShortcutLabel() ) {
      std::string sname = shortcuts[i]->getName();
      if( sname=="Q1" || sname=="Q3" || sname=="Q4" || sname=="Q6" ) {
        Tools::convert( sname.substr(1), l );
        break;
      }
    }
  }

  const std::string prodLab=getShortcutLabel() + "_prod";
  if( l>0 ) {
    std::string vargs="";
    std::string svargs="";
    std::string comma="ARG=";
    std::string sargs = "ARG=" + matLab;
    auto stringDiv=[&](const std::string& num, char realImm) {
      return specB + "_"+realImm+"mn-" + num + ": CUSTOM ARG="
             + specB + "_sp."+realImm+"m-" + num + ","
             + specB + "_denom FUNC=x/y PERIODIC=NO";
    };
    for(int i=-l; i<=l; ++i) {
      std::string num = getMomentumSymbol(i);
      if( !plumed.getActionSet().selectWithLabel<ActionWithValue*>(specB + "_rmn-" + num) ) {
        readInputLine( stringDiv(num,'r') );
      }
      if( !plumed.getActionSet().selectWithLabel<ActionWithValue*>(specB + "_imn-" + num) ) {
        readInputLine( stringDiv(num,'i') );
      }
      vargs += comma +  specB + "_rmn-" + num + "," + specB + "_imn-" + num;
      svargs += comma + prodLab + "." + specB + "_rmn-" + num
                + "," + prodLab + "." + specB + "_imn-" + num;
      comma=",";
      sargs += "," + specB + "_rmn-" + num + "," + specB  + "_imn-" + num;
    }

    const std::string vstackLab=getShortcutLabel() + "_vstack";
    const std::string vpstackLab=getShortcutLabel() + "_vpstack";
    readInputLine(vstackLab + ": VSTACK " + vargs );
    readInputLine(prodLab + ": MATRIX_VECTOR_PRODUCT " + sargs + doUSEGPU);
    readInputLine(vpstackLab + ": VSTACK " + svargs );
    std::string twolplusone;
    Tools::convert( 2*(2*l+1), twolplusone );

    const std::string lonesLab=getShortcutLabel() + "_lones";
    readInputLine(lonesLab + ": ONES SIZE=" + twolplusone );

    const std::string unormLab=getShortcutLabel() + "_unorm";
    readInputLine(unormLab + ": OUTER_PRODUCT ARG=" +coordLab + "," + lonesLab );

    const std::string avLab=getShortcutLabel() + "_av";
    readInputLine(avLab + ": CUSTOM ARG=" + vpstackLab + "," + vstackLab + ","
                  + unormLab + " FUNC=(x+y)/(1+z) PERIODIC=NO");

    const std::string av2Lab=getShortcutLabel() + "_av2";
    readInputLine(av2Lab + ": CUSTOM ARG=" +avLab + " FUNC=x*x PERIODIC=NO");

    const std::string twoLab=getShortcutLabel() + "_2";
    readInputLine( twoLab+": MATRIX_VECTOR_PRODUCT ARG=" +av2Lab + "," + lonesLab
                   + doUSEGPU);
    readInputLine( getShortcutLabel() + ": CUSTOM "
                   "ARG=" + twoLab + " FUNC=sqrt(x) PERIODIC=NO");
  } else {
    readInputLine( prodLab + ": MATRIX_VECTOR_PRODUCT "
                   "ARG=" +matLab + ","+ specB
                   + doUSEGPU);
    readInputLine( getShortcutLabel() + ": CUSTOM "
                   "ARG=" + prodLab + "," + specA + "," +coordLab
                   + " FUNC=(x+y)/(1+z) PERIODIC=NO");
  }
  multicolvar::MultiColvarShortcuts::expandFunctions( getShortcutLabel(),
      getShortcutLabel(),
      "",
      keymap, this );
}

std::string LocalAverage::getMomentumSymbol( const int& m ) const {
  if( m<0 ) {
    std::string num;
    Tools::convert( -1*m, num );
    return "n" + num;
  } else if( m>0 ) {
    std::string num;
    Tools::convert( m, num );
    return "p" + num;
  }
  return "0";
}

}
}
