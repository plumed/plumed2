/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2023 The plumed team
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

//+PLUMEDOC ANALYSIS KL_ENTROPY
/*
Calculate the KL entropy of a distribution

This shortcut was implemented in order to make the implementation of CVs like those described in the paper that is cited below straightforward.
An example input that can be used to implement one of the CVs that is introduced in that paper is shown below:

```plumed
#SETTINGS INPUTFILES=regtest/gridtools/rt-weights-integral/kde.grid

# These two commands compute a set of weights that allow us to determine if each distance in our system is a bond
d1: DISTANCE ATOMS1=1,2 ATOMS2=1,3 ATOMS3=1,4 ATOMS4=1,5 ATOMS5=2,3 ATOMS6=2,4 ATOMS7=2,5 ATOMS8=3,4 ATOMS9=3,5 ATOMS10=4,5
d1lt: LESS_THAN ARG=d1 SWITCH={RATIONAL D_0=2.0 R_0=0.5 D_MAX=5.0}

# These commands caluclate the angle between the z axis and the bond directions
d1c: DISTANCE ATOMS1=2,1 ATOMS2=3,1 ATOMS3=4,1 ATOMS4=5,1 ATOMS5=3,2 ATOMS6=4,2 ATOMS7=5,2 ATOMS8=4,3 ATOMS9=5,3 ATOMS10=5,4 COMPONENTS
d2: COMBINE ARG=d1c.x,d1c.y,d1c.z POWERS=2,2,2 PERIODIC=NO
aa: MATHEVAL ARG=d1c.z,d2 FUNC=acos(x/sqrt(y)) PERIODIC=NO

# These commands compute the torsion angle between the bond direction and the positive x direction around the z axis
dd0: FIXEDATOM AT=0,0,0
ddx: FIXEDATOM AT=1,0,0
ddz: FIXEDATOM AT=0,0,1

tt: TORSION ...
   VECTORA1=2,1 VECTORB1=ddx,dd0 AXIS1=ddz,dd0
   VECTORA2=3,1 VECTORB2=ddx,dd0 AXIS2=ddz,dd0
   VECTORA3=4,1 VECTORB3=ddx,dd0 AXIS3=ddz,dd0
   VECTORA4=5,1 VECTORB4=ddx,dd0 AXIS4=ddz,dd0
   VECTORA5=3,2 VECTORB5=ddx,dd0 AXIS5=ddz,dd0
   VECTORA6=4,2 VECTORB6=ddx,dd0 AXIS6=ddz,dd0
   VECTORA7=5,2 VECTORB7=ddx,dd0 AXIS7=ddz,dd0
   VECTORA8=4,3 VECTORB8=ddx,dd0 AXIS8=ddz,dd0
   VECTORA9=5,3 VECTORB9=ddx,dd0 AXIS9=ddz,dd0
   VECTORA10=5,4 VECTORB10=ddx,dd0 AXIS10=ddz,dd0
...

# We now compute our instaneous normalised histogram of the bond directions.  Notice that the weights are used here so that we do not consider
# non-bonded atoms
hu: KDE VOLUMES=d1lt ARG=aa,tt GRID_BIN=20,20 GRID_MIN=0,-pi GRID_MAX=pi,pi BANDWIDTH=0.2,0.2
de: SUM ARG=d1lt PERIODIC=NO
h: CUSTOM ARG=hu,de FUNC=x/y PERIODIC=NO

# And lastly compute the KL_ENTROPY
kl: KL_ENTROPY ARG=h REFERENCE=regtest/gridtools/rt-weights-integral/kde.grid VALUE=h
PRINT ARG=kl FILE=colvar
```

The CV `kl` that is defined in the input above describes the disribution of bond directions that is observed in the structure.  To compute it we first define the orientations
of all the bonds in our structure relative to the lab frame by computing the angle between each bond direction and the z axis and the torsional angle around the z axis between
our bond vector and the positive x direction.  In other words, we compute two of the three [spherical polar coordinates](https://en.wikipedia.org/wiki/Spherical_coordinate_system)
for each of the bond directions in our crystal structure.

A normalised histogram that describes the instaneous distribution for these bond directions is then computed by using the [KDE](KDE.md) action. The final scalar value for the CV is computed by
determining the [Kullbeck-Leibler divergence](https://en.wikipedia.org/wiki/Kullback–Leibler_divergence) between this instaneous distribution and a reference distribution that is provided in input.

Notice that this action will also works if you use KDEs computed using the [SPHERICAL_KDE](SPHERICAL_KDE.md) action in input as is illustrated in the following input:

```plumed
#SETTINGS INPUTFILES=regtest/gridtools/rt-kldiv/averageSp_X240.grid

# These two commands compute a set of weights that allow us to determine if each distance in our system is a bond
d1: DISTANCE ATOMS1=1,2 ATOMS2=1,3 ATOMS3=1,4 ATOMS4=1,5 ATOMS5=2,3 ATOMS6=2,4 ATOMS7=2,5 ATOMS8=3,4 ATOMS9=3,5 ATOMS10=4,5
d1lt: LESS_THAN ARG=d1 SWITCH={RATIONAL D_0=2.0 R_0=0.5 D_MAX=5.0}

# These commands compute the directors of the bonds between the atoms
d1c: DISTANCE ATOMS1=2,1 ATOMS2=3,1 ATOMS3=4,1 ATOMS4=5,1 ATOMS5=3,2 ATOMS6=4,2 ATOMS7=5,2 ATOMS8=4,3 ATOMS9=5,3 ATOMS10=5,4 COMPONENTS
d1l: CUSTOM ARG=d1c.x,d1c.y,d1c.z FUNC=sqrt(x*x+y*y+z*z) PERIODIC=NO
d1nx: CUSTOM ARG=d1c.x,d1l FUNC=x/y PERIODIC=NO
d1ny: CUSTOM ARG=d1c.y,d1l FUNC=x/y PERIODIC=NO
d1nz: CUSTOM ARG=d1c.z,d1l FUNC=x/y PERIODIC=NO

# Now compute our instaneous normalised histogram of the bond directions
hu: SPHERICAL_KDE HEIGHTS=d1lt ARG=d1nx,d1ny,d1nz GRID_BIN=400 CONCENTRATION=100.0
de: SUM ARG=d1lt PERIODIC=NO
h: CUSTOM ARG=hu,de FUNC=x/y PERIODIC=NO

# And lastly compute the KL_ENTROPY
kl: KL_ENTROPY ARG=h REFERENCE=regtest/gridtools/rt-kldiv/averageSp_X240.grid VALUE=av_Sp_X
PRINT ARG=kl FILE=colvar
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

class KLEntropy : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit KLEntropy(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(KLEntropy,"KL_ENTROPY")

void KLEntropy::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords(keys);
  keys.add("compulsory","ARG","the grid that you wish to use in the KL entropy calculation");
  keys.add("compulsory","REFERENCE","a file containing the reference density in grid format");
  keys.add("compulsory","VALUE","the name of the value that should be read from the grid");
  keys.setValueDescription("scalar","the value of the KL-Entropy");
  keys.addDOI("10.1021/acs.jctc.7b01027");
  keys.needsAction("REFERENCE_GRID");
  keys.needsAction("CUSTOM");
  keys.needsAction("INTEGRATE_GRID");
}

KLEntropy::KLEntropy( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // Reference grid object
  std::string ref_str, val_str, input_g;
  parse("VALUE",val_str);
  parse("REFERENCE",ref_str);
  parse("ARG",input_g);
  readInputLine( getShortcutLabel() + "_ref: REFERENCE_GRID VALUE=" + val_str + " FILE=" + ref_str );
  // Compute KL divergence
  if( input_g=="") {
    plumed_merror("could not find ARG keyword in input to KL_ENTROPY");
  }
  readInputLine( getShortcutLabel()  + "_kl: CUSTOM ARG=" + input_g + "," + getShortcutLabel() + "_ref FUNC=y*log(y/(0.5*(x+y))) PERIODIC=NO");
  // Compute integral
  readInputLine( getShortcutLabel() + ": INTEGRATE_GRID ARG=" + getShortcutLabel() + "_kl PERIODIC=NO");
}

}
}
