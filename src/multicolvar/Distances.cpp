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
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "MultiColvarShortcuts.h"
#include <string>
#include <cmath>

//+PLUMEDOC MCOLVAR DISTANCES
/*
Calculate the distances between multiple piars of atoms

__This shortcut action allows you to calculate function of the distribution of distances and reproduces the syntax in older PLUMED versions.
We think that the syntax in newer versions is more flexible and easier to use and that you should try to avoid using this old syntax.  If you look at the example inputs below you can
see how the new syntax can be used to replace this old one. We would strongly encourage you to use the newer syntax.__

The following input tells plumed to calculate the distances between atoms 3 and 5 and
between atoms 1 and 2 and to print the shorter of these two distances.

```plumed
d1: DISTANCES ATOMS1=3,5 ATOMS2=1,2 LOWEST
PRINT ARG=d1.lowest
```

The following input tells plumed to calculate the distances between atoms 3 and 5 and between atoms 1 and 2
and then to calculate the number of these distances that are less than 0.1 nm.  The number of distances
less than 0.1nm is then printed to a file.

```plumed
d1: DISTANCES ATOMS1=3,5 ATOMS2=1,2 LESS_THAN={RATIONAL R_0=0.1}
PRINT ARG=d1.lessthan
```

The following input tells plumed to calculate all 3 distances between atoms 1, 2 and 3 (i.e. the distances between atoms
1 and 2, atoms 1 and 3 and atoms 2 and 3).  The average of these distances is then calculated in two ways.

```plumed
d1: DISTANCES GROUP=1-3 MEAN SUM
d1mean: CUSTOM ARG=d1.sum FUNC=x/3 PERIODIC=NO
PRINT ARG=d1.mean,d1mean
```

The following input tells plumed to calculate all the distances between the atoms in GROUPA and the atoms in GROUPB.
In other words the distances between atoms 1 and 2 and the distance between atoms 1 and 3.  The number of distances
more than 0.1 is then printed to a file.

```plumed
d1: DISTANCES GROUPA=1 GROUPB=2,3 MORE_THAN={RATIONAL R_0=0.1}
PRINT ARG=d1.morethan
```

## Calculating minimum distances

There are various way for calculating the minimum distance between two groups of atoms.  The simplest one
is the one we saw above; namely:

```plumed
d1: DISTANCES GROUPA=1-10 GROUPB=11-20 LOWEST
PRINT ARG=d1.lowest
```

A disadvantage of this approach is that the derivatives of the resulting value are not continuous. To be clear,
the value of the shortest distance is continuous when you use this approach so you will likely be fine if you use
the quantity output by the command above as an input in a bias. However, if you would like to use a quantity with
continuous derivatives for the minimum distance you can use the following input:

```plumed
d1: DISTANCES GROUPA=1-10 GROUPB=11-20 MIN={BETA=500.}
PRINT ARG=d1.min FILE=colvar STRIDE=10
```

In order to ensure that the minimum value has continuous derivatives we use the following function:

$$
d_{min} = \frac{\beta}{ \log \sum_i \exp\left( \frac{\beta}{d_i} \right) }
$$

where $\beta$ is a user specified parameter.

This input is used rather than a separate MINDIST colvar so that the same routine and the same input style can be
used to calculate minimum coordination numbers (see [COORDINATIONNUMBER](COORDINATIONNUMBER.md)), minimum
angles (see [ANGLES](ANGLES.md)) and many other variables. If you expand the input above you will notice that
the DISTANCES command is a shortcut.  By expanding this shortcut you can see how PLUMED calculates the minimum distance.

An alternative implementation of the minimum distance is provided through the ALT_MIN keyword that is demonstrated below:

```plumed
d1: DISTANCES GROUPA=1-10 GROUPB=11-20 ALT_MIN={BETA=500.}
PRINT ARG=d1.altmin FILE=colvar STRIDE=10
```

The minimum value output by this command also has continuous derivatives and is evaluated using the following function:

$$
d_{min} = -\frac{1}{\beta} \log \left( \sum_i \exp(-\beta d_i) \right)
$$

## Calculating maximum distances

There are two ways for calculating the maximum distance between two groups of atom. The simplest one is:

```plumed
d1: DISTANCES GROUPA=1-10 GROUPB=11-20 HIGHEST
PRINT ARG=d1.highest
```

As with LOWEST above a disadvantage of this approach is that the derivatives of the resulting value are not continuous.
If you would like to use a quantity with continuous derivatives for the maximum distance you can use the following input:

```plumed
d1: DISTANCES GROUPA=1-10 GROUPB=11-20 MAX={BETA=500.}
PRINT ARG=d1.max
```

The maximum value output by this command is evaluated using the following function:

$$
d_{max} = \beta \log \left( \sum_i \exp\left( \frac{x}{\beta} \right) \right)
$$

## Other parts of the old multicolvar syntax

The ways of calculating the minimum and maximum distance described above are part of plumed 2's multicolvar functionality. In older
versions of PLUMED these special actions
allowed you to calculate multiple functions of a distribution of simple collective variables.  As an example you
can calculate the number of distances less than 1.0, the minimum distance, the number of distances more than
2.0, the number of distances between 1.0 and 2.0 and histogram that esimates the distribution of distances by using the following command:

```plumed
d1: DISTANCES ...
 GROUPA=1-10 GROUPB=11-20
 LESS_THAN={RATIONAL R_0=1.0}
 MORE_THAN={RATIONAL R_0=2.0}
 BETWEEN={GAUSSIAN LOWER=1.0 UPPER=2.0}
 HISTOGRAM={GAUSSIAN LOWER=0.0 UPPER=4.0 NBINS=10}
 MIN={BETA=500.}
...
PRINT ARG=d1.* FILE=colvar STRIDE=10
```

A calculation performed this way is fast because the expensive part of the calculation - the calculation of all the distances - is only done once per step.
If you expand the shortcut in the input above you can see how this result is achieved using a combination of many actions.  We recommend you use these combinations of
actions rather than the shortcut when using newer versions of PLUMED as you will have a more control over what you are doing.  This newer input syntax is easier to
understand especially when you start evaluating more complicated CVs.


*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR XDISTANCES
/*
Calculate the x components of the vectors connecting one or many pairs of atoms.

__This shortcut action allows you to calculate functions of the distribution of the x components of the vectors connecting pairs of atoms and
reproduces the syntax in older PLUMED versions.
If you look at the example inputs below you can
see how the new syntax operates. We would strongly encourage you to use the newer syntax as it offers greater flexibility.__

The following input tells plumed to calculate the x-component of the vector connecting atom 3 to atom 5 and
the x-component of the vector connecting atom 1 to atom 2.  The minimum of these two quantities is then
printed

```plumed
d1: XDISTANCES ATOMS1=3,5 ATOMS2=1,2 MIN={BETA=0.1}
PRINT ARG=d1.min
```

The following input tells plumed to calculate the x-component of the vector connecting atom 3 to atom 5 and
the x-component of the vector connecting atom 1 to atom 2.  The number of values that are
less than 0.1nm is then printed to a file.

```plumed
d1: XDISTANCES ATOMS1=3,5 ATOMS2=1,2 LESS_THAN={RATIONAL R_0=0.1}
PRINT ARG=d1.lessthan
```

The following input tells plumed to calculate the x-components of all the distinct vectors that can be created
between atoms 1, 2 and 3 (i.e. the vectors between atoms 1 and 2, atoms 1 and 3 and atoms 2 and 3).
The average of these quantities is then calculated.

```plumed
d1: XDISTANCES GROUP=1-3 MEAN
PRINT ARG=d1.mean
```

The following input tells plumed to calculate all the vectors connecting the the atoms in GROUPA to the atoms in GROUPB.
In other words the vector between atoms 1 and 2 and the vector between atoms 1 and 3.  The number of values
more than 0.1 is then printed to a file.

```plumed
d1: XDISTANCES GROUPA=1 GROUPB=2,3 MORE_THAN={RATIONAL R_0=0.1}
PRINT ARG=d1.morethan
```

*/
//+ENDPLUMEDOC


//+PLUMEDOC MCOLVAR YDISTANCES
/*
Calculate the y components of the vectors connecting one or many pairs of atoms.

__This shortcut action allows you to calculate functions of the distribution of the y components of the vectors connecting pairs of atoms and
reproduces the syntax in older PLUMED versions.
If you look at the example inputs below you can
see how the new syntax operates. We would strongly encourage you to use the newer syntax as it offers greater flexibility.__

The following input tells plumed to calculate the y-component of the vector connecting atom 3 to atom 5 and
the y-component of the vector connecting atom 1 to atom 2.  The minimum of these two quantities is then
printed

```plumed
d1: YDISTANCES ATOMS1=3,5 ATOMS2=1,2 MIN={BETA=0.1}
PRINT ARG=d1.min
```

The following input tells plumed to calculate the y-component of the vector connecting atom 3 to atom 5 and
the y-component of the vector connecting atom 1 to atom 2.  The number of values that are
less than 0.1nm is then printed to a file.

```plumed
d1: YDISTANCES ATOMS1=3,5 ATOMS2=1,2 LESS_THAN={RATIONAL R_0=0.1}
PRINT ARG=d1.lessthan
```

The following input tells plumed to calculate the y-components of all the distinct vectors that can be created
between atoms 1, 2 and 3 (i.e. the vectors between atoms 1 and 2, atoms 1 and 3 and atoms 2 and 3).
The average of these quantities is then calculated.

```plumed
d1: YDISTANCES GROUP=1-3 MEAN
PRINT ARG=d1.mean
```

The following input tells plumed to calculate all the vectors connecting the the atoms in GROUPA to the atoms in GROUPB.
In other words the vector between atoms 1 and 2 and the vector between atoms 1 and 3.  The number of values
more than 0.1 is then printed to a file.

```plumed
d1: YDISTANCES GROUPA=1 GROUPB=2,3 MORE_THAN={RATIONAL R_0=0.1}
PRINT ARG=d1.morethan
```

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR ZDISTANCES
/*
Calculate the z components of the vectors connecting one or many pairs of atoms.

__This shortcut action allows you to calculate functions of the distribution of the z components of the vectors connecting pairs of atoms and
reproduces the syntax in older PLUMED versions.
If you look at the example inputs below you can
see how the new syntax operates. We would strongly encourage you to use the newer syntax as it offers greater flexibility.__

The following input tells plumed to calculate the z-component of the vector connecting atom 3 to atom 5 and
the z-component of the vector connecting atom 1 to atom 2.  The minimum of these two quantities is then
printed

```plumed
d1: ZDISTANCES ATOMS1=3,5 ATOMS2=1,2 MIN={BETA=0.1}
PRINT ARG=d1.min
```

The following input tells plumed to calculate the z-component of the vector connecting atom 3 to atom 5 and
the z-component of the vector connecting atom 1 to atom 2.  The number of values that are
less than 0.1nm is then printed to a file.

```plumed
d1: ZDISTANCES ATOMS1=3,5 ATOMS2=1,2 LESS_THAN={RATIONAL R_0=0.1}
PRINT ARG=d1.lessthan
```

The following input tells plumed to calculate the z-components of all the distinct vectors that can be created
between atoms 1, 2 and 3 (i.e. the vectors between atoms 1 and 2, atoms 1 and 3 and atoms 2 and 3).
The average of these quantities is then calculated.

```plumed
d1: ZDISTANCES GROUP=1-3 MEAN
PRINT ARG=d1.mean
```

The following input tells plumed to calculate all the vectors connecting the the atoms in GROUPA to the atoms in GROUPB.
In other words the vector between atoms 1 and 2 and the vector between atoms 1 and 3.  The number of values
more than 0.1 is then printed to a file.

```plumed
d1: ZDISTANCES GROUPA=1 GROUPB=2,3 MORE_THAN={RATIONAL R_0=0.1}
PRINT ARG=d1.morethan
```

*/
//+ENDPLUMEDOC


namespace PLMD {
namespace multicolvar {

class Distances : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit Distances(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(Distances,"DISTANCES")
PLUMED_REGISTER_ACTION(Distances,"XDISTANCES")
PLUMED_REGISTER_ACTION(Distances,"YDISTANCES")
PLUMED_REGISTER_ACTION(Distances,"ZDISTANCES")

void Distances::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.setDeprecated("DISTANCE");
  keys.add("atoms-1","GROUP","Calculate the distance between each distinct pair of atoms in the group");
  keys.add("atoms-2","GROUPA","Calculate the distances between all the atoms in GROUPA and all "
           "the atoms in GROUPB. This must be used in conjunction with GROUPB.");
  keys.add("atoms-2","GROUPB","Calculate the distances between all the atoms in GROUPA and all the atoms "
           "in GROUPB. This must be used in conjunction with GROUPA.");
  keys.add("numbered","ATOMS","the pairs of atoms that you would like to calculate the angles for");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.addFlag("COMPONENTS",false,"calculate the x, y and z components of the distance separately and store them as label.x, label.y and label.z");
  keys.addFlag("SCALED_COMPONENTS",false,"calculate the a, b and c scaled components of the distance separately and store them as label.a, label.b and label.c");
  keys.addDeprecatedFlag("LOWMEM","");
  keys.reset_style("ATOMS","atoms");
  MultiColvarShortcuts::shortcutKeywords( keys );
  keys.add("atoms","ORIGIN","calculate the distance of all the atoms specified using the ATOMS keyword from this point");
  keys.add("numbered","LOCATION","the location at which the CV is assumed to be in space");
  keys.reset_style("LOCATION","atoms");
  keys.setValueDescription("vector","the DISTANCES between the each pair of atoms that were specified");
  keys.addOutputComponent("x","COMPONENTS","vector","the x-components of the distance vectors");
  keys.addOutputComponent("y","COMPONENTS","vector","the y-components of the distance vectors");
  keys.addOutputComponent("z","COMPONENTS","vector","the z-components of the distance vectors");
  keys.needsAction("GROUP");
  keys.needsAction("DISTANCE");
  keys.needsAction("CENTER");
}

Distances::Distances(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // Create distances
  bool lowmem;
  parseFlag("LOWMEM",lowmem);
  if( lowmem ) {
    warning("LOWMEM flag is deprecated and is no longer required for this action");
  }
  std::string dline = getShortcutLabel() + ": DISTANCE";
  bool nopbc;
  parseFlag("NOPBC",nopbc);
  if( nopbc ) {
    dline += " NOPBC";
  }
  if( getName()=="DISTANCES" ) {
    bool comp;
    parseFlag("COMPONENTS",comp);
    if( comp ) {
      dline += " COMPONENTS";
    }
    bool scomp;
    parseFlag("SCALED_COMPONENTS",scomp);
    if( scomp ) {
      dline += " SCALED_COMPONENTS";
    }
  } else {
    dline += " COMPONENTS";
  }
  // Parse origin
  std::string ostr;
  parse("ORIGIN",ostr);
  std::string num;
  if( ostr.length()>0 ) {
    // Parse atoms
    std::vector<std::string> afstr;
    MultiColvarShortcuts::parseAtomList("ATOMS",afstr,this);
    for(unsigned i=0; i<afstr.size(); ++i) {
      Tools::convert( i+1, num );
      dline += " ATOMS" + num + "=" + ostr + "," + afstr[i];
    }
  } else {
    std::vector<std::string> grp;
    MultiColvarShortcuts::parseAtomList("GROUP",grp,this);
    std::vector<std::string> grpa;
    MultiColvarShortcuts::parseAtomList("GROUPA",grpa,this);
    if( grp.size()>0 ) {
      if( grpa.size()>0 ) {
        error("should not be using GROUPA in tandem with GROUP");
      }
      unsigned n=0;
      for(unsigned i=1; i<grp.size(); ++i) {
        for(unsigned j=0; j<i; ++j) {
          Tools::convert( n+1, num );
          n++;
          dline += " ATOMS" + num + "=" + grp[i] + "," + grp[j];
        }
      }
    } else if( grpa.size()>0 ) {
      std::vector<std::string> grpb;
      MultiColvarShortcuts::parseAtomList("GROUPB",grpb,this);
      if( grpb.size()==0 ) {
        error("found GROUPA but no corresponding GROUPB");
      }
      for(unsigned i=0; i<grpa.size(); ++i) {
        for(unsigned j=0; j<grpb.size(); ++j) {
          Tools::convert( i*grpb.size() + j + 1, num );
          dline += " ATOMS" + num + "=" + grpa[i] + "," + grpb[j];
        }
      }
    } else {
      std::string grpstr = getShortcutLabel() + "_grp: GROUP ATOMS=";
      for(unsigned i=1;; ++i) {
        std::string atstring;
        parseNumbered("ATOMS",i,atstring);
        if( atstring.length()==0 ) {
          break;
        }
        std::string locstr;
        parseNumbered("LOCATION",i,locstr);
        if( locstr.length()==0 ) {
          Tools::convert( i, num );
          readInputLine( getShortcutLabel() + "_vatom" + num + ": CENTER ATOMS=" + atstring );
          if( i==1 ) {
            grpstr += getShortcutLabel() + "_vatom" + num;
          } else {
            grpstr += "," + getShortcutLabel() + "_vatom" + num;
          }
        } else {
          if( i==1 ) {
            grpstr += locstr;
          } else {
            grpstr += "," + locstr;
          }
        }
        Tools::convert( i, num );
        dline += " ATOMS" + num + "=" + atstring;
      }
      readInputLine( grpstr );
    }
  }
  readInputLine( dline );
  // Add shortcuts to label
  if( getName()=="DISTANCES" ) {
    MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel(), "", this );
  } else if( getName()=="XDISTANCES" ) {
    MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel() + ".x", "", this );
  } else if( getName()=="YDISTANCES" ) {
    MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel() + ".y", "", this );
  } else if( getName()=="ZDISTANCES" ) {
    MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel() + ".z", "", this );
  }
}

}
}
