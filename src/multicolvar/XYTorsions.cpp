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
#include "MultiColvarShortcuts.h"

//+PLUMEDOC COLVAR XYTORSIONS
/*
Calculate the torsional angle around the x axis between an arbitrary vector and the positive y direction

__As you can see if you expand the inputs below, you can achieve what this shortcut action does by using [TORSION](TORSION.md) together with [CUSTOM](CUSTOM.md),
[BETWEEN](BETWEEN.md), [LESS_THAN](LESS_THAN.md), [SUM](SUM.md) and [MEAN](MEAN.md).  We strongly encourage you to use these actions instead as using them will provide
you with a clearer understanding of the equations you are using.__

The following input tells plumed to calculate the angle around the x direction between the positive y-axis and the vector connecting atom 3 to atom 5 and
the angle around the x direction between the positive y axis and the vector connecting atom 1 to atom 2.  The average of these two quantities is then output

```plumed
d1: XYTORSIONS ATOMS1=3,5 ATOMS2=1,2 MEAN
PRINT ARG=d1_mean
```

Notice that this command is a shortcut. You can thus learn more about how to use PLUMED by examining the expanded version of the input above.

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR XZTORSIONS
/*
Calculate the torsional angle around the x axis between an arbitrary vector and the positive z direction

__As you can see if you expand the inputs below, you can achieve what this shortcut action does by using [TORSION](TORSION.md) together with [CUSTOM](CUSTOM.md),
[BETWEEN](BETWEEN.md), [LESS_THAN](LESS_THAN.md), [SUM](SUM.md) and [MEAN](MEAN.md).  We strongly encourage you to use these actions instead as using them will provide
you with a clearer understanding of the equations you are using.__

The following input tells plumed to calculate the angle around the x direction between the positive z-axis and the vector connecting atom 3 to atom 5 and
the angle around the x direction between the positive z axis and the vector connecting atom 1 to atom 2.  The average of these two quantities is then output

```plumed
d1: XZTORSIONS ATOMS1=3,5 ATOMS2=1,2 MEAN
PRINT ARG=d1_mean
```

Notice that this command is a shortcut. You can thus learn more about how to use PLUMED by examining the expanded version of the input above.

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR YXTORSIONS
/*
Calculate the torsional angle around the y axis between an arbitrary vector and the positive x direction

__As you can see if you expand the inputs below, you can achieve what this shortcut action does by using [TORSION](TORSION.md) together with [CUSTOM](CUSTOM.md),
[BETWEEN](BETWEEN.md), [LESS_THAN](LESS_THAN.md), [SUM](SUM.md) and [MEAN](MEAN.md).  We strongly encourage you to use these actions instead as using them will provide
you with a clearer understanding of the equations you are using.__

The following input tells plumed to calculate the angle around the y direction between the positive x-axis and the vector connecting atom 3 to atom 5 and
the angle around the y direction between the positive x axis and the vector connecting atom 1 to atom 2.  The average of these two quantities is then output

```plumed
d1: YXTORSIONS ATOMS1=3,5 ATOMS2=1,2 MEAN
PRINT ARG=d1_mean
```

Notice that this command is a shortcut. You can thus learn more about how to use PLUMED by examining the expanded version of the input above.

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR YZTORSIONS
/*
Calculate the torsional angle around the y axis between an arbitrary vector and the positive z direction

__As you can see if you expand the inputs below, you can achieve what this shortcut action does by using [TORSION](TORSION.md) together with [CUSTOM](CUSTOM.md),
[BETWEEN](BETWEEN.md), [LESS_THAN](LESS_THAN.md), [SUM](SUM.md) and [MEAN](MEAN.md).  We strongly encourage you to use these actions instead as using them will provide
you with a clearer understanding of the equations you are using.__

The following input tells plumed to calculate the angle around the y direction between the positive z-axis and the vector connecting atom 3 to atom 5 and
the angle around the y direction between the positive z axis and the vector connecting atom 1 to atom 2.  The average of these two quantities is then output

```plumed
d1: YZTORSIONS ATOMS1=3,5 ATOMS2=1,2 MEAN
PRINT ARG=d1_mean
```

Notice that this command is a shortcut. You can thus learn more about how to use PLUMED by examining the expanded version of the input above.

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR ZXTORSIONS
/*
Calculate the torsional angle around the z axis between an arbitrary vector and the positive x direction

__As you can see if you expand the inputs below, you can achieve what this shortcut action does by using [TORSION](TORSION.md) together with [CUSTOM](CUSTOM.md),
[BETWEEN](BETWEEN.md), [LESS_THAN](LESS_THAN.md), [SUM](SUM.md) and [MEAN](MEAN.md).  We strongly encourage you to use these actions instead as using them will provide
you with a clearer understanding of the equations you are using.__

The following input tells plumed to calculate the angle around the z direction between the positive x-axis and the vector connecting atom 3 to atom 5 and
the angle around the z direction between the positive x axis and the vector connecting atom 1 to atom 2.  The average of these two quantities is then output

```plumed
d1: ZXTORSIONS ATOMS1=3,5 ATOMS2=1,2 MEAN
PRINT ARG=d1_mean
```

Notice that this command is a shortcut. You can thus learn more about how to use PLUMED by examining the expanded version of the input above.

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR ZYTORSIONS
/*
Calculate the torsional angle around the z axis between an arbitrary vector and the positive y direction

__As you can see if you expand the inputs below, you can achieve what this shortcut action does by using [TORSION](TORSION.md) together with [CUSTOM](CUSTOM.md),
[BETWEEN](BETWEEN.md), [LESS_THAN](LESS_THAN.md), [SUM](SUM.md) and [MEAN](MEAN.md).  We strongly encourage you to use these actions instead as using them will provide
you with a clearer understanding of the equations you are using.__

The following input tells plumed to calculate the angle around the z direction between the positive y-axis and the vector connecting atom 3 to atom 5 and
the angle around the z direction between the positive y-axis and the vector connecting atom 1 to atom 2.  The average of these two quantities is then output

```plumed
d1: ZYTORSIONS ATOMS1=3,5 ATOMS2=1,2 MEAN
PRINT ARG=d1_mean
```

Notice that this command is a shortcut. You can thus learn more about how to use PLUMED by examining the expanded version of the input above.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

class XYTorsions : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit XYTorsions(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(XYTorsions,"XYTORSIONS")
PLUMED_REGISTER_ACTION(XYTorsions,"XZTORSIONS")
PLUMED_REGISTER_ACTION(XYTorsions,"YXTORSIONS")
PLUMED_REGISTER_ACTION(XYTorsions,"YZTORSIONS")
PLUMED_REGISTER_ACTION(XYTorsions,"ZXTORSIONS")
PLUMED_REGISTER_ACTION(XYTorsions,"ZYTORSIONS")

void XYTorsions::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("numbered","ATOMS","the pairs of atoms that you would like to calculate the angles for");
  keys.reset_style("ATOMS","atoms");
  MultiColvarShortcuts::shortcutKeywords( keys );
  keys.setValueDescription("vector","the angle between the vector connecting each pair of atoms and the the positive X/Y/Z direction around the X/Y/Z axis");
  keys.needsAction("FIXEDATOM");
  keys.needsAction("TORSION");
  keys.setDeprecated("TORSION");
}

XYTorsions::XYTorsions(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  std::string vdir = getShortcutLabel() + "_vec2," + getShortcutLabel() + "_origin";
  std::string adir = getShortcutLabel() + "_axis," + getShortcutLabel() + "_origin";
  // Create action for position of origin
  readInputLine( getShortcutLabel() + "_origin: FIXEDATOM AT=0,0,0");
  if( getName()=="XYTORSIONS" ) {
    readInputLine( getShortcutLabel() + "_axis: FIXEDATOM AT=1,0,0");
    readInputLine( getShortcutLabel() + "_vec2: FIXEDATOM AT=0,1,0");
  }
  if( getName()=="XZTORSIONS" ) {
    readInputLine( getShortcutLabel() + "_axis: FIXEDATOM AT=1,0,0");
    readInputLine( getShortcutLabel() + "_vec2: FIXEDATOM AT=0,0,1");
  }
  if( getName()=="YXTORSIONS" ) {
    readInputLine( getShortcutLabel() + "_axis: FIXEDATOM AT=0,1,0");
    readInputLine( getShortcutLabel() + "_vec2: FIXEDATOM AT=1,0,0");
  }
  if( getName()=="YZTORSIONS" ) {
    readInputLine( getShortcutLabel() + "_axis: FIXEDATOM AT=0,1,0");
    readInputLine( getShortcutLabel() + "_vec2: FIXEDATOM AT=0,0,1");
  }
  if( getName()=="ZXTORSIONS" ) {
    readInputLine( getShortcutLabel() + "_axis: FIXEDATOM AT=0,0,1");
    readInputLine( getShortcutLabel() + "_vec2: FIXEDATOM AT=1,0,0");
  }
  if( getName()=="ZYTORSIONS" ) {
    readInputLine( getShortcutLabel() + "_axis: FIXEDATOM AT=0,0,1");
    readInputLine( getShortcutLabel() + "_vec2: FIXEDATOM AT=0,1,0");
  }

  // Now create action to compute all torsions
  std::string torsions_str = getShortcutLabel() + ": TORSION";
  for(unsigned i=1;; ++i) {
    std::string atstring;
    parseNumbered("ATOMS",i,atstring);
    if( atstring.length()==0 ) {
      break;
    }
    std::string num;
    Tools::convert( i, num );
    torsions_str += " VECTORA" + num + "=" + atstring + " VECTORB" + num + "=" + vdir + " AXIS" + num + "=" + adir;
  }
  readInputLine( torsions_str );
  // Add shortcuts to label
  MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel(), "", this );
}

}
}
