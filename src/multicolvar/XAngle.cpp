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
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "MultiColvarShortcuts.h"

//+PLUMEDOC COLVAR XANGLES
/*
Calculate the angle between an arbitrary vector and the positive x direction

__As you can see if you expand the inputs below, you can achieve what this shortcut action does by using [DISTANCE](DISTANCE.md) together with [CUSTOM](CUSTOM.md),
[BETWEEN](BETWEEN.md), [LESS_THAN](LESS_THAN.md), [SUM](SUM.md) and [MEAN](MEAN.md).  We strongly encourage you to use these actions instead as using them will provide
you with a clearer understanding of the equations you are using.__

The following input tells plumed to calculate the angles between the x-axis and the vector connecting atom 3 to atom 5 and between the x-axis
and the vector connecting atom 1 to atom 2.  The minimum of these two quantities is then output

```plumed
d1: XANGLES ATOMS1=3,5 ATOMS2=1,2 MIN={BETA=0.1}
PRINT ARG=d1.min FILE=colvar
```

Notice that this command is a shortcut. You can thus learn more about how to use PLUMED by examining the expanded version of the input above.

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR YANGLES
/*
Calculate the angle between an arbitrary vector and the positive y direction

__As you can see if you expand the inputs below, you can achieve what this shortcut action does by using [DISTANCE](DISTANCE.md) together with [CUSTOM](CUSTOM.md),
[BETWEEN](BETWEEN.md), [LESS_THAN](LESS_THAN.md), [SUM](SUM.md) and [MEAN](MEAN.md).  We strongly encourage you to use these actions instead as using them will provide
you with a clearer understanding of the equations you are using.__

The following input tells plumed to calculate the angles between the y-axis and the vector connecting atom 3 to atom 5 and between the y-axis
and the vector connecting atom 1 to atom 2.  The minimum of these two quantities is then output

```plumed
d1: YANGLES ATOMS1=3,5 ATOMS2=1,2 MIN={BETA=0.1}
PRINT ARG=d1.min FILE=colvar
```

Notice that this command is a shortcut. You can thus learn more about how to use PLUMED by examining the expanded version of the input above.

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR ZANGLES
/*
Calculate the angle between an arbitrary vector and the positive z direction

__As you can see if you expand the inputs below, you can achieve what this shortcut action does by using [DISTANCE](DISTANCE.md) together with [CUSTOM](CUSTOM.md),
[BETWEEN](BETWEEN.md), [LESS_THAN](LESS_THAN.md), [SUM](SUM.md) and [MEAN](MEAN.md).  We strongly encourage you to use these actions instead as using them will provide
you with a clearer understanding of the equations you are using.__

The following input tells plumed to calculate the angles between the z-axis and the vector connecting atom 3 to atom 5 and between the z-axis
and the vector connecting atom 1 to atom 2.  The minimum of these two quantities is then output

```plumed
d1: ZANGLES ATOMS1=3,5 ATOMS2=1,2 MIN={BETA=0.1}
PRINT ARG=d1.min FILE=colvar
```

Notice that this command is a shortcut. You can thus learn more about how to use PLUMED by examining the expanded version of the input above.

*/
//+ENDPLUMEDOC


namespace PLMD {
namespace multicolvar {

class XAngle : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit XAngle(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(XAngle,"XANGLES")
PLUMED_REGISTER_ACTION(XAngle,"YANGLES")
PLUMED_REGISTER_ACTION(XAngle,"ZANGLES")

void XAngle::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("numbered","ATOMS","the pairs of atoms that you would like to calculate the angles for");
  keys.reset_style("ATOMS","atoms");
  MultiColvarShortcuts::shortcutKeywords( keys );
  keys.needsAction("DISTANCE");
  keys.needsAction("COMBINE");
  keys.needsAction("CUSTOM");
  keys.setDeprecated("DISTANCE");
}

XAngle::XAngle(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // Create distances
  std::string dline = getShortcutLabel() + "_dists: DISTANCE COMPONENTS";
  for(unsigned i=1;; ++i) {
    std::string atstring;
    parseNumbered("ATOMS",i,atstring);
    if( atstring.length()==0 ) {
      break;
    }
    std::string num;
    Tools::convert( i, num );
    dline += " ATOMS" + num + "=" + atstring;
  }
  readInputLine( dline );
  // Normalize the vectors
  readInputLine( getShortcutLabel() + "_norm2: COMBINE ARG=" + getShortcutLabel() + "_dists.x" + "," + getShortcutLabel() + "_dists.y," + getShortcutLabel() + "_dists.z POWERS=2,2,2 PERIODIC=NO");
  readInputLine( getShortcutLabel() + "_norm: CUSTOM ARG=" + getShortcutLabel() + "_norm2 FUNC=sqrt(x) PERIODIC=NO");
  readInputLine( getShortcutLabel() + "_norm_x: CUSTOM ARG=" + getShortcutLabel() + "_dists.x," + getShortcutLabel() + "_norm FUNC=x/y PERIODIC=NO");
  readInputLine( getShortcutLabel() + "_norm_y: CUSTOM ARG=" + getShortcutLabel() + "_dists.y," + getShortcutLabel() + "_norm FUNC=x/y PERIODIC=NO");
  readInputLine( getShortcutLabel() + "_norm_z: CUSTOM ARG=" + getShortcutLabel() + "_dists.z," + getShortcutLabel() + "_norm FUNC=x/y PERIODIC=NO");
  // Now compute the angles with matheval
  if( getName()=="XANGLES" ) {
    readInputLine( getShortcutLabel() + "_ang: CUSTOM FUNC=acos(x) PERIODIC=NO ARG=" + getShortcutLabel() + "_norm_x");
  }
  if( getName()=="YANGLES" ) {
    readInputLine( getShortcutLabel() + "_ang: CUSTOM FUNC=acos(x) PERIODIC=NO ARG=" + getShortcutLabel() + "_norm_y");
  }
  if( getName()=="ZANGLES" ) {
    readInputLine( getShortcutLabel() + "_ang: CUSTOM FUNC=acos(x) PERIODIC=NO ARG=" + getShortcutLabel() + "_norm_z");
  }
  // Add shortcuts to label
  MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel() + "_ang", "", this );
}

}
}
