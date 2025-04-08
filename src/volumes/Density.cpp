/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2023 The plumed team
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

//+PLUMEDOC MCOLVAR DENSITY
/*
Depreciated command that is bascially equivalant to GROUP.

Here is an example but Plase don't use this anymore.  Use [GROUP](GROUP.md) instead.

```plumed
g1: DENSITY SPECIES=1-100
DUMPATOMS ATOMS=g1 FILE=group.xyz
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace volumes {

class Density : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit Density(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(Density,"DENSITY")

void Density::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.setDeprecated("GROUP");
  keys.add("compulsory","SPECIES","the atoms in the group");
  keys.setValueDescription("atoms","indices for the specified group of atoms");
  keys.needsAction("ONES");
  keys.needsAction("GROUP");
}

Density::Density(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  std::string atoms;
  parse("SPECIES",atoms);
  warning("This action has been depracated.  Look at the log to see how the same result is achieved with the new syntax");
  readInputLine( getShortcutLabel() + ": GROUP ATOMS=" + atoms);
}

}
}
