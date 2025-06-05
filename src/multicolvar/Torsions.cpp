/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2020 The plumed team
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

//+PLUMEDOC MCOLVAR TORSIONS
/*
Calculate whether or not a set of torsional angles are within a particular range.

__This shortcut action allows you to calculate function of a distribution of torsions and reproduces the syntax in older PLUMED versions.
If you look at the example inputs below you can
see how the new syntax operates. We would strongly encourage you to use the newer syntax as it offers greater flexibility.__

The following provides an example of the input for the TORSIONS command

```plumed
ab: TORSIONS ...
  ATOMS1=168,170,172,188
  ATOMS2=170,172,188,190
  ATOMS3=188,190,192,230
  BETWEEN={GAUSSIAN LOWER=0 UPPER=pi SMEAR=0.1}
...
PRINT ARG=ab.* FILE=colvar STRIDE=10
```

The input above calculates how many of torsion angles for the three groups of atoms that have been specified are between 0 and $\pi$.

Writing out the atoms involved in all the torsion angles in this way can be rather tedious. Thankfully if you are working with protein you
can avoid this by using the [MOLINFO](MOLINFO.md) command.  PLUMED uses the pdb file that you provide to this command to learn
about the topology of the protein molecule.  This means that you can specify torsion angles using the following syntax:

```plumed
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO MOLTYPE=protein STRUCTURE=regtest/basic/rt32/helix.pdb
TORSIONS ...
ATOMS1=@phi-3
ATOMS2=@psi-3
ATOMS3=@phi-4
BETWEEN={GAUSSIAN LOWER=0 UPPER=pi SMEAR=0.1}
LABEL=ab
... TORSIONS
PRINT ARG=ab.* FILE=colvar STRIDE=10
```

Here, `@phi-3` tells plumed that you would like to calculate the $\phi$ angle in the third residue of the protein.
Similarly `@psi-4` tells plumed that you want to calculate the $\psi$ angle of the fourth residue of the protein.


*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

class Torsions : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit Torsions(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(Torsions,"TORSIONS")

void Torsions::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  MultiColvarShortcuts::shortcutKeywords( keys );
  keys.needsAction("TORSION");
  keys.add("atoms","ATOMS","the four atoms involved in the torsional angle");
  keys.setValueDescription("vector","the TORSION for each set of three atoms that were specified");
  keys.setDeprecated("TORSION");
}

Torsions::Torsions(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  log.printf("Action TORSION\n");
  log.printf("  with label %s \n", getShortcutLabel().c_str() );
  std::map<std::string,std::string> keymap;
  MultiColvarShortcuts::readShortcutKeywords( keymap, this );
  readInputLine( getShortcutLabel() + ": TORSION " + convertInputLineToString() );
  // Add shortcuts to label
  MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel(), "", keymap, this );
}

}
}
