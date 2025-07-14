/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#include "MultiColvarShortcuts.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

namespace PLMD {
namespace multicolvar {

//+PLUMEDOC COLVAR DIHCOR
/*
Measures the degree of similarity between dihedral angles.

This colvar calculates the following quantity.

$$
s = \frac{1}{2} \sum_i \left[ 1 + \cos( \phi_i - \psi_i ) \right]
$$

where the $\phi_i$ and $\psi_i$ values are the instantaneous values for the TORSION angles of interest.

You can see an example input for the DIHCOR action below

```plumed
dih: DIHCOR ...
  ATOMS1=1,2,3,4,5,6,7,8
  ATOMS2=5,6,7,8,9,10,11,12
...
PRINT ARG=dih FILE=colvar STRIDE=10
```

In the above input we are calculating the correlation between the torsion angle involving atoms 1, 2, 3 and 4 and the torsion angle
involving atoms 5, 6, 7 and 8.	This is then added to the correlation between the torsion angle involving atoms 5, 6, 7 and 8 and the
correlation angle involving atoms 9, 10, 11 and 12.

Writing out the atoms involved in all the torsion angles in this way can be rather tedious. Thankfully if you are working with protein you
can avoid this by using the [MOLINFO](MOLINFO.md) command.  PLUMED uses the pdb file that you provide to this command to learn
about the topology of the protein molecule.  This means that you can specify torsion angles using the following syntax:

```plumed
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO MOLTYPE=protein STRUCTURE=regtest/basic/rt32/helix.pdb
dih: DIHCOR ...
ATOMS1=@phi-3,@psi-3
ATOMS2=@psi-3,@phi-4
ATOMS3=@phi-4,@psi-4
...
PRINT ARG=dih FILE=colvar STRIDE=10
```

Here, `@phi-3` tells plumed that you would like to calculate the $\phi$ angle in the third residue of the protein.
Similarly `@psi-4` tells plumed that you want to calculate the $\psi$ angle of the fourth residue of the protein.

Notice, last of all, that if you want not to reassemble the atoms that have been broken by the periodic boundary conditions using a procedure
like that outlined in [WHOLEMOLECULES](WHOLEMOLECULES.md) you can add a NOPBC as shown below:

```plumed
dih: DIHCOR ...
  ATOMS1=1,2,3,4,5,6,7,8
  ATOMS2=5,6,7,8,9,10,11,12
  NOPBC
...
PRINT ARG=dih FILE=colvar STRIDE=10
```

*/
//+ENDPLUMEDOC

// We have a little helper class here to ensure that we actually do what is required by this action
class Dihcor : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  Dihcor(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(Dihcor,"DIHCOR")

void Dihcor::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.needsAction("DIHEDRAL_CORRELATION");
  keys.needsAction("SUM");
  keys.add("atoms","ATOMS","the set of 8 atoms that are being used each of the dihedral correlation values");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.setValueDescription("scalar","the sum of all the dihedral correlations");
}

Dihcor::Dihcor(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  readInputLine( getShortcutLabel() +"_data: DIHEDRAL_CORRELATION " + convertInputLineToString() );
  readInputLine( getShortcutLabel() + ": SUM ARG=" + getShortcutLabel() + "_data PERIODIC=NO");
}

}
}
