/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2023 The plumed team
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
#include "Position.h"
#include "ColvarShortcut.h"
#include "MultiColvarTemplate.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR POSITION
/*
Calculate the components of the position of an atom or atoms.

To print the position of atom one to a file you can use an input like this:

```plumed
p: POSITION ATOM=1
PRINT ARG=p.x,p.y,p.z FILE=colvar
```

To print the position of four atoms you can use an input like this:

```plumed
p: POSITION ATOMS=1-4
PRINT ARG=p.x,p.y,p.z FILE=colvar
```

The three output values, p.x, p.y and p.z, here are all four dimensional vectors.  Furthermore, if you wish
to use a procedure akin to that described in the documentation for [WHOLEMOLECULES](WHOLEMOLECULES.md) to ensure
that molecules are reassembled if they are broken by periodic boundary conditions you can use the `WHOLEMOLECULES`
flag as shown below:

```plumed
p: POSITION ATOMS=1-4 WHOLEMOLECULES
PRINT ARG=p.x,p.y,p.z FILE=colvar
```

This is used in the [CENTER](CENTER.md) when we compute centers with arbitrary weights for shortcuts.

## Periodic boundary conditions

!!! warning ""

    Notice that single components will not have the proper periodicity!

If you need the values to be consistent through PBC you can use SCALED_COMPONENTS flag as shown below:

```plumed
p: POSITION ATOM=1 SCALED_COMPONENTS
PRINT ARG=p.a,p.b,p.c FILE=colvar
```

This flag ensures that values are output that are, by construction, in the -0.5,0.5 domain. The output is
similar to the equivalent flag for [DISTANCE](DISTANCE.md).

If it also important to note that by default the positions are output by calculating minimal image distance between
the instantaneous position of the atoms and the position of the origin, $(0.0,0.0,0.0)$.
This behaviour can be changed by using NOPBC as shown below, which will output the position of the atom directly.

```plumed
p: POSITION ATOM=1 NOPBC
PRINT ARG=p.x,p.y,p.z FILE=colvar
```

!!! warning ""

    This variable should be used with extreme care since it allows you to easily get in troubles.
    It can be only be used if the Hamiltonian is not invariant for translation (i.e. there are other absolute positions which are biased, e.g. by position restraints)
    and cell size and shapes are fixed through the simulation.

If you are not in this situation and still want to use the absolute position of an atom you should first fix the reference frame
by using [FIT_TO_TEMPLATE](FIT_TO_TEMPLATE.md) as shown in the example below

```plumed
#SETTINGS INPUTFILES=regtest/basic/rt63/align.pdb
# align to a template
FIT_TO_TEMPLATE REFERENCE=regtest/basic/rt63/align.pdb
p: POSITION ATOM=3
PRINT ARG=p.x,p.y,p.z
```

*/
//+ENDPLUMEDOC

typedef Position<double> PositionD;
typedef ColvarShortcut<PositionD> PositionShortcut;
PLUMED_REGISTER_ACTION(PositionShortcut,"POSITION")
PLUMED_REGISTER_ACTION(PositionD,"POSITION_SCALAR")
typedef MultiColvarTemplate<PositionD> PositionMulti;
PLUMED_REGISTER_ACTION(PositionMulti,"POSITION_VECTOR")


}
}



