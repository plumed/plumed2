/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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

#include "Distance.h"
#include "ColvarShortcut.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR DISTANCE
/*
Calculate the distance/s between pairs of atoms.

The following example illustrates how this action can be used to calculate and print the distance between atom 1
and atom 2.

```plumed
d: DISTANCE ATOMS=1,2
PRINT ARG=d FILE=colvar
```

By default the distance is computed in a way that takes periodic
boundary conditions in account. This behavior can be changed by using the NOPBC flag.
Furthermore, if you wish to calculate the vector connecting a pair of atoms you can use the
`COMPONENTS` flag as shown below:

```plumed
d: DISTANCE ATOMS=1,2 COMPONENTS
PRINT ARG=d.x,d.y,d.z FILE=colvar
```

Alternatively, you can calculate the components projected on the lattice vector by using the `SCALED_COMPONENTS`
flag as shown below;

```plumed
d: DISTANCE ATOMS=1,2 SCALED_COMPONENTS
PRINT ARG=d.a,d.b,d.c FILE=colvar
```

The advantage of using `SCALED_COMPONENTS` over `COMPONENTS` is that the a, b and c variables
that are calculated when `SCALED_COMPONENTS` is employed have the proper periodicity. This feature is useful
if you wish to study the motion of a molecule across a membrane.

You can also use this command to calculate multiple indistinguishable distances or vectors with a single
line of PLUMED input.  For example, the following input calculates and outputs the distances between four
pairs of atoms:

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8
PRINT ARG=d FILE=colvar
```

By a similar token, the following input outputs three four dimensional vectors that contain the x, y and z
components of the vectors connecting the four atoms:

```plumed
d: DISTANCE COMPONENTS ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8
PRINT ARG=d.x,d.y,d.z FILE=colvar
```

You can also replace COMPONENTS with SCALED_COMPONENTS in the above input and obtain the projects of these vectors
on the lattice vectors.

## Managing periodic boundary conditions

When using the DISTANCE command to calculate the end-to-end distance for a large polymer you need to ensure that you
are managing PBCs correctly.  This problems that can occur with these calculations are explained at length in the
early parts of the document that is referenced in the bibliography. Notice, however, that the input
provides an example of an input that could be used to compute the end-to-end distance for a polymer
of 100 atoms and keeps it at a value around 5.

```plumed
WHOLEMOLECULES ENTITY0=1-100
e2e: DISTANCE ATOMS=1,100 NOPBC
RESTRAINT ARG=e2e KAPPA=1 AT=5
```

Notice that NOPBC is used here so as to ensure that the distance is calculated correctely even if the end-to-end distance is larger than half the simulation
Also notice that, since many MD codes break molecules across cell boundary, it might be necessary to
use the [WHOLEMOLECULES](WHOLEMOLECULES.md) keyword (also notice that it should be _before_ distance). The list of atoms provided to [WHOLEMOLECULES](WHOLEMOLECULES.md)
here contains all the atoms between 1 and 100. Strictly speaking, this
is not necessary. If you know for sure that atoms with difference in
the index say equal to 10 are _not_ going to be farther than half cell
you can e.g. use

```plumed
WHOLEMOLECULES ENTITY0=1,10,20,30,40,50,60,70,80,90,100
e2e: DISTANCE ATOMS=1,100 NOPBC
RESTRAINT ARG=e2e KAPPA=1 AT=5
```

Just be sure that the ordered list provide to [WHOLEMOLECULES](WHOLEMOLECULES.md) has the following
properties:

- Consecutive atoms should be closer than half-cell throughout the entire simulation.
- Atoms required later for the distance (e.g. 1 and 100) should be included in the list

The following example shows how to take periodicity into account when computing the z-component of a distance

```plumed
# this is a center of mass of a large group
c: COM ATOMS=1-100
# this is the distance between atom 101 and the group
d: DISTANCE ATOMS=c,101 COMPONENTS
# this makes a new variable, dd, equal to d and periodic, with domain -10,10
# this is the right choise if e.g. the cell is orthorombic and its size in
# z direction is 20.
dz: COMBINE ARG=d.z PERIODIC=-10,10
# metadynamics on dd
METAD ARG=dz SIGMA=0.1 HEIGHT=0.1 PACE=200
```

You can use the same input even if the DISTANCE command is calculating the vectors connecting multiple pairs of atoms.
However, using SCALED_COMPONENTS ensures this problem does not arise because these variables are always periodic
with domain (-0.5,+0.5).

*/
//+ENDPLUMEDOC

typedef ColvarShortcut<Distance> DistanceShortcut;
PLUMED_REGISTER_ACTION(DistanceShortcut,"DISTANCE")
PLUMED_REGISTER_ACTION(Distance,"DISTANCE_SCALAR")
typedef MultiColvarTemplate<Distance> DistanceMulti;
PLUMED_REGISTER_ACTION(DistanceMulti,"DISTANCE_VECTOR")

} //namespace colvar
} //namespace PLMD
