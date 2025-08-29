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
#include "Plane.h"
#include "ColvarShortcut.h"
#include "core/ActionRegister.h"
#include "MultiColvarTemplate.h"
#include "tools/Pbc.h"

//+PLUMEDOC COLVAR PLANE
/*
Calculate the plane perpendicular to two vectors in order to represent the orientation of a planar molecule.

To calculate the orientation of the plane connecting atoms 1, 2 and 3 you use an input like this:

```plumed
p: PLANE ATOMS=1,2,3
PRINT ARG=p.x,p.y,p.z FILE=colvar
```

The three components, p.x, p.y and p.z, output by the PLANE action here are the x, y and z components of the normal
vector to the plane that is obtained by taking the cross product between the vector connecting atoms 1 and 2 and
the vector connecting atoms 2 and 3.  Notice that the default here is to evaluate these two vectors in a way that takes
any periodic boundary conditions (PBC) into account. If you wish to disregard the PBC you can use the NOPBC flag as shown in the following input:

```plumed
p: PLANE ATOMS=1,2,3 NOPBC
PRINT ARG=p.x,p.y,p.z FILE=colvar
```

To calculate the cross product of the vector connecting atoms 1 and 2 and the vector connecting atoms 3 and 4 you use
an input like this:

```plumed
p: PLANE ATOMS=1,2,3,4
PRINT ARG=p.x,p.y,p.z FILE=colvar
```

If you have multiple molecules and wish to determine the orientations of the planes containing all them with one line of PLUMED
input you can use an input like this:

```plumed
p: PLANE ATOMS1=1,2,3 ATOMS2=4,5,6 ATOMS3=7,8,9 ATOMS4=10,11,12
PRINT ARG=p.x,p.y,p.z FILE=colvar
```

The output from this command consists of 3 vectors with 4 components. These vectors, p.x, p.y and p.z, contain the x, y and z components
of the normals to the planes of the molecules.  Commands similar to this are useful for variables that can be used to monitor
nucleation of molecular crystals such as [SMAC](SMAC.md).

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace colvar {

typedef ColvarShortcut<Plane> PlaneShortcut;
PLUMED_REGISTER_ACTION(PlaneShortcut,"PLANE")
PLUMED_REGISTER_ACTION(Plane,"PLANE_SCALAR")
typedef MultiColvarTemplate<Plane> PlaneMulti;
PLUMED_REGISTER_ACTION(PlaneMulti,"PLANE_VECTOR")


}
}

