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
#include "Angle.h"
#include "ColvarShortcut.h"
#include "MultiColvarTemplate.h"
#include "core/ActionRegister.h"
#include "tools/Angle.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR ANGLE
/*
Calculate one or multiple angle/s.

The following input instructs PLUMED to calculate and print the angle between the vector
connecting atom 2 and atom 1 and the vector connecting atom 2 and atom 3.

```plumed
a1: ANGLE ATOMS=1,2,3
PRINT ARG=a1 FILE=colvar
```

In other words, the angle that is output by the input above is calculated as:

$$
\theta=\arccos\left(\frac{ r_{21}\cdot r_{23}}{
| r_{21}| | r_{23}|}\right)
$$

Here $r_{ij}$ is the vector connecting the $i$th and $j$th atoms, which by default is evaluated
in a way that takes periodic boundary conditions into account. If you wish to disregard the PBC you
can use the NOPBC flag as shown in the following input:

```plumed
a1: ANGLE ATOMS=1,2,3 NOPBC
PRINT ARG=a1 FILE=colvar
```

If the NOPBC flag is not included any sets of atoms that are broken by the periodic boundaries are made whole
using a procedure that is the same as that described in the documentation for [WHOLEMOLECULES](WHOLEMOLECULES.md).

We can also instruct PLUMED to calculate the angle between the vectors connecting
atoms 1 and atom 2 and atoms 3 and atom 4 by using the following input:

```plumed
a2: ANGLE ATOMS=1,2,3,4
PRINT ARG=a2 FILE=colvar

The angle in this input is calculated using:

$$
\theta=\arccos\left(\frac{ r_{21}\cdot r_{34}}{
| r_{21}| | r_{34}|}\right)
$$

Notice that angles defined in this way are non-periodic variables - their values must lie in between 0 and $\pi$.

You can specify multiple sets of three or four atoms to calculate vectors of angles as illustrated in the following
input which instructs PLUMED to calculate and output three angles:

```plumed
a3: ANGLE ATOMS1=1,2,3 ATOMS2=4,5,6 ATOMS3=7,8,9
PRINT ARG=a3 FILE=colvar
```

It is common to assume when using this feature that all the angles being computed are indistinguishable
so it makes sense to perform the same series of operations on every element of the output vector. The input
file below is more approrpriate if the angles are not indistinguishable:

```plumed
a4: ANGLE ATOMS=1,2,3
a5: ANGLE ATOMS=4,5,6
a6: ANGLE ATOMS=7,8,9
PRINT ARG=a4,a5,a6 FILE=colvar
```

*/
//+ENDPLUMEDOC

typedef ColvarShortcut<Angle> AngleShortcut;
PLUMED_REGISTER_ACTION(AngleShortcut,"ANGLE")
PLUMED_REGISTER_ACTION(Angle,"ANGLE_SCALAR")
typedef MultiColvarTemplate<Angle> AngleMulti;
PLUMED_REGISTER_ACTION(AngleMulti,"ANGLE_VECTOR")

} //namespace colvar
} //namespace PLMD

