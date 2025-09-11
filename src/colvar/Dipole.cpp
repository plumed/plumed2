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
#include "Dipole.h"
#include "ColvarShortcut.h"
#include "MultiColvarTemplate.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR DIPOLE
/*
Calculate the dipole moment for a group of atoms.

The following input tells plumed to calculate the dipole for the group of atoms containing
the atoms from 1-10 and print it every 5 steps

```plumed
d: DIPOLE GROUP=1-10
PRINT FILE=output STRIDE=5 ARG=d
```

The output value from this input is a scalar that tells you the total magnitude of the dipole vector.  If you would
like to access the dipole vector directly you can use the command:

```plumed
d: DIPOLE GROUP=1-10 COMPONENTS
PRINT ARG=d.* FILE=output
```

This command will output three values d.x, d.y and d.z, which are the x, y and z components of the dipole respectively.


You can calculate three instinguishable dipoles using a single DIPOLE command by using an input like the one below:

```plumed
d: DIPOLE GROUP1=1-10 GROUP2=11-20 GROUP3=21-30
PRINT ARG=d FILE=output
```

The output, d, here is a three dimensional vector.  The first element of this vector is the magnitude of the dipole for
atoms 1-10, the second is the magnitude of the dipole for atoms 11-20 and the third is the magnitude of the dipole for
atoms 21-30.  You can also obtain vector components for the three dipoles above by using the following input:

```plumed
d: DIPOLE COMPONENTS GROUP1=1-10 GROUP2=11-20 GROUP3=21-30
PRINT ARG=d.x,d.y,d.z FILE=output
```

The output from the DIPOLE command now consists of three three dimensional vectors called d.x, d.y and d.z that contain the
x, y and z components of the three dipoles respectively.

A final important thing to note is that in all the commands above the default is to use a procedure akin to that used in [WHOLEMOLECULES](WHOLEMOLECULES.md) to ensure
that the sets of atoms that are specified to each GROUP keyword are not broken by the periodic
boundary conditions.  If you would like to turn this off for any reason you add the NOPBC in your input file as shown
below:

```plumed
d: DIPOLE GROUP=1-10 NOPBC
PRINT FILE=output STRIDE=5 ARG=d
```

!!! caution "behaviour for non charge neutral groups"

    If the total charge Q of any of the specified groups is non zero, then a charge Q/N will be subtracted from every atom,
    where N is the number of atoms. This implies that the dipole (which for a charged system depends
    on the position) is computed on the geometric center of the group.

*/
//+ENDPLUMEDOC

typedef Dipole<double> DipoleD;
typedef ColvarShortcut<DipoleD> DipoleShortcut;
PLUMED_REGISTER_ACTION(DipoleShortcut,"DIPOLE")
PLUMED_REGISTER_ACTION(DipoleD,"DIPOLE_SCALAR")
typedef MultiColvarTemplate<DipoleD> DipoleMulti;
PLUMED_REGISTER_ACTION(DipoleMulti,"DIPOLE_VECTOR")

}
}
