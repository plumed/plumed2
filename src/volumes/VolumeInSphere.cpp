/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#include "VolumeInSphere.h"
#include "core/ActionRegister.h"
#include "VolumeShortcut.h"

//+PLUMEDOC VOLUMES INSPHERE
/*
This quantity can be used to calculate functions of the distribution of collective variables for the atoms that lie in a particular, user-specified part of of the cell.

This action can be used to calculate whether each of the atoms are within a particular part of the simulation box or not as illustrated by the following example:

```plumed
f: FIXEDATOM AT=0,0,0
a: INSPHERE ATOMS=1-100 CENTER=f RADIUS={GAUSSIAN R_0=1.5}
PRINT ARG=a FILE=colvar
```

The 100 elements of the vector `a` that is returned from the INSPHERE action in the above input are calculated using:

$$
w(x_i,y_i,z_i) = \sigma\left( \sqrt{ x_i^2 + y_i^2 + z_i^2}  \right)
$$

In this expression $x_i$, $y_i$ and $z_i$ are the components of the vector that connects the $i$th atom that was specified using the `ATOMS` keyword to the atom that was specified using the `CENTER` keyword and
$\sigma$ is one of the switching functions that is described that in the documentation for the action [LESS_THAN](LESS_THAN.md).  In other words,
$w(x_i,y_i,z_i)$ is 1 if atom $i$ is within a sphere with a radial extent that is determined by the parameters of the specified switching function
and zero otherwise.

If you want to caluclate and print the number of atoms in the sphere of interest you can use an input like the one shown below:

```plumed
f: FIXEDATOM AT=0,0,0
a: INSPHERE ATOMS=1-100 CENTER=f RADIUS={GAUSSIAN R_0=1.5}
s: SUM ARG=a PERIODIC=NO
PRINT ARG=s FILE=colvar
```

If by constrast you want to calculate and print the number of atoms that are not in the sphere of interest you OUTSIDE flag as shown below

```plumed
f: FIXEDATOM AT=0,0,0
a: INSPHERE ...
  ATOMS=1-100 CENTER=f
  RADIUS={GAUSSIAN R_0=1.5}
  OUTSIDE
...
s: SUM ARG=a PERIODIC=NO
PRINT ARG=s FILE=colvar
```

!!! note ""

    You can also calculate the average values of symmetry functions in the sphere of interest by using inputs similar to those described the documentation for the [AROUND](AROUND.md)
    action. In other words, you can swap out AROUND actions for an INSPHERE actions.  Also as with [AROUND](AROUND.md), earlier versions of PLUMED used a different syntax for doing these types of calculations, which can
    still be used with this new version of the command.  We strongly recommend using the newer syntax but if you are interested in the
    old syntax you can find more information in the old syntax section of the documentation for [AROUND](AROUND.md).

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace volumes {

typedef ActionVolume<VolumeInSphere> Vols;
PLUMED_REGISTER_ACTION(Vols,"INSPHERE_CALC")
char glob_sphere[] = "INSPHERE";
typedef VolumeShortcut<glob_sphere> VolumeInSphereShortcut;
PLUMED_REGISTER_ACTION(VolumeInSphereShortcut,"INSPHERE")
} // namespace PLMD
} // namespace volumes
