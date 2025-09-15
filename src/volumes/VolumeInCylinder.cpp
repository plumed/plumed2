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
#include "VolumeInCylinder.h"
#include "core/ActionRegister.h"
#include "tools/Pbc.h"
#include "tools/SwitchingFunction.h"
#include "ActionVolume.h"
#include "tools/HistogramBead.h"
#include "VolumeShortcut.h"

//+PLUMEDOC VOLUMES INCYLINDER
/*
This quantity can be used to calculate functions of the distribution of collective variables for the atoms that lie in a particular, user-specified part of of the cell.

This action can be used to calculate whether each of the atoms are within a particular part of the simulation box or not as illustrated by the following example:

```plumed
f: FIXEDATOM AT=0,0,0
a: INCYLINDER ...
   ATOMS=1-100 CENTER=f DIRECTION=Z
   RADIUS={TANH R_0=1.5} SIGMA=0.1
   LOWER=-1.0 UPPER=1.0
   KERNEL=gaussian
...
PRINT ARG=a FILE=colvar
```

The 100 elements of the vector `a` that is returned from the INCYLINDER action in the above input are calculated using:

$$
w(x_i,y_i,z_i) = s\left( r_{xy} \right) \int_{zl}^{zu} \textrm{d}z K\left( \frac{z - z_i}{\sigma} \right)
$$

where

$$
r_{xy} = \sqrt{ x_i^2 + y_i^2 }
$$

In these expressions $K$ is one of the kernel functions described in the documentation for the function [BETWEEN](BETWEEN.md), $\sigma$ is a bandwidth parameter and the limits
for the integrals are the values specified using the keywords `LOWER` and `UPPER`. $x_i$, $y_i$ and $z_i$ are then the components
of the vector that connects the $i$th atom that was specified using the `ATOMS` keyword to the atom that was specified using the `CENTER` keyword.  In other words,
$w(x_i,y_i,z_i)$ is 1 if the atom is within a cylinder that is centered on the $z$ axis that has an extent along the $z$ direction around the position of atom `f` that
is determined by the values specified by `LOWER` and `UPPER` keywords.  The radial extent of this cylinder is determined by the parameters of the switching function that is
specified using the `RADIUS` keyword.

If you want to caluclate and print the number of atoms in the cylinder of interest you can use an input like the one shown below:

```plumed
f: FIXEDATOM AT=0,0,0
a: INCYLINDER ...
  ATOMS=1-100 CENTER=f DIRECTION=X
  RADIUS={TANH R_0=1.5} SIGMA=0.1
  LOWER=-1.0 UPPER=1.0
...
s: SUM ARG=a PERIODIC=NO
PRINT ARG=s FILE=colvar
```

Alternatively, if you want to calculate an print the number of atoms that are not in the cylinder of interest you can use the OUTSIDE flag as shown below:

```plumed
f: FIXEDATOM AT=0,0,0
a: INCYLINDER ...
  ATOMS=1-100 CENTER=f DIRECTION=X
  RADIUS={TANH R_0=1.5} SIGMA=0.1
  LOWER=-1.0 UPPER=1.0
  OUTSIDE
...
s: SUM ARG=a PERIODIC=NO
PRINT ARG=s FILE=colvar
```

Notice that now the cylinder is centered on the `x` axis rather than the `z` axis as we have changed the input for the `DIRECTION` keyword.

!!! note ""

    You can also calculate the average values of symmetry functions in the cylinder of interest by using inputs similar to those described the documentation for the [AROUND](AROUND.md)
    action. In other words, you can swap out AROUND actions for an INCLYLINDER actions. Also as with [AROUND](AROUND.md), earlier versions of PLUMED used a different syntax for doing these
    types of calculations, which can still be used with this new version of the command.  We strongly recommend using the newer syntax but if you are interested in the
    old syntax you can find more information in the old syntax section of the documentation for [AROUND](AROUND.md).


*/
//+ENDPLUMEDOC

namespace PLMD {
namespace volumes {


typedef ActionVolume<VolumeInCylinder> Volc;
PLUMED_REGISTER_ACTION(Volc,"INCYLINDER_CALC")
char glob_cylinder[] = "INCYLINDER";
typedef VolumeShortcut<glob_cylinder> VolumeInCylinderShortcut;
PLUMED_REGISTER_ACTION(VolumeInCylinderShortcut,"INCYLINDER")

}
}
