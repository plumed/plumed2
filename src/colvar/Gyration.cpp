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
#include "Gyration.h"
#include "ColvarShortcut.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR GYRATION
/*
Calculate the radius of gyration, or other properties related to it.

The radius of gyration is calculated using:

$$
s_{\rm Gyr}=\Big ( \frac{\sum_i^{n} w_i \vert {r}_i -{r}_{\rm COM} \vert ^2 }{\sum_i^{n} w_i} \Big)^{1/2}
$$

with the position of the center ${r}_{\rm COM}$ given by:

$$
{r}_{\rm COM}=\frac{\sum_i^{n} {r}_i\ w_i }{\sum_i^{n} w_i}
$$

In these expressions $r_i$ indicates the position of atom $i$ and $w_i$ is the weight for atom $i$.  The following input
shows how you can calculate the expressions for a set of atoms by using PLUMED:

```plumed
w: CONSTANT VALUES=1,2,2,3,4
g: GYRATION ATOMS=1-5 TYPE=RADIUS WEIGHTS=w
PRINT ARG=g FILE=colvar
```

In the above input the weights in the expressions above are set equal to the elemnts of a constant vector, `w`.  A more common
approach is to use the masses of the atoms, which you can do using any one of the three inputs below:

```plumed
g: GYRATION ATOMS=1-5 TYPE=RADIUS MASS_WEIGHTED
# This input is equivalent
g2: GYRATION ATOMS=1-5 TYPE=RADIUS WEIGHTS=@Masses
# As is this one
g3: GYRATION ATOMS=1-5 TYPE=RADIUS MASS
# The following print statement thus outputs the same number three times
PRINT ARG=g,g2,g3 FILE=colvar
```

Alternatively, you can set all the input weights equal to one as is done in the input below;

```plumed
g: GYRATION ATOMS=1-5 TYPE=RADIUS
PRINT ARG=g FILE=colvar
```

Notice that in the first input above GYRATION is a shortcut and that in the second two inputs it is not.  This is because there
is a faster (non-shortcutted) version that can be used for these two more-commonplace inputs.  The shortcut input is still useful
as it allows one to do less comonplace analyses such as this one where a clustering is performed and the radius of gyration for the
largest cluster is determined.

```plumed
c1: CONTACT_MATRIX GROUP=1-100 SWITCH={RATIONAL R_0=0.1 D_MAX=0.5}
dfs: DFSCLUSTERING ARG=c1
w: CLUSTER_WEIGHTS CLUSTERS=dfs CLUSTER=1
g: GYRATION ATOMS=1-100 WEIGHTS=w TYPE=RADIUS
```

In calculating the gyration radius you first calculate the [GYRATION_TENSOR](GYRATION_TENSOR.md). Instad of computing the radius from this
tensor you can instead compute the trace of the tensor by using the following input:

```plumed
w: CONSTANT VALUES=1,2,2,3,4
g: GYRATION ATOMS=1-5 TYPE=TRACE WEIGHTS=w
PRINT ARG=g FILE=colvar
```

Notice, that when you compute the gyration tensor in the above calculation it is unormalized.  If you want to calculate any other
quantity from the unormalized gyration tensor you can add the UNORMALIZED flag as shown below, which computes the unormalized radius
of gyration.

```plumed
w: CONSTANT VALUES=1,2,2,3,4
g: GYRATION ATOMS=1-5 UNORMALIZED TYPE=RADIUS WEIGHTS=w
PRINT ARG=g FILE=colvar
```

You can also calculate the largest, middle and smallest principal moments of the normalized tensor:

```plumed
w: CONSTANT VALUES=1,2,2,3,4
# This computes the largest principal moment
l: GYRATION ATOMS=1-5 TYPE=GTPC_1 WEIGHTS=w
# This computes the middle principal moment
m: GYRATION ATOMS=1-5 TYPE=GTPC_2 WEIGHTS=w
# This computes the smallest principal moment
s: GYRATION ATOMS=1-5 TYPE=GTPC_3 WEIGHTS=w
PRINT ARG=l,m,s FILE=colvar
```

From these principal moments you can calculate the Asphericiry, Acylindricity, or the Relative Shape Anisotropy by using the following input:

```plumed
w: CONSTANT VALUES=1,2,2,3,4
# This computes the Asphericiry
l: GYRATION ATOMS=1-5 TYPE=ASPHERICITY WEIGHTS=w
# This computes the Acylindricity
m: GYRATION ATOMS=1-5 TYPE=ACYLINDRICITY WEIGHTS=w
# This computes the Relative Shape Anisotropy
s: GYRATION ATOMS=1-5 TYPE=KAPPA2 WEIGHTS=w
PRINT ARG=l,m,s FILE=colvar
```

Lastly, you can calculate the largest, middle and smallest principal radii of gyration by using the following input:

```plumed
w: CONSTANT VALUES=1,2,2,3,4
# This computes the largest principal radius of gyration
l: GYRATION ATOMS=1-5 TYPE=RGYR_1 WEIGHTS=w
# This computes the middle principal radius of gyration
m: GYRATION ATOMS=1-5 TYPE=RGYR_2 WEIGHTS=w
# This computes the smallest principal radius of gyration
s: GYRATION ATOMS=1-5 TYPE=RGYR_3 WEIGHTS=w
PRINT ARG=l,m,s FILE=colvar
```

By expanding the shortcuts in the inputs above or by reading the papers cited below you can see how these quantities are computed from the [GYRATION_TENSOR](GYRATION_TENSOR.md).
Notice, however, that as with the radius if you use the masses or a vector of ones as weights a faster implemention of this colvar that
calculates these quantities directly will be used.

## A note on periodic boundary conditions

Calculating the radius of gyration is normally used to determine the shape of a molecule so all the specified atoms
would normally be part of the same molecule.  When computing this CV it is important to ensure that the periodic boundaries
are calculated correctly.  There are two ways that you can manage periodic boundary conditions when using this action.  The
first and simplest is to reconstruct the molecule similarly to the way that [WHOLEMOLECULES](WHOLEMOLECULES.md) operates.
This reconstruction of molecules has been done automatically since PLUMED 2.2.  If for some reason you want to turn it off
you can use the NOPBC flag as shown below:

```plumed
g: GYRATION ATOMS=1-5 TYPE=RADIUS MASS_WEIGHTED NOPBC
PRINT ARG=g FILE=colvar
```

An alternative approach to handling PBC is to use the PHASES keyword.  This keyword instructs PLUMED to use the PHASES option
when computing the position of the center using the [CENTER](CENTER.md) command.  Distances of atoms from this center are then
computed using PBC as usual. The example shown below shows you how to use this option

```plumed
g: GYRATION ATOMS=1-5 TYPE=RADIUS MASS_WEIGHTED PHASES
PRINT ARG=g FILE=colvar
```

## Working with multiple molecules

If you are working with a system that contains multiple identical molecules and you wish to calculate gyration radii for all the molecules you can do this with multiple GYRATION actions like this

```plumed
g1: GYRATION ATOMS=1-6 TYPE=RADIUS
g2: GYRATION ATOMS=7-12 TYPE=RADIUS
g3: GYRATION ATOMS=13-18 TYPE=RADIUS
PRINT ARG=g1,g2,g3 FILE=colvar
```

Alternatively, you can use a single line command such as the one below:

```plumed
g: GYRATION ATOMS1=1-6 ATOMS2=7-12 ATOMS3=13-18 TYPE=RADIUS
PRINT ARG=g FILE=colvar
```

The value outputted by above input is a three dimensional vector that contains the radii of gyration computed for atoms 1-6, 7-12 and 13-18.  If you want to calculate a vector that contains multiple instances of one of the other
quantities that have been defined on this page you change the input to the TYPE keyword.  Similarly, if you want to use the masses of the atoms for the weights at the top of this page you can use the MASS_WEIGHTED flag.

The single line option that computes the vector is preferable to the multi-line option that computes multiple scalars because the calculation is parallelized in this second case. Furthermore, there will be fewer virtual
function calls if you use the second option.  This second calculation will thus likely run faster than the first.  However, you should note that if you are using the WEIGHTS or PHASES keywords you can only compute scalars.
There is currently no option to use WEIGHTS and PHASES and multiple ATOMS keywords.

*/
//+ENDPLUMEDOC

typedef Gyration<double> GyrationD;
PLUMED_REGISTER_ACTION(GyrationD,"GYRATION_SCALAR")
typedef MultiColvarTemplate<GyrationD> GyrationMulti;
PLUMED_REGISTER_ACTION(GyrationMulti,"GYRATION_VECTOR")

}
}
