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
#include "VolumeAround.h"
#include "core/ActionRegister.h"
#include "tools/Pbc.h"
#include "tools/HistogramBead.h"
#include "VolumeShortcut.h"

//+PLUMEDOC VOLUMES AROUND
/*
This quantity can be used to calculate functions of the distribution of collective variables for the atoms that lie in a particular, user-specified part of of the cell.

This action can be used to calculate whether each of the atoms are within a particular part of the simulation box or not as illustrated by the following example:

```plumed
f: FIXEDATOM AT=0,0,0
a: AROUND ...
  ATOMS=1-100 ORIGIN=f
  SIGMA=0.2 KERNEL=gaussian
  XLOWER=-1.0 XUPPER=1.0
  YLOWER=-1.0 YUPPER=1.0
  ZLOWER=-1.0 ZUPPER=1.0
...
PRINT ARG=a FILE=colvar
```

The 100 elements of the vector `a` that is returned from the AROUND action in the above input are calculated using:

$$
w(x_i,y_i,z_i) = \int_{xl}^{xu} \int_{yl}^{yu} \int_{zl}^{zu} \textrm{d}x\textrm{d}y\textrm{d}z K\left( \frac{x - x_i}{\sigma} \right)K\left( \frac{y - y_i}{\sigma} \right)K\left( \frac{z - z_i}{\sigma} \right)
$$

where $K$ is one of the kernel functions described in the documentation for the function [BETWEEN](BETWEEN.md), $\sigma$ is a bandwidth parameter and the limits
for the integrals are the values specified using the keywords XLOWER, XUPPER, YLOWER, YUPPER, YUPPER, ZLOWER and ZUPPER.  $x_i$, $y_i$ and $z_i$ are then the components
of the vector that connects the $i$th atom that was specified using the ATOMS keyword to the atom that was specified using the ORIGIN keyword.  In other words,
$w(x_i,y_i,z_i)$ is 1 if the atom is within a rectangular box that is centered on the atom that is specified as the origin and zero otherwise.

If instead of calculating if the atoms are inside this box you want to calculate if they are outside this box you can use the following input:

```plumed
f: FIXEDATOM AT=0,0,0
a: AROUND ...
  ATOMS=1-100 ORIGIN=f
  SIGMA=0.2 KERNEL=gaussian
  XLOWER=-1.0 XUPPER=1.0
  YLOWER=-1.0 YUPPER=1.0
  ZLOWER=-1.0 ZUPPER=1.0
  OUTSIDE
...
PRINT ARG=a FILE=colvar
```

The 100 elements of the vector `a` that is returned from the AROUND action in the above input are calculated using:

$$
v(x_i,y_i,z_i) = 1 - w(x_i,y_i,z_i)
$$

## Calculating the number of atoms in a particular part of the box

Lets suppose that you want to calculate how many atoms are in have an $x$ coordinate that is between -1.0 and 1.0. You can do this using the following PLUMED input:

```plumed
f: FIXEDATOM AT=0,0,0
a: AROUND ATOMS=1-100 ORIGIN=f SIGMA=0.2 XLOWER=-1.0 XUPPER=1.0
s: SUM ARG=a PERIODIC=NO
PRINT ARG=s FILE=colvar
```

In this example the components of `a` are calculated as:

$$
w(x_i,y_i,z_i) = \int_{xl}^{xu} \textrm{d}x K\left( \frac{x - x_i}{\sigma} \right)
$$

as the YLOWER, YUPPER, YUPPER, ZLOWER and ZUPPER flags have not been included.  The [SUM](SUM.md) command then adds together all the elements of the vector `a` to calculate the total number of atoms in the region
of the box that is of interest.

## Calculating the average value for an order parameter in a particular part of the box

Suppose that you have calculated a vector of order parameters that can be assigned to a particular point in three dimensional space.
The symmetry functions in the [symfunc](module_symfunc.md) module are examples of order parameters that satisfy this criteria. You can use
the AROUND command to calculate the average value of the symmetry function in a particular part of the box as follows:

```plumed
f: FIXEDATOM AT=0,0,0
a: AROUND ...
  ATOMS=1-100 ORIGIN=f SIGMA=0.2
  XLOWER=-1.0 XUPPER=1.0
  YLOWER=-1.0 YUPPER=1.0
...

c: COORDINATIONNUMBER SPECIES=1-100 SWITCH={RATIONAL R_0=1.0} MASK=a
p: CUSTOM ARG=c,a FUNC=x*y PERIODIC=NO
n: SUM ARG=p PERIODIC=NO
d: SUM ARG=a PERIODIC=NO
av: CUSTOM ARG=n,d FUNC=x/y PERIODIC=NO
PRINT ARG=av FILE=colvar
```

The final quantity `av` here is:

$$
\overline{s}_{\tau} = \frac{ \sum_i c_i w(x_i,y_i,z_i) }{ \sum_i w(x_i,y_i,z_i) }
$$

where $c_i$ are the coordination numbers and $w_i$ is:

$$
w(x_i,y_i,z_i) = \int_{xl}^{xu} \int_{yl}^{yu} \textrm{d}x \textrm{d}y K\left( \frac{x - x_i}{\sigma} \right) K\left( \frac{y - y_i}{\sigma} \right)
$$

## Old syntax

In earlier versions of PLUMED the syntax for the calculation in the previous section is as follows:

```plumed
c: COORDINATIONNUMBER SPECIES=1-100 SWITCH={RATIONAL R_0=1.0}
f: FIXEDATOM AT=0,0,0
a: AROUND ...
  DATA=c ORIGIN=f SIGMA=0.2
  XLOWER=-1.0 XUPPER=1.0
  YLOWER=-1.0 YUPPER=1.0
  MEAN
...
PRINT ARG=a.mean FILE=colvar
```

This old syntax still works but __we highly recommend you use the newer syntax__ as it is easlier to understand,
more flexible and calculations with this newer input will run faster. You will notice that AROUND in the input above
is a shortcut that expands to a longer input
that is similar to that given in the previous section.  The main difference is that the order of the AROUND
and [COORDINATIONNUMBER](COORDINATIONNUMBER.md) actions is reversed in the new syntax.  The reason this reversal is necessary
is that the vector output by AROUND must be passed as as a MASK action to the [COORDINATIONNUMBER](COORDINATIONNUMBER.md)
action in order to optimize code performance.  Passing the vector from AROUND as a MASK in coordination number ensures that
we only calculate the coordination numbers for those atomms that are in the region of interest.  We thus avoid a lot of computational
expense that would otherwise be associated with calculating coordination numbers for atoms that are not within the region of
interest and would thus make no difference to the final average that we are calculating.

The old syntax also allowed you to compute the sum of the coordination numbers in the region of interest using an input like the one shown below:

```plumed
c: COORDINATIONNUMBER SPECIES=1-100 SWITCH={RATIONAL R_0=1.0}
f: FIXEDATOM AT=0,0,0
a: AROUND ...
  DATA=c ORIGIN=f SIGMA=0.2
  XLOWER=-1.0 XUPPER=1.0
  YLOWER=-1.0 YUPPER=1.0
  SUM
...
PRINT ARG=a.sum FILE=colvar
```

The final CV that is computed here is:

$$
\overline{s}_{\tau} = \sum_i c_i w(x_i,y_i,z_i)
$$

the equivalent input with the new syntax is thus:

```plumed
f: FIXEDATOM AT=0,0,0
a: AROUND ...
  ATOMS=1-100 ORIGIN=f SIGMA=0.2
  XLOWER=-1.0 XUPPER=1.0
  YLOWER=-1.0 YUPPER=1.0
...

c: COORDINATIONNUMBER SPECIES=1-100 SWITCH={RATIONAL R_0=1.0} MASK=a
p: CUSTOM ARG=c,a FUNC=x*y PERIODIC=NO
n: SUM ARG=p PERIODIC=NO
PRINT ARG=n FILE=colvar
```

That old syntax also allowed you to accumulate quantities such as:

$$
\overline{s}_{\tau} = \sum_i f(c_i) w(x_i,y_i,z_i)
$$

where $f$ could be one of the switching functions discussed in the documentation for [LESS_THAN](LESS_THAN.md),
one of the reverse switching functions discussed in the documentation for [MORE_THAN](MORE_THAN.md) or one of the
two sided switching functions discussed in the documentation for [BETWEEN](BETWEEN.md). An example of an old input
that computes all three of three types of sum is shown below:

```plumed
c: COORDINATIONNUMBER SPECIES=1-100 SWITCH={RATIONAL R_0=1.0}
f: FIXEDATOM AT=0,0,0
a: AROUND ...
  DATA=c ORIGIN=f SIGMA=0.2
  XLOWER=-1.0 XUPPER=1.0
  YLOWER=-1.0 YUPPER=1.0
  LESS_THAN={RATIONAL R_0=3}
  MORE_THAN={RATIONAL R_0=6}
  BETWEEN={GAUSSIAN LOWER=3 UPPER=6 SMEAR=0.5}
...
PRINT ARG=a.lessthan,a.between,a.morethan FILE=colvar
```

With the new syntax we can achieve the same result using the following input:

```plumed
f: FIXEDATOM AT=0,0,0
a: AROUND ...
  ATOMS=1-100 ORIGIN=f SIGMA=0.2
  XLOWER=-1.0 XUPPER=1.0
  YLOWER=-1.0 YUPPER=1.0
...

c: COORDINATIONNUMBER SPECIES=1-100 SWITCH={RATIONAL R_0=1.0} MASK=a
# This part does the LESS_THAN={RATIONAL R_0=3}
lt: LESS_THAN ARG=c SWITCH={RATIONAL R_0=3}
wlt: CUSTOM ARG=a,lt FUNC=x*y PERIODIC=NO
lessthan: SUM ARG=wlt PERIODIC=NO
# This part does the BETWEEN={GAUSSIAN LOWER=3 UPPER=6 SMEAR=0.5}
bt: BETWEEN ARG=c SWITCH={GAUSSIAN LOWER=3 UPPER=6 SMEAR=0.5}
wbt: CUSTOM ARG=a,bt FUNC=x*y PERIODIC=NO
between: SUM ARG=wbt PERIODIC=NO
# This part does the MORE_THAN={RATIONAL R_0=6}
mt: MORE_THAN ARG=c SWITCH={RATIONAL R_0=6}
wmt: CUSTOM ARG=a,mt FUNC=x*y PERIODIC=NO
morethan: SUM ARG=wmt PERIODIC=NO
PRINT ARG=lessthan,between,morethan FILE=colvar
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace volumes {
typedef ActionVolume<VolumeAround> Vola;
PLUMED_REGISTER_ACTION(Vola,"AROUND_CALC")
char glob_around[] = "AROUND";
typedef VolumeShortcut<glob_around> VolumeAroundShortcut;
PLUMED_REGISTER_ACTION(VolumeAroundShortcut,"AROUND")
}
}
