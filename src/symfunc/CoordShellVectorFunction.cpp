/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2020 The plumed team
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
#include "CoordinationNumbers.h"
#include "multicolvar/MultiColvarShortcuts.h"
#include "core/ActionWithValue.h"
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

#include <string>
#include <cmath>

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVAR COORDINATION_SHELL_FUNCTION
/*
Calculate an arbitrary function of all the bond vectors in the first coordination sphere of an atom

This shortcut allows you to calculate the sum for an arbitrary function of the bond vectors that connect an atom to each of its neighbours.
In other words, this action allows you to compute the following:

$$
s_i = \sum_{i \ne j} \sigma(r_{ij}) f(x_{ij}, y_{ij}, z_{ij}, r_{ij}) )
$$

In this expression, $x_{ij}, y_{ij}, z_{ij}$ are the components of the vector connecting atoms $i$ and $j$ and $r_{ij}$ is the magnitude of this vector.
$\sigma(r_{ij})$ is then a switching function that ensures that the aritrary function $f$ is only evaluated for if atom $j$ is within a certain cutoff distance
of atom $i$.

Below you can see a simple example that shows how this action can be used in practise.

```plumed
cv: COORDINATION_SHELL_FUNCTION SPECIES=1-64 D_0=3.0 R_0=1.5 FUNCTION=((x+y+z)/r)^3+((x-y-z)/r)^3+((-x+y-z)/r)^3+((-x-y+z)/r)^3
PRINT ARG=cv FILE=colvar
```

The above input calculates 64 $s_i$ values - one $s_i$ values for each of the atoms specified using the SPECIES keyword.  These 64 numbers are then output to a file.
The function of x, y, z and r to be evaluated is specified using the FUNCTION keyword.  Obviously, if your function does not depend on all four of these variables
they can be excluded from your function.

In the input above we use a rational [switching function](LESS_THAN.md) with the parameters above. We would recommend using SWITCH syntax
rather than the syntax above when giving the parameters for the switching function as you can then use any of the switching functions described
in the documentation for [LESS_THAN](LESS_THAN.md).  More importantly, however, using this syntax allows you to set the D_MAX parameter for the
switching function as demonstrated below:

```plumed
cv: COORDINATION_SHELL_FUNCTION ...
   SPECIES=1-64 SWITCH={RATIONAL D_0=3.0 R_0=1.5 D_MAX=7.0}
   FUNCTION=((x+y+z)/r)^3+((x-y-z)/r)^3+((-x+y-z)/r)^3+((-x-y+z)/r)^3
...
PRINT ARG=cv FILE=colvar
```

Setting the `D_MAX` can substantially improve PLUMED performance as it turns on the linked list algorithm that is discussed in the optimisation details part
of the documentation for [CONTACT_MATRIX](CONTACT_MATRIX.md).

## Rotating the reference frame

As discussed in [this paper](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.92.180102) it is sometimes useful to rotate the bond vectors before computing the
arbitrary function $f$ in the above expression.  To perform such rotations you use the PHI, THETA and PSI keywords.  The $x_{ij}, y_{ij}$ and $z_{ij}$ values that enter $f$ in the
expression above are then calculated as:

$$
\left(
\begin{matrix}
x_{ij} \\
y_{ij} \\
z_{ij}
\end{matrix}
\right) =
\left(
\begin{matrix}
\cos(\psi)\cos(\phi) - \cos(\theta)\sin(\phi)\sin(\psi) & \cos(\psi)*\sin(\phi)+\cos(\theta)*\cos(\phi)*\sin(\psi) & \sin(\psi)*sin(\theta) \\
-\sin(\psi)*\cos(\phi)-\cos(\theta)*\sin(\phi)*\cos(\psi) & -\sin(\psi)*\sin(\phi)+\cos(\theta)*\cos(\phi)*std::cos(\psi), & \cos(\psi)*\sin(\theta) \\
\sin(\theta)*\sin(\phi) & \sin(\theta)*\cos(\phi) & \cos(\theta)
\end{matrix}
\right)
\left(
\begin{matrix}
x_{ij}' \\
y_{ij}' \\
z_{ij}'
\end{matrix}
\right)
$$

$x_{ij}', y_{ij}'$ and $z_{ij}'$ in this expression are the bond vectors that are calculated in the lab frame.  The matrix in the above expression is thus just a
[rotation matrix](https://en.wikipedia.org/wiki/Rotation_matrix) that converts the lab frame to some frame of interest.

## Working with two types of atom

If you would like to calculate a function of the bonds connecting the atoms in GROUPA to the atoms in GROUPB you can use an input like the one
shown below:

```plumed
d: COORDINATION_SHELL_FUNCTION ...
   SPECIESA=1-64 SPECIESB=65-200
   FUNCTION=((x+y+z)/r)^3+((x-y-z)/r)^3+((-x+y-z)/r)^3+((-x-y+z)/r)^3
   SWITCH={RATIONAL D_0=3.0 R_0=1.5 D_MAX=6.0}
...
lt: MORE_THAN ARG=d SWITCH={RATIONAL R_0=0.5}
s: SUM ARG=lt PERIODIC=NO
PRINT ARG=s FILE=colv
```

## The MASK keyword

You can use the MASK keyword with this action in the same way that it is used with [COORDINATIONNUMBER](COORDINATIONNUMBER.md).  This keyword thus expects a vector in
input, which tells PLUMED the atoms for which you do not need to calculate the function.  As illustrated below, this is useful if you are using functionality
from the [volumes module](module_volumes.md) to calculate the average value of the TETRAHEDRAL parameter for only those atoms that lie in a certain part of the simulation box.

```plumed
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=1-400 CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# Calculate the tetrahedral parameter of the atoms
d: COORDINATION_SHELL_FUNCTION ...
   SPECIES=1-400 MASK=sphere
   FUNCTION=((x+y+z)/r)^3+((x-y-z)/r)^3+((-x+y-z)/r)^3+((-x-y+z)/r)^3
   SWITCH={RATIONAL D_0=3.0 R_0=1.5 D_MAX=6.0}
...
# Multiply fccubic parameters numbers by sphere vector
prod: CUSTOM ARG=d,sphere FUNC=x*y PERIODIC=NO
# Sum of coordination numbers for atoms that are in the sphere of interest
numer: SUM ARG=prod PERIODIC=NO
# Number of atoms that are in sphere of interest
denom: SUM ARG=sphere PERIODIC=NO
# Average coordination number for atoms in sphere of interest
av: CUSTOM ARG=numer,denom FUNC=x/y PERIODIC=NO
# And print out final CV to a file
PRINT ARG=av FILE=colvar STRIDE=1
```

This input calculate the average value of the function for only those atoms that are within a spherical region that is centered on the point
$(2.5,2.5,2.5)$.

## Deprecated syntax

More information on the deprecated keywords that are given below is available in the documentation for the [DISTANCES](DISTANCES.md) command.

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR COORDINATION_SHELL_AVERAGE
/*
Calculate an arbitrary function of all the bond vectors in the first coordination sphere of an atom and take an average

This shortcut allows you to calculate the average for an arbitrary function of the bond vectors that connect an atom to each of its neighbours.
In other words, this action allows you to compute the following:

$$
s_i = \frac{\sum_{i \ne j} \sigma(r_{ij}) f(x_{ij}, y_{ij}, z_{ij}, r_{ij}) )}{\sum_{i \ne j} \sigma(r_{ij})}
$$

In this expression, $x_{ij}, y_{ij}, z_{ij}$ are the components of the vector connecting atoms $i$ and $j$ and $r_{ij}$ is the magnitude of this vector.
$\sigma(r_{ij})$ is then a switching function that ensures that the aritrary function $f$ is only evaluated for if atom $j$ is within a certain cutoff distance
of atom $i$.

Below you can see a simple example that shows how this action can be used in practise.

```plumed
cv: COORDINATION_SHELL_AVERAGE SPECIES=1-64 D_0=3.0 R_0=1.5 FUNCTION=((x+y+z)/r)^3+((x-y-z)/r)^3+((-x+y-z)/r)^3+((-x-y+z)/r)^3
PRINT ARG=cv FILE=colvar
```

The above input calculates 64 $s_i$ values - one $s_i$ values for each of the atoms specified using the SPECIES keyword.  These 64 numbers are then output to a file.
The function of x, y, z and r to be evaluated is specified using the FUNCTION keyword.  Obviously, if your function does not depend on all four of these variables
they can be excluded from your function.

In the input above we use a rational [switching function](LESS_THAN.md) with the parameters above. We would recommend using SWITCH syntax
rather than the syntax above when giving the parameters for the switching function as you can then use any of the switching functions described
in the documentation for [LESS_THAN](LESS_THAN.md).  More importantly, however, using this syntax allows you to set the D_MAX parameter for the
switching function as demonstrated below:

```plumed
cv: COORDINATION_SHELL_AVERAGE ...
   SPECIES=1-64 SWITCH={RATIONAL D_0=3.0 R_0=1.5 D_MAX=7.0}
   FUNCTION=((x+y+z)/r)^3+((x-y-z)/r)^3+((-x+y-z)/r)^3+((-x-y+z)/r)^3
...
PRINT ARG=cv FILE=colvar
```

Setting the `D_MAX` can substantially improve PLUMED performance as it turns on the linked list algorithm that is discussed in the optimisation details part
of the documentation for [CONTACT_MATRIX](CONTACT_MATRIX.md).

Notice that you can you can rotate the bond vectors before computing the
arbitrary function $f$ in the above expression as is discussed in the documentation for [COORDINATION_SHELL_FUNCTION](COORDINATION_SHELL_FUNCTION.md)

## Working with two types of atom

If you would like to calculate a function of the bonds connecting the atoms in GROUPA to the atoms in GROUPB you can use an input like the one
shown below:

```plumed
d: COORDINATION_SHELL_AVERAGE ...
   SPECIESA=1-64 SPECIESB=65-200
   FUNCTION=((x+y+z)/r)^3+((x-y-z)/r)^3+((-x+y-z)/r)^3+((-x-y+z)/r)^3
   SWITCH={RATIONAL D_0=3.0 R_0=1.5 D_MAX=6.0}
...
lt: MORE_THAN ARG=d SWITCH={RATIONAL R_0=0.5}
s: SUM ARG=lt PERIODIC=NO
PRINT ARG=s FILE=colv
```

## The MASK keyword

You can use the MASK keyword with this action in the same way that it is used with [COORDINATIONNUMBER](COORDINATIONNUMBER.md).  This keyword thus expects a vector in
input, which tells PLUMED the atoms for which you do not need to calculate the function.  As illustrated below, this is useful if you are using functionality
from the [volumes module](module_volumes.md) to calculate the average value of the TETRAHEDRAL parameter for only those atoms that lie in a certain part of the simulation box.

```plumed
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=1-400 CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# Calculate the function for the atoms
d: COORDINATION_SHELL_AVERAGE ...
   SPECIES=1-400 MASK=sphere
   FUNCTION=((x+y+z)/r)^3+((x-y-z)/r)^3+((-x+y-z)/r)^3+((-x-y+z)/r)^3
   SWITCH={RATIONAL D_0=3.0 R_0=1.5 D_MAX=6.0}
...
# Multiply fccubic parameters numbers by sphere vector
prod: CUSTOM ARG=d,sphere FUNC=x*y PERIODIC=NO
# Sum of coordination numbers for atoms that are in the sphere of interest
numer: SUM ARG=prod PERIODIC=NO
# Number of atoms that are in sphere of interest
denom: SUM ARG=sphere PERIODIC=NO
# Average coordination number for atoms in sphere of interest
av: CUSTOM ARG=numer,denom FUNC=x/y PERIODIC=NO
# And print out final CV to a file
PRINT ARG=av FILE=colvar STRIDE=1
```

This input calculate the average value of the function for only those atoms that are within a spherical region that is centered on the point
$(2.5,2.5,2.5)$.

## Deprecated syntax

More information on the deprecated keywords that are given below is available in the documentation for the [DISTANCES](DISTANCES.md) command.

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR SIMPLECUBIC
/*
Calculate whether or not the coordination spheres of atoms are arranged as they would be in a simple cubic structure.

This shortcut is an example of a [COORDINATION_SHELL_AVERAGE](COORDINATION_SHELL_AVERAGE.md),
which we can use to measure how similar the environment around atom $i$ is to a simple cubic structure.  We perform this comparison by evaluating
the following quantity:

$$
s_i = \frac{ \sum_{i \ne j} \sigma(r_{ij}) \left[ \frac{ x_{ij}^4 + y_{ij}^4 + z_{ij}^4 }{r_{ij}^4} \right] }{ \sum_{i \ne j} \sigma(r_{ij}) }
$$

In this expression $x_{ij}$, $y_{ij}$ and $z_{ij}$ are the $x$, $y$ and $z$ components of the vector connecting atom $i$ to
atom $j$ and $r_{ij}$ is the magnitude of this vector.  $\sigma(r_{ij})$ is a switching function that acts on the distance between atom $i$ and atom $j$ and its inclusion in the numerator and the denominator of the above expression as well as the fact that we are summing
over all of the other atoms in the system ensures that we are calculating an average
of the function of $x_{ij}$, $y_{ij}$ and $z_{ij}$ for the atoms in the first coordination sphere around atom $i$.
This quantity is once again a multicolvar so you can compute it for multiple atoms using a single PLUMED action and then compute
the average value for the atoms in your system, the number of atoms that have an $s_i$ value that is more that some target and
so on.  Notice also that you can rotate the reference frame if you are using a non-standard unit cell.

The following input tells plumed to calculate the simple cubic parameter for the atoms 1-100 with themselves.
The mean value is then calculated.

```plumed
sc: SIMPLECUBIC SPECIES=1-100 R_0=1.0
sc_mean: MEAN ARG=sc PERIODIC=NO
PRINT ARG=sc_mean FILE=colvar
```

In the input above we use a rational [switching function](LESS_THAN.md) with the parameters above. We would recommend using SWITCH syntax
rather than the syntax above when giving the parameters for the switching function as you can then use any of the switching functions described
in the documentation for [LESS_THAN](LESS_THAN.md).  More importantly, however, using this syntax allows you to set the D_MAX parameter for the
switching function as demonstrated below:

```plumed
d: SIMPLECUBIC SPECIES=1-64 SWITCH={RATIONAL D_0=3.0 R_0=1.5 D_MAX=6.0}
dm: MEAN ARG=d PERIODIC=NO
PRINT ARG=dm FILE=colv
```

Setting the `D_MAX` can substantially improve PLUMED performance as it turns on the linked list algorithm that is discussed in the optimisation details part
of the documentation for [CONTACT_MATRIX](CONTACT_MATRIX.md).

## Working with two types of atom

If you would like to calculate whether the atoms in GROUPB are arranged around the atoms in GROUPA as they in an tetrahedral structure you use an input like the one
shown below:

```plumed
d: SIMPLECUBIC SPECIESA=1-64 SPECIESB=65-200 SWITCH={RATIONAL D_0=3.0 R_0=1.5 D_MAX=6.0}
lt: MORE_THAN ARG=d SWITCH={RATIONAL R_0=0.5}
s: SUM ARG=lt PERIODIC=NO
PRINT ARG=s FILE=colv
```

## The MASK keyword

You can use the MASK keyword with this action in the same way that it is used with [COORDINATIONNUMBER](COORDINATIONNUMBER.md).  This keyword thus expects a vector in
input, which tells SIMPLECUBIC that it is safe not to calculate the SIMPLECUBIC parameter for some of the atoms.  As illustrated below, this is useful if you are using functionality
from the [volumes module](module_volumes.md) to calculate the average value of the SIMPLECUBIC parameter for only those atoms that lie in a certain part of the simulation box.

```plumed
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=1-400 CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# Calculate the tetrahedral parameter of the atoms
cc: SIMPLECUBIC ...
  SPECIES=1-400 MASK=sphere
  SWITCH={RATIONAL D_0=3.0 R_0=1.5 D_MAX=6.0}
...
# Multiply fccubic parameters numbers by sphere vector
prod: CUSTOM ARG=cc,sphere FUNC=x*y PERIODIC=NO
# Sum of coordination numbers for atoms that are in the sphere of interest
numer: SUM ARG=prod PERIODIC=NO
# Number of atoms that are in sphere of interest
denom: SUM ARG=sphere PERIODIC=NO
# Average coordination number for atoms in sphere of interest
av: CUSTOM ARG=numer,denom FUNC=x/y PERIODIC=NO
# And print out final CV to a file
PRINT ARG=av FILE=colvar STRIDE=1
```

This input calculate the average value of the SIMPLECUBIC parameter for only those atoms that are within a spherical region that is centered on the point
$(2.5,2.5,2.5)$.

## Deprecated syntax

More information on the deprecated keywords that are given below is available in the documentation for the [DISTANCES](DISTANCES.md) command.

The following command illustrates how this deprecated syntax works. It tells plumed to look at the ways atoms 1-100 are within 3.0 are arranged about atoms
from 101-110.  The number of simple cubic parameters that are greater than 0.8 is then output

```plumed
SIMPLECUBIC SPECIESA=101-110 SPECIESB=1-100 R_0=3.0 MORE_THAN={RATIONAL R_0=0.8 NN=6 MM=12 D_0=0}
```

Notice that you can you can rotate the bond vectors before computing the
function in the above expression as is discussed in the documentation for [COORDINATION_SHELL_FUNCTION](COORDINATION_SHELL_FUNCTION.md)

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR TETRAHEDRAL
/*
Calculate the degree to which the environment about ions has a tetrahedral order.

This shortcut is an example of a [COORDINATION_SHELL_AVERAGE](COORDINATION_SHELL_AVERAGE.md),
which we can use to measure the degree to which the atoms in the first coordination shell around any atom, $i$ is
is arranged like a tetrahedron.  We perform this comparison by evaluating the following function.

$$
 s(i) = \frac{1}{\sum_j \sigma( r_{ij} )} \sum_j \sigma( r_{ij} )\left[ \frac{(x_{ij} + y_{ij} + z_{ij})^3}{r_{ij}^3} +
                                                                        \frac{(x_{ij} - y_{ij} - z_{ij})^3}{r_{ij}^3} +
                                                                        \frac{(-x_{ij} + y_{ij} - z_{ij})^3}{r_{ij}^3} +
                                                                        \frac{(-x_{ij} - y_{ij} + z_{ij})^3}{r_{ij}^3} \right]
$$

Here $r_{ij}$ is the magnitude of the vector connecting atom $i$ to atom $j$ and $x_{ij}$, $y_{ij}$ and $z_{ij}$
are its three components.  The function  $\sigma( r_{ij} )$ is a switching function that acts on the distance between
atoms $i$ and $j$.  The parameters of this function should be set so that the function is equal to one
when atom $j$ is in the first coordination sphere of atom $i$ and is zero otherwise.

The following command calculates the average value of the TETRAHEDRAL parameter for a set of 64 atoms all of the same type
and outputs this quantity to a file called colvar.

```plumed
tt: TETRAHEDRAL SPECIES=1-64 D_0=3.0 R_0=1.5 NN=6 MM=12
tt_mean: MEAN ARG=tt PERIODIC=NO
PRINT ARG=tt_mean FILE=colvar
```

In the input above we use a rational [switching function](LESS_THAN.md) with the parameters above. We would recommend using SWITCH syntax
rather than the syntax above when giving the parameters for the switching function as you can then use any of the switching functions described
in the documentation for [LESS_THAN](LESS_THAN.md).  More importantly, however, using this syntax allows you to set the D_MAX parameter for the
switching function as demonstrated below:

```plumed
d: TETRAHEDRAL SPECIES=1-64 SWITCH={RATIONAL D_0=3.0 R_0=1.5 D_MAX=6.0}
dm: MEAN ARG=d PERIODIC=NO
PRINT ARG=dm FILE=colv
```

Setting the `D_MAX` can substantially improve PLUMED performance as it turns on the linked list algorithm that is discussed in the optimisation details part
of the documentation for [CONTACT_MATRIX](CONTACT_MATRIX.md).

## Working with two types of atom

If you would like to calculate whether the atoms in GROUPB are arranged around the atoms in GROUPA as they in an tetrahedral structure you use an input like the one
shown below:

```plumed
d: TETRAHEDRAL SPECIESA=1-64 SPECIESB=65-200 SWITCH={RATIONAL D_0=3.0 R_0=1.5 D_MAX=6.0}
lt: MORE_THAN ARG=d SWITCH={RATIONAL R_0=0.5}
s: SUM ARG=lt PERIODIC=NO
PRINT ARG=s FILE=colv
```

## The MASK keyword

You can use the MASK keyword with this action in the same way that it is used with [COORDINATIONNUMBER](COORDINATIONNUMBER.md).  This keyword thus expects a vector in
input, which tells TETRAHEDRAL that it is safe not to calculate the TETRAHEDRAL parameter for some of the atoms.  As illustrated below, this is useful if you are using functionality
from the [volumes module](module_volumes.md) to calculate the average value of the TETRAHEDRAL parameter for only those atoms that lie in a certain part of the simulation box.

```plumed
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=1-400 CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# Calculate the tetrahedral parameter of the atoms
cc: TETRAHEDRAL ...
  SPECIES=1-400 MASK=sphere
  SWITCH={RATIONAL D_0=3.0 R_0=1.5 D_MAX=6.0}
...
# Multiply fccubic parameters numbers by sphere vector
prod: CUSTOM ARG=cc,sphere FUNC=x*y PERIODIC=NO
# Sum of coordination numbers for atoms that are in the sphere of interest
numer: SUM ARG=prod PERIODIC=NO
# Number of atoms that are in sphere of interest
denom: SUM ARG=sphere PERIODIC=NO
# Average coordination number for atoms in sphere of interest
av: CUSTOM ARG=numer,denom FUNC=x/y PERIODIC=NO
# And print out final CV to a file
PRINT ARG=av FILE=colvar STRIDE=1
```

This input calculate the average value of the TETRAHEDRAL parameter for only those atoms that are within a spherical region that is centered on the point
$(2.5,2.5,2.5)$.

## Deprecated syntax

More information on the deprecated keywords that are given below is available in the documentation for the [DISTANCES](DISTANCES.md) command.

The following command illustrates how this deprecated syntax works.  It calculates the number of TETRAHEDRAL parameters that are greater than 0.8 in a set of 10 atoms.
In this calculation it is assumed that there are two atom types A and B and that the first coordination sphere of the
10 atoms of type A contains atoms of type B.  The formula above is thus calculated for ten different A atoms and within
it the sum over $j$ runs over 40 atoms of type B that could be in the first coordination sphere.

```plumed
tt: TETRAHEDRAL SPECIESA=1-10 SPECIESB=11-40 SWITCH={RATIONAL D_0=1.3 R_0=0.2} MORE_THAN={RATIONAL R_0=0.8}
PRINT ARG=tt.* FILE=colvar
```

Notice that you can you can rotate the bond vectors before computing the
function in the above expression as is discussed in the documentation for [COORDINATION_SHELL_FUNCTION](COORDINATION_SHELL_FUNCTION.md)

*/
//+ENDPLUMEDOC

class CoordShellVectorFunction : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit CoordShellVectorFunction(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(CoordShellVectorFunction,"FCCUBIC")
PLUMED_REGISTER_ACTION(CoordShellVectorFunction,"TETRAHEDRAL")
PLUMED_REGISTER_ACTION(CoordShellVectorFunction,"SIMPLECUBIC")
PLUMED_REGISTER_ACTION(CoordShellVectorFunction,"COORDINATION_SHELL_FUNCTION")
PLUMED_REGISTER_ACTION(CoordShellVectorFunction,"COORDINATION_SHELL_AVERAGE")

void CoordShellVectorFunction::registerKeywords( Keywords& keys ) {
  CoordinationNumbers::shortcutKeywords( keys );
  if( keys.getDisplayName()=="COORDINATION_SHELL_FUNCTION" || keys.getDisplayName()=="COORDINATION_SHELL_AVERAGE" ) {
    keys.add("compulsory","FUNCTION","the function of the bond vectors that you would like to evaluate");
  }
  keys.add("compulsory","PHI","0.0","The Euler rotational angle phi");
  keys.add("compulsory","THETA","0.0","The Euler rotational angle theta");
  keys.add("compulsory","PSI","0.0","The Euler rotational angle psi");
  if( keys.getDisplayName()=="FCCUBIC" ) {
    keys.add("compulsory","ALPHA","3.0","The alpha parameter of the angular function that is used for FCCUBIC");
  }
  keys.addDeprecatedFlag("LOWMEM","");
  keys.setValueDescription("vector","the symmetry function for each of the specified atoms");
  keys.needsAction("CONTACT_MATRIX");
  keys.needsAction("FCCUBIC_FUNC");
  keys.needsAction("CUSTOM");
  keys.needsAction("ONES");
  keys.needsAction("MATRIX_VECTOR_PRODUCT");
  keys.addFlag("USEGPU",false,"run part of this calculation on the GPU");
  keys.addDOI("10.1103/PhysRevB.81.125416");
  keys.addDOI("10.1103/PhysRevB.92.180102");
  keys.addDOI("10.1063/1.4997180");
  keys.addDOI("10.1063/1.5134461");
}

CoordShellVectorFunction::CoordShellVectorFunction(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
#define createLabel(name) const std::string name##Lab = getShortcutLabel()+"_"#name;
  bool lowmem;
  parseFlag("LOWMEM",lowmem);
  if( lowmem ) {
    warning("LOWMEM flag is deprecated and is no longer required for this action");
  }
  bool usegpu;
  parseFlag("USEGPU",usegpu);
  const std::string doUSEGPU = usegpu?" USEGPU":"";
  std::string sp_str;
  parse("SPECIES",sp_str);
  std::string specA;
  parse("SPECIESA",specA);
  std::string specB;
  parse("SPECIESB",specB);
  createLabel(mat);
  if( sp_str.length()>0 || specA.length()>0 ) {
    CoordinationNumbers::expandMatrix( true, getShortcutLabel(),  sp_str, specA, specB, this );
  } else {
    error("found no input atoms use SPECIES/SPECIESA");
  }
  double phi;
  parse("PHI",phi);
  double theta;
  parse("THETA",theta);
  double psi;
  parse("PSI",psi);
  std::vector<std::string> rotelements(9);
  std::string xvec = matLab + ".x", yvec = matLab + ".y", zvec = matLab + ".z";
  if( phi!=0 || theta!=0 || psi!=0 ) {
    Tools::convert( std::cos(psi)*std::cos(phi)
                    - std::cos(theta)*std::sin(phi)*std::sin(psi),
                    rotelements[0] );
    Tools::convert( std::cos(psi)*std::sin(phi)
                    + std::cos(theta)*std::cos(phi)*std::sin(psi),
                    rotelements[1] );
    Tools::convert( std::sin(psi)*std::sin(theta), rotelements[2] );

    Tools::convert( -std::sin(psi)*std::cos(phi)
                    - std::cos(theta)*std::sin(phi)*std::cos(psi),
                    rotelements[3] );
    Tools::convert( -std::sin(psi)*std::sin(phi)
                    + std::cos(theta)*std::cos(phi)*std::cos(psi),
                    rotelements[4] );
    Tools::convert( std::cos(psi)*std::sin(theta), rotelements[5] );

    Tools::convert( std::sin(theta)*std::sin(phi), rotelements[6] );
    Tools::convert( -std::sin(theta)*std::cos(phi), rotelements[7] );
    Tools::convert( std::cos(theta), rotelements[8] );
    createLabel(xrot);
    createLabel(yrot);
    createLabel(zrot);
    readInputLine( xrotLab + ": CUSTOM ARG=" + matLab + ".x," + matLab + ".y," + matLab + ".z "
                   "FUNC=" + rotelements[0] + "*x+" + rotelements[1] + "*y+" + rotelements[2] + "*z "
                   "PERIODIC=NO");
    readInputLine( yrotLab + ": CUSTOM ARG=" + matLab + ".x," + matLab + ".y," + matLab + ".z "
                   "FUNC=" + rotelements[3] + "*x+" + rotelements[4] + "*y+" + rotelements[5] + "*z "
                   "PERIODIC=NO");
    readInputLine( zrotLab + ": CUSTOM ARG=" + matLab + ".x," + matLab + ".y," + matLab + ".z "
                   "FUNC=" + rotelements[6] + "*x+" + rotelements[7] + "*y+" + rotelements[8] + "*z "
                   "PERIODIC=NO");
  }
  // Calculate FCC cubic function from bond vectors
  createLabel(vfunc)
  createLabel(r);
  if( getName()=="FCCUBIC" ) {
    std::string alpha;
    parse("ALPHA",alpha);
    readInputLine( vfuncLab + ": FCCUBIC_FUNC "
                   "ARG=" + xvec + "," + yvec + "," + zvec+ " ALPHA=" + alpha);
  } else if( getName()=="TETRAHEDRAL" ) {
    readInputLine( rLab + ": CUSTOM ARG=" + xvec + "," + yvec + "," + zvec + " "
                   "PERIODIC=NO FUNC=sqrt(x*x+y*y+z*z)");
    readInputLine( vfuncLab + ": CUSTOM "
                   "ARG=" + xvec + "," + yvec + "," + zvec + "," + rLab
                   + " VAR=x,y,z,r PERIODIC=NO "
                   "FUNC=((x+y+z)/r)^3+((x-y-z)/r)^3+((-x+y-z)/r)^3+((-x-y+z)/r)^3" );
  } else if( getName()=="SIMPLECUBIC" ) {
    readInputLine( rLab + ": CUSTOM ARG=" + xvec + "," + yvec + "," + zvec
                   + " PERIODIC=NO FUNC=sqrt(x*x+y*y+z*z)");
    readInputLine( vfuncLab + ": CUSTOM "
                   "ARG=" + xvec + "," + yvec + "," + zvec + "," + rLab
                   + " VAR=x,y,z,r PERIODIC=NO FUNC=(x^4+y^4+z^4)/(r^4)" );
  } else {
    std::string myfunc;
    parse("FUNCTION",myfunc);
    if( myfunc.find("r")!=std::string::npos ) {
      readInputLine( rLab + ": CUSTOM ARG=" + xvec + "," + yvec + "," + zvec
                     + " PERIODIC=NO FUNC=sqrt(x*x+y*y+z*z)");
      readInputLine( vfuncLab + ": CUSTOM "
                     "ARG=" + xvec + "," + yvec + "," + zvec + "," + rLab
                     + " VAR=x,y,z,r PERIODIC=NO FUNC=" + myfunc );
    } else {
      readInputLine( vfuncLab + ": CUSTOM ARG=" + xvec + "," + yvec + "," + zvec
                     + " PERIODIC=NO FUNC=" + myfunc );
    }
  }
  // Hadamard product of function above and weights
  createLabel(wvfunc);
  readInputLine( wvfuncLab + ": CUSTOM ARG=" + vfuncLab + "," + matLab + ".w "
                 "FUNC=x*y PERIODIC=NO");
  // And coordination numbers
  ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( matLab );
  plumed_assert( av && av->getNumberOfComponents()>0 && (av->copyOutput(0))->getRank()==2 );
  std::string size;
  Tools::convert( (av->copyOutput(0))->getShape()[1], size );
  createLabel(ones);
  readInputLine( onesLab + ": ONES SIZE=" + size );
  readInputLine( getShortcutLabel() + ": MATRIX_VECTOR_PRODUCT "
                 "ARG=" + wvfuncLab + "," + onesLab + doUSEGPU);
  std::string olab=getShortcutLabel();
  if( getName()!="COORDINATION_SHELL_FUNCTION" ) {
    olab = getShortcutLabel() + "_n";
    // Calculate coordination numbers for denominator
    createLabel(denom);
    readInputLine( denomLab + ": MATRIX_VECTOR_PRODUCT "
                   "ARG=" + matLab + ".w," + onesLab + doUSEGPU);
    // And normalise
    readInputLine( olab + ": CUSTOM "
                   "ARG=" + getShortcutLabel() + "," + denomLab + " FUNC=x/y PERIODIC=NO");
  }
  // And expand the functions
  std::map<std::string,std::string> keymap;
  multicolvar::MultiColvarShortcuts::readShortcutKeywords( keymap, this );
  multicolvar::MultiColvarShortcuts::expandFunctions( getShortcutLabel(),
      olab, "", keymap, this );
#undef createLabel
}

}
}

