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
#include "MatrixProduct.h"
#include "core/ActionRegister.h"
#include "core/AccelerableShortcut.h"

//+PLUMEDOC MCOLVAR MATRIX_PRODUCT
/*
Calculate the product of two matrices

This action allows you to [multiply](https://en.wikipedia.org/wiki/Matrix_multiplication) two matrices. The following input shows an
example where two contact matrices are multiplied together.

```plumed
c1: CONTACT_MATRIX GROUP=1-100 SWITCH={RATIONAL R_0=0.1}
c2: CONTACT_MATRIX GROUP=1-100 SWITCH={RATIONAL R_0=0.2}
m: MATRIX_PRODUCT ARG=c1,c2
PRINT ARG=m FILE=colvar
```

The functionality in this action is useful for claculating the relative orientations of large numbers of molecules.  For example in
his input the [DISTANCE](DISTANCE.md) command is used to calculate the orientation of a collection of molecules.  We then can then use the [VSTACK](VSTACK.md), [TRANSPOSE](TRANSPOSE.md) and the
MATRIX_PRODUCT commands to calculate the dot products between all these vectors

```plumed
# Calculate the vectors connecting these three pairs of atoms
d: DISTANCE COMPONENTS ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6
# Construct a matrix that contains all the components of the vectors calculated
v: VSTACK ARG=d.x,d.y,d.z
# Transpose v
vT: TRANSPOSE ARG=v
# And now calculate the 3x3 matrix of dot products
m: MATRIX_PRODUCT ARG=v,vT
# And output the matrix product to a file
PRINT ARG=m FILE=colvar
```

## The MASK keyword

We can use MATRIX_PRODUCT to calculate a measure of local order as shown below:

```plumed
# Calculate the vectors connecting these three pairs of atoms
d: DISTANCE COMPONENTS ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6
# Normalize the distance vectors
dm: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6
dx: CUSTOM ARG=d.x,dm FUNC=x/y PERIODIC=NO
dy: CUSTOM ARG=d.y,dm FUNC=x/y PERIODIC=NO
dz: CUSTOM ARG=d.z,dm FUNC=x/y PERIODIC=NO
# Construct a matrix that contains all the directors of the vectors calculated
v: VSTACK ARG=dx,dy,dz
# Transpose v
vT: TRANSPOSE ARG=v
# Calculate a contact matrix between pairs of atoms
cmap: CONTACT_MATRIX GROUP=1,3,5 SWITCH={RATIONAL R_0=0.2 D_MAX=0.5}
# Calculate the matrix of dot products
prod: MATRIX_PRODUCT ARG=v,vT MASK=cmap
# Take the product of the two matrices we computed above
p: CUSTOM ARG=cmap,prod FUNC=x*y PERIODIC=NO
# And compute the average of the dot product for the molecules in the
# first coordination sphere around each of the atoms
ones: ONES SIZE=3
numer: MATRIX_VECTOR_PRODUCT ARG=p,ones
denom: MATRIX_VECTOR_PRODUCT ARG=cmap,ones
order: CUSTOM ARG=numer,denom FUNC=x/y PERIODIC=NO
# Print the order parameter values
PRINT ARG=order FILE=colvar
```

In the above input the [DISTANCE](DISTANCE.md) action is used to calculate the orientations of some molecules (you would normally use more distances than the three
that are used in this input when doing this sort of calculation).  We then construct a matrix called `v` that contains the directors of the vectors that define the
orientation of these molecules.  The final vector of values `order` that is output here measures the average for the dot product between the orientations of the molecules
and the molecules in their first coordination sphere.

Ideass similar to this are used in the calculation of [LOCAL_Q3](LOCAL_Q3.md), [LOCAL_Q4](LOCAL_Q4.md) and [LOCAL_Q6](LOCAL_Q6.md).  Notice, that it is very important to
use `MASK` keyword in the MATRIX_PRODUCT action when you are doing this type calculation.  Using this keyword ensures that PLUMED does not calculate the $(i,j)$ matrix of
the matrix `prod` if the corresponding element in `cmap` is zero.  Using this keyword thus ensures that we can avoid a large number of unecessary calculations and improves
the performance of the calculation considerably.

## Calculating angles

MATRIX_PRODUCT can also be used to calculate a matrix of angles between bonds as shown below:

```plumed
# Calculate the directors for a set of vectors
d: DISTANCE COMPONENTS ATOMS1=1,2 ATOMS2=1,3 ATOMS3=1,4 ATOMS4=1,5 ATOMS5=1,6
dm: DISTANCE ATOMS1=1,2 ATOMS2=1,3 ATOMS3=1,4 ATOMS4=1,5 ATOMS5=1,6
dx: CUSTOM ARG=d.x,dm FUNC=x/y PERIODIC=NO
dy: CUSTOM ARG=d.y,dm FUNC=x/y PERIODIC=NO
dz: CUSTOM ARG=d.z,dm FUNC=x/y PERIODIC=NO
# Construct a matrix that contains all the directors of the vectors calculated
v: VSTACK ARG=dx,dy,dz
# Transpose v
vT: TRANSPOSE ARG=v
# Calculate the matrix of dot products between the input directors
dpmat: MATRIX_PRODUCT ELEMENTS_ON_DIAGONAL_ARE_ZERO ARG=v,vT
# And calculate the angles
angles: CUSTOM ARG=dpmat FUNC=acos(x) PERIODIC=NO
# Print the matrix of angles
PRINT ARG=angles FILE=colvar
```

Notice that we have to use the `ELEMENTS_ON_DIAGONAL_ARE_ZERO` flag here to avoid numerical issues in the calculation.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace matrixtools {

typedef MatrixTimesMatrix<MatrixProduct> mtimes;
PLUMED_REGISTER_ACTION(mtimes,"MATRIX_PRODUCT_CPU")
typedef AccelerableShortcut<mtimes> shortcut;
PLUMED_REGISTER_ACTION(shortcut,"MATRIX_PRODUCT")


}
}
