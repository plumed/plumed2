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
#include "MatrixDissimilarities.h"
#include "core/ActionRegister.h"
#include "core/AccelerableShortcut.h"

//+PLUMEDOC ANALYSIS DISSIMILARITIES
/*
Calculate the matrix of dissimilarities between a trajectory of atomic configurations.

This action allows you to calculate a dissimilarity matrix, which is a square matrix tht describes the
pairwise distinction between a set of objects. This action is routinely used in dimensionality reduction
calculations. The example shown below shows how we might get a matrix of dissimilarities upon which we can run a
dimensionality reduction calculation.

```plumed
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
# Collect the values of the two distances and store them for later analysis
ff: COLLECT_FRAMES ARG=d1,d2 STRIDE=1
# Transpose the matrix that is collected above
ff_dataT: TRANSPOSE ARG=ff_data
# Now compute all the dissimilarities between the collected frames
ss1: DISSIMILARITIES ARG=ff_data,ff_dataT
DUMPVECTOR ARG=ss1 FILE=mymatrix.dat
```

Some dimensionality reduction algorithms take the squares of the dissimilarities rather than the dissimilarities.
To calculate these quantities using PLUMED you use the SQUARED keyword as shown below:

```plumed
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
# Collect the values of the two distances and store them for later analysis
ff: COLLECT_FRAMES ARG=d1,d2 STRIDE=1
# Transpose the matrix that is collected above
ff_dataT: TRANSPOSE ARG=ff_data
# Now compute all the dissimilarities between the collected frames
ss1: DISSIMILARITIES ARG=ff_data,ff_dataT SQUARED
DUMPVECTOR ARG=ss1 FILE=mymatrix.dat
```

## The MASK keyword

The MASK keyword for DISSIMILARITIES works in the same way that it does for [MATRIX_PRODUCT](MATRIX_PRODUCT.md). This keyword thus
expects a matrix in input and will only evaluate elements the elements of the dissimilarity matrix whose corresponding elements in
the matrix that was input to MASK are non-zero. The following input shows an example of how this might be used in an example input
that has not been used in any production calculation.

```plumed
# Calculate and store three vectors connecting three pairs of atoms
d1: DISTANCE COMPONENTS ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6
ff1: VSTACK ARG=d1.x,d1.y,d1.z
# Calculate the magnitude of the vector and transform this magnitude using a
# switching function
d1m: CUSTOM ARG=d1.x,d1.y,d1.z FUNC=sqrt(x*x+y*y+z*z) PERIODIC=NO
sw1: LESS_THAN ARG=d1m SWITCH={RATIONAL R_0=0.2 D_MAX=0.5}
# Now perform the same operations for a different set of three atoms
d2: DISTANCE COMPONENTS ATOMS1=7,8 ATOMS2=9,10 ATOMS3=11,12
ff2: VSTACK ARG=d2.x,d2.y,d2.z
ff2T: TRANSPOSE ARG=ff2
d2m: CUSTOM ARG=d2.x,d2.y,d2.z FUNC=sqrt(x*x+y*y+z*z) PERIODIC=NO
sw2: LESS_THAN ARG=d2m SWITCH={RATIONAL R_0=0.2 D_MAX=0.5}
# Use OUTER_PRODUCT to work out which pairs of atoms are less than a certain cutoff
swmat: OUTER_PRODUCT ARG=sw1,sw2
# Now calculate the dissimilarities between the two sets of vectors for those pairs of atoms that are within a certain cutoff
d: DISSIMILARITIES ARG=ff1,ff2T MASK=swmat
mat: CUSTOM ARG=swmat,d FUNC=x*y PERIODIC=NO
# And compute the average dissimilarity
numer: SUM ARG=mat PERIODIC=NO
denom: SUM ARG=swmat PERIODIC=NO
v: CUSTOM ARG=numer,denom FUNC=x/y PERIODIC=NO
PRINT ARG=v FILE=colvar
```

We have not provided a more relevant example because we haven't really used this functionality for anything. If you have a better example
for a way of using this functionality please get in touch so we can update this part of the manual with your example.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace matrixtools {

typedef MatrixTimesMatrix<Dissimilarities> dissims;
PLUMED_REGISTER_ACTION(dissims,"DISSIMILARITIES_CPU")
typedef AccelerableShortcut<dissims> shortcut;
PLUMED_REGISTER_ACTION(shortcut,"DISSIMILARITIES")


}
}
