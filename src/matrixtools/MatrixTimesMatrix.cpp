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
#ifdef __PLUMED_HAS_OPENACC
#define __PLUMED_USE_OPENACC 1
#endif //__PLUMED_HAS_OPENACC
#include "MatrixTimesMatrix.h"
#include "core/ActionRegister.h"

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

struct MatrixProduct {
  static void registerKeywords( Keywords& keys );
  void setup( MatrixTimesMatrix<MatrixProduct>* action, const Value* myval ) {}
  static void calculate( bool noderiv,
                         const MatrixProduct& actdata,
                         const InputVectors& vectors,
                         MatrixElementOutput& output );
#ifdef __PLUMED_USE_OPENACC
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1])
  }
  void removeFromACCDevice() const {
#pragma acc exit data delete(this[0:1])
  }
#endif //__PLUMED_USE_OPENACC
};

typedef MatrixTimesMatrix<MatrixProduct> mtimes;
PLUMED_REGISTER_ACTION(mtimes,"MATRIX_PRODUCT")

void MatrixProduct::registerKeywords( Keywords& keys ) {
  keys.setValueDescription("matrix","the product of the two input matrices");
}

void MatrixProduct::calculate( bool noderiv,
                               const MatrixProduct& actdata,
                               const InputVectors& vectors,
                               MatrixElementOutput& output ) {
  std::size_t pp = vectors.arg1.size();
  for(unsigned i=0; i<vectors.nelem; ++i) {
    output.values[0] += vectors.arg1[i]*vectors.arg2[i];
    output.derivs[0][i] = vectors.arg2[i];
    output.derivs[0][pp+i] = vectors.arg1[i];
  }
}

struct Dissimilarities {
  bool squared{false};
  bool periodic{false};
  double min{0};
  double max{0};
  double max_minus_min{0};
  double inv_max_minus_min{0};
  static void registerKeywords( Keywords& keys );
  void setup( MatrixTimesMatrix<Dissimilarities>* action,
              const Value* myval );
  static void calculate( bool noderiv,
                         const Dissimilarities& actdata,
                         const InputVectors& vectors,
                         MatrixElementOutput& output );
#ifdef __PLUMED_USE_OPENACC
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1], squared, periodic, min, max, \
                              max_minus_min, inv_max_minus_min)

  }
  void removeFromACCDevice() const {
#pragma acc exit data delete(inv_max_minus_min, max_minus_min, max, \
                             min, periodic, squared, this[0:1])
  }
#endif //__PLUMED_USE_OPENACC
};

typedef MatrixTimesMatrix<Dissimilarities> dissims;
PLUMED_REGISTER_ACTION(dissims,"DISSIMILARITIES")

void Dissimilarities::registerKeywords( Keywords& keys ) {
  keys.addFlag("SQUARED",false,"calculate the squares of the dissimilarities (this option cannot be used with MATRIX_PRODUCT)");
  keys.setValueDescription("matrix","the dissimilarity matrix");
}

void Dissimilarities::setup( MatrixTimesMatrix<Dissimilarities>* action, const Value* myval ) {
  action->parseFlag("SQUARED",squared);
  if( squared ) {
    action->log.printf("  calculating the squares of the dissimilarities \n");
  }
  periodic = myval->isPeriodic();
  if( periodic ) {
    std::string strmin, strmax;
    myval->getDomain( strmin, strmax );
    Tools::convert( strmin, min );
    Tools::convert( strmax, max );
    max_minus_min = max - min;
    plumed_assert( max_minus_min>0 );
    inv_max_minus_min=1.0/max_minus_min;
  }
}

void Dissimilarities::calculate( bool noderiv, const Dissimilarities& actdata, const InputVectors& vectors, MatrixElementOutput& output ) {
  std::size_t pp = vectors.arg1.size();
  for(unsigned i=0; i<vectors.nelem; ++i) {
    double tmp1 = vectors.arg1[i] - vectors.arg2[i];
    if( actdata.periodic ) {
      tmp1 = tmp1*actdata.inv_max_minus_min;
      tmp1 = Tools::pbc(tmp1);
      tmp1 = tmp1*actdata.max_minus_min;
    }
    output.values[0] += tmp1*tmp1;
    output.derivs[0][i] = 2*tmp1;
    output.derivs[0][pp+i] = -2*tmp1;
  }
  if( actdata.squared ) {
    return ;
  }
  output.values[0] = sqrt( output.values[0] );

  for(unsigned i=0; i<vectors.nelem; ++i) {
    output.derivs[0][i] = output.derivs[0][i] / (2*output.values[0]);
    output.derivs[0][pp+i] = output.derivs[0][vectors.nelem+i] / (2*output.values[0]);
  }
}

}
}
