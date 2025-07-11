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
