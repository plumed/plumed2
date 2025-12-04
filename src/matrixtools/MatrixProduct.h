
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
#ifndef __PLUMED_matrixtools_MatrixProduct_h
#define __PLUMED_matrixtools_MatrixProduct_h
#include "MatrixTimesMatrix.h"

namespace PLMD {
namespace matrixtools {

struct MatrixProduct {
  static void registerKeywords( Keywords& keys );
  template <typename PTM>
  void setup( MatrixTimesMatrix<MatrixProduct,PTM>* action,
              const Value* myval ) {}
  static void calculate( bool noderiv,
                         const MatrixProduct& actdata,
                         const InputVectors& vectors,
                         MatrixElementOutput& output );
#ifdef __PLUMED_HAS_OPENACC
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1])
  }
  void removeFromACCDevice() const {
#pragma acc exit data delete(this[0:1])
  }
#endif //__PLUMED_HAS_OPENACC
};


void MatrixProduct::registerKeywords( Keywords& keys ) {
  keys.setValueDescription("matrix","the product of the two input matrices");
  keys.addFlag("ELEMENTS_ON_DIAGONAL_ARE_ZERO",false,"set all diagonal elements to zero");
}

void MatrixProduct::calculate(  const bool noderiv,
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
} //namespace matrixtools
} //PLMD
#endif //__PLUMED_matrixtools_MatrixProduct_h
