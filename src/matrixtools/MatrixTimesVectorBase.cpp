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
#include "MatrixTimesVectorBase.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace matrixtools {

class MatrixTimesVectorRowSums {
public:
  static void performTask( const MatrixTimesVectorInput& input, MatrixTimesVectorOutput& output );
  static void gatherVectorForces( const MatrixTimesVectorForceInput& fin, View<double,helpers::dynamic_extent>& fout ) {}
};

typedef MatrixTimesVectorBase<MatrixTimesVectorRowSums> mycr;
PLUMED_REGISTER_ACTION(mycr,"MATRIX_VECTOR_PRODUCT_ROWSUMS")

void MatrixTimesVectorRowSums::performTask( const MatrixTimesVectorInput& input, MatrixTimesVectorOutput& output ) {
  for(unsigned i=0; i<input.rowlen; ++i) {
    output.values[0] += input.matrow[i];
    if( input.noderiv ) {
      continue;
    }
    output.matrow_deriv[i] = 1;
  }
}

class MatrixTimesVectorProper {
public:
  static void performTask( const MatrixTimesVectorInput& input, MatrixTimesVectorOutput& output );
  static void gatherVectorForces( const MatrixTimesVectorForceInput& fin, View<double,helpers::dynamic_extent>& fout );
};

typedef MatrixTimesVectorBase<MatrixTimesVectorProper> mycp;
PLUMED_REGISTER_ACTION(mycp,"MATRIX_VECTOR_PRODUCT_PROPER")

void MatrixTimesVectorProper::performTask( const MatrixTimesVectorInput& input, MatrixTimesVectorOutput& output ) {
  for(unsigned i=0; i<input.rowlen; ++i) {
    std::size_t ind = input.indices[i];
    output.values[0] += input.matrow[i]*input.vector[ind];
    if( input.noderiv ) {
      continue;
    }
    output.matrow_deriv[i] = input.vector[ind];
    output.vector_deriv[i] = input.matrow[i];
  }
}

void MatrixTimesVectorProper::gatherVectorForces( const MatrixTimesVectorForceInput& fin, View<double,helpers::dynamic_extent>& fout ) {
  for(unsigned i=0; i<fin.rowlen; ++i) {
    fout[fin.indices[i]] += fin.force*fin.deriv[i];
  }
}

}
}
