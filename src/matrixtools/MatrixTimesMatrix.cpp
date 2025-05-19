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
#include "MatrixTimesMatrix.h"
#include "core/ActionRegister.h"

//+PLUMEDOC MCOLVAR MATRIX_PRODUCT
/*
Calculate the product of two matrices

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC ANALYSIS DISSIMILARITIES
/*
Calculate the matrix of dissimilarities between a trajectory of atomic configurations.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace matrixtools {

class MatrixProduct {
public:
  static void registerKeywords( Keywords& keys );
  void setup( MatrixTimesMatrix<MatrixProduct>* action, const Value* myval ) {}
  static void calculate( bool noderiv, const MatrixProduct& actdata, const InputVectors& vectors, MatrixElementOutput& output );
};

typedef MatrixTimesMatrix<MatrixProduct> mtimes;
PLUMED_REGISTER_ACTION(mtimes,"MATRIX_PRODUCT")

void MatrixProduct::registerKeywords( Keywords& keys ) {
  keys.setValueDescription("matrix","the product of the two input matrices");
}

void MatrixProduct::calculate( bool noderiv, const MatrixProduct& actdata, const InputVectors& vectors, MatrixElementOutput& output ) {
  std::size_t pp = vectors.arg1.size();
  for(unsigned i=0; i<vectors.nelem; ++i) {
    output.values[0] += vectors.arg1[i]*vectors.arg2[i];
    output.derivs[0][i] = vectors.arg2[i];
    output.derivs[0][pp+i] = vectors.arg1[i];
  }
}

class Dissimilarities {
public:
  bool squared{false};
  bool periodic{false};
  double min{0}, max{0}, max_minus_min{0}, inv_max_minus_min{0};
  static void registerKeywords( Keywords& keys );
  void setup( MatrixTimesMatrix<Dissimilarities>* action, const Value* myval );
  static void calculate( bool noderiv, const Dissimilarities& actdata, const InputVectors& vectors, MatrixElementOutput& output );
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
