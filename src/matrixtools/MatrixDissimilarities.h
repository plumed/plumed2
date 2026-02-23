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
#ifndef __PLUMED_matrixtools_MatrixDissimilarities_h
#define __PLUMED_matrixtools_MatrixDissimilarities_h
#include "MatrixTimesMatrix.h"

namespace PLMD {
namespace matrixtools {
struct Dissimilarities {
  typedef void isDissimilarities;
  bool squared{false};
  bool periodic{false};
  double min{0};
  double max{0};
  double max_minus_min{0};
  double inv_max_minus_min{0};
  static void registerKeywords( Keywords& keys );
  template <typename PTM>
  void setup( MatrixTimesMatrix<Dissimilarities,PTM>* action,
              const Value* myval );
  static void calculate( bool noderiv,
                         const Dissimilarities& actdata,
                         const InputVectors& vectors,
                         MatrixElementOutput& output );
#ifdef __PLUMED_HAS_OPENACC
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1], squared, periodic, min, max, \
                              max_minus_min, inv_max_minus_min)

  }
  void removeFromACCDevice() const {
#pragma acc exit data delete(inv_max_minus_min, max_minus_min, max, \
                             min, periodic, squared, this[0:1])
  }
#endif //__PLUMED_HAS_OPENACC
};


void Dissimilarities::registerKeywords( Keywords& keys ) {
  keys.setDisplayName("DISSIMILARITY");
  keys.addFlag("SQUARED",false,"calculate the squares of the dissimilarities (this option cannot be used with MATRIX_PRODUCT)");
  keys.setValueDescription("matrix","the dissimilarity matrix");
}

template <typename PTM>
void Dissimilarities::setup( MatrixTimesMatrix<Dissimilarities, PTM>* action,
                             const Value* myval ) {
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

void Dissimilarities::calculate( const bool noderiv,
                                 const Dissimilarities& actdata,
                                 const InputVectors& vectors,
                                 MatrixElementOutput& output ) {
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
} //namespace matrixtools
} //PLMD
#endif //__PLUMED_matrixtools_MatrixDissimilarities_h
