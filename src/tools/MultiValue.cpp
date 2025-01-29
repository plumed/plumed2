/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2023 The plumed team
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
#include "MultiValue.h"
#include "Tools.h"

namespace PLMD {

MultiValue::MultiValue( const size_t& nvals, const size_t& nder, const size_t& nmat, const size_t& maxcol, const size_t& nbook ):
  task_index(0),
  task2_index(0),
  values(nvals),
  nderivatives(nder),
  derivatives(nvals*nder),
  hasderiv(nvals*nder,false),
  tmpval(0),
  nactive(nvals),
  active_list(nvals*nder),
  tmpder(nder),
  atLeastOneSet(false),
  vector_call(false),
  nindices(0),
  nsplit(0),
  nmatrix_cols(maxcol),
  matrix_row_stash(nmat*maxcol,0),
  matrix_force_stash(nder*nmat),
  matrix_bookeeping(nbook,0),
  matrix_row_nderivatives(nmat,0),
  matrix_row_derivative_indices(nmat) {
  for(unsigned i=0; i<nmat; ++i) {
    matrix_row_derivative_indices[i].resize( nder );
  }
  // This is crap that will be deleted in future
  std::vector<unsigned> myind( nder );
  for(unsigned i=0; i<nder; ++i) {
    myind[i]=i;
  }
}

void MultiValue::resize( const size_t& nvals, const size_t& nder, const size_t& nmat, const size_t& maxcol, const size_t& nbook ) {
  values.resize(nvals);
  nderivatives=nder;
  derivatives.resize( nvals*nder );
  hasderiv.resize(nvals*nder,false);
  nactive.resize(nvals);
  active_list.resize(nvals*nder);
  nmatrix_cols=maxcol;
  matrix_row_stash.resize(nmat*maxcol,0);
  matrix_force_stash.resize(nmat*nder,0);
  matrix_bookeeping.resize(nbook, 0);
  matrix_row_nderivatives.resize(nmat,0);
  matrix_row_derivative_indices.resize(nmat);
  atLeastOneSet=false;
  for(unsigned i=0; i<nmat; ++i) {
    matrix_row_derivative_indices[i].resize( nder );
  }
  // All crap from here onwards
  tmpder.resize( nder );
  std::vector<unsigned> myind( nder );
  for(unsigned i=0; i<nder; ++i) {
    myind[i]=i;
  }
}

void MultiValue::clearAll() {
  for(unsigned i=0; i<values.size(); ++i) {
    values[i]=0;
  }
  // Clear matrix row
  std::fill( matrix_row_stash.begin(), matrix_row_stash.end(), 0 );
  // Clear matrix derivative indices
  std::fill( matrix_row_nderivatives.begin(), matrix_row_nderivatives.end(), 0 );
  // Clear matrix forces
  std::fill(matrix_force_stash.begin(),matrix_force_stash.end(),0);
  if( !atLeastOneSet ) {
    return;
  }
  for(unsigned i=0; i<values.size(); ++i) {
    clearDerivatives(i);
  }
  atLeastOneSet=false;
}

void MultiValue::clearDerivatives( const unsigned& ival ) {
  values[ival]=0;
  if( !atLeastOneSet ) {
    return;
  }
  unsigned base=ival*nderivatives;
  for(unsigned i=0; i<nactive[ival]; ++i) {
    unsigned k = base+active_list[base+i];
    derivatives[k]=0.;
    hasderiv[k]=false;
  }
  nactive[ival]=0;
#ifndef NDEBUG
  for(unsigned i=0; i<nderivatives; ++i) {
    if( hasderiv[base+i] ) {
      std::string num1, num2;
      Tools::convert(ival,num1);
      Tools::convert(i,num2);
      plumed_merror("FAILING TO CLEAR VALUE " + num1 + " DERIVATIVE " + num2 + " IS PROBLEMATIC");
    }
    // As this is debugging code we hard code an escape here otherwise this check is really expensive
    if( i>1000 ) {
      return;
    }
  }
#endif
}

}
