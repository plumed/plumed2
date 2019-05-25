/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2019 The plumed team
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

MultiValue::MultiValue( const unsigned& nvals, const unsigned& nder, const unsigned ncols, const unsigned nmat, const unsigned nfder ):
  task_index(0),
  values(nvals),
  nderivatives(nder),
  derivatives(nvals*nder),
  hasderiv(nvals*nder,false),
  nactive(nvals),
  active_list(nvals*nder),
  atLeastOneSet(false),
  nmatrix_cols(ncols),
  nmat_force(nfder),
  rerunning_matrix(false),
  matrix_element_nind(nmat),
  matrix_element_indices(ncols*nmat),
  matrix_element_stash(ncols*nmat),
  matrix_force_stash(nfder*nmat),
  vector_call(false),
  nindices(0),
  nfblock(0),
  nsplit(0),
  mat_nindices(nmat,0),
  mat_indices(nmat),
  tmp_atoms(2),
  symfunc_index(0),
  symfunc_tmp_derivs(nvals)
{
  for(unsigned i=0; i<nmat; ++i) mat_indices[i].resize( nder );
}

void MultiValue::resize( const unsigned& nvals, const unsigned& nder, const unsigned& ncols, const unsigned& nmat ) {
  values.resize(nvals); nderivatives=nder; derivatives.resize( nvals*nder );
  hasderiv.resize(nvals*nder,false); nactive.resize(nvals); active_list.resize(nvals*nder);
  nmatrix_cols=ncols; matrix_element_nind.resize(nmat); matrix_element_indices.resize(ncols*nmat);
  matrix_element_stash.resize(ncols*nmat); nindices=0; mat_nindices.resize(nmat,0); mat_indices.resize(nmat);
  for(unsigned i=0; i<nmat; ++i) mat_indices[i].resize( nder );
  symfunc_tmp_derivs.resize( nvals ); atLeastOneSet=false;
}

void MultiValue::clearAll() {
  for(unsigned i=0; i<values.size(); ++i) values[i]=0;
  // Clear matrix indices
  for(unsigned i=0; i<mat_nindices.size(); ++i) mat_nindices[i]=0;
  for(unsigned i=0; i<matrix_element_nind.size(); ++i) matrix_element_nind[i]=0;
  if( !atLeastOneSet ) return;
  for(unsigned i=0; i<values.size(); ++i) clear(i);
  atLeastOneSet=false;
}

void MultiValue::clear( const unsigned& ival ) {
  values[ival]=0;
  if( !atLeastOneSet ) return;
  unsigned base=ival*nderivatives;
  for(unsigned i=0; i<nactive[ival]; ++i) {
    unsigned k = base+active_list[base+i]; derivatives[k]=0.; hasderiv[k]=false;
  }
  nactive[ival]=0;
#ifndef NDEBUG
  for(unsigned i=0; i<nderivatives; ++i) {
    if( hasderiv[base+i] ) {
      std::string num1, num2;
      Tools::convert(ival,num1); Tools::convert(i,num2);
      plumed_merror("FAILING TO CLEAR VALUE " + num1 + " DERIVATIVE " + num2 + " IS PROBLEMATIC");
    }
  }
#endif
}

}
