/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2017 The plumed team
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

namespace PLMD {

MultiValue::MultiValue( const unsigned& nvals, const unsigned& nder ):
  task_index(0),
  values(nvals),
  nderivatives(nder),
  derivatives(nvals*nder),
  hasderiv(nvals*nder,false),
  nactive(nvals),
  active_list(nvals*nder),
  atLeastOneSet(false)
{
}

void MultiValue::resize( const unsigned& nvals, const unsigned& nder ) {
  values.resize(nvals); nderivatives=nder; derivatives.resize( nvals*nder );
  hasderiv.resize(nvals*nder,false); nactive.resize(nvals); active_list.resize(nvals*nder); 
  atLeastOneSet=false;
}

void MultiValue::clearAll() {
  if( !atLeastOneSet ) return;
  for(unsigned i=0; i<values.size(); ++i) clear(i);
  atLeastOneSet=false;
}

void MultiValue::clear( const unsigned& ival ) {
  values[ival]=0; 
  unsigned base=ival*nderivatives;
  for(unsigned i=0; i<nactive[ival]; ++i){
     unsigned k = base+active_list[base+i]; derivatives[k]=0.; hasderiv[k]=false; 
  }
  nactive[ival]=0;
#ifndef NDEBUG
  for(unsigned i=0; i<nderivatives;++i) plumed_dbg_assert( hasderiv[base+i]==false ); 
#endif
}

}
