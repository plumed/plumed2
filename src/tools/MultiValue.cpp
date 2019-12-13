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

namespace PLMD {

MultiValue::MultiValue( const unsigned& nvals, const unsigned& nder ):
  values(nvals),
  nderivatives(nder),
  derivatives(nvals*nder),
  tmpval(0),
  tmpder(nder),
  atLeastOneSet(false)
{
  std::vector<unsigned> myind( nder );
  for(unsigned i=0; i<nder; ++i) myind[i]=i;
  hasDerivatives.createIndexListFromVector( myind );
}

void MultiValue::resize( const unsigned& nvals, const unsigned& nder ) {
  values.resize(nvals); nderivatives=nder; derivatives.resize( nvals*nder );
  tmpder.resize( nder ); hasDerivatives.clear(); std::vector<unsigned> myind( nder );
  for(unsigned i=0; i<nder; ++i) myind[i]=i;
  hasDerivatives.createIndexListFromVector( myind );
  atLeastOneSet=false;
}

void MultiValue::clearAll() {
  if( atLeastOneSet && !hasDerivatives.updateComplete() ) hasDerivatives.updateActiveMembers();
  for(unsigned i=0; i<values.size(); ++i) clear(i);
  clearTemporyDerivatives(); hasDerivatives.deactivateAll(); atLeastOneSet=false;
}

void MultiValue::clear( const unsigned& ival ) {
  values[ival]=0;
  unsigned base=ival*nderivatives, ndert=hasDerivatives.getNumberActive();
  for(unsigned i=0; i<ndert; ++i) derivatives[ base+hasDerivatives[i] ]=0.;
}

void MultiValue::clearTemporyDerivatives() {
  unsigned ndert=hasDerivatives.getNumberActive(); tmpval=0.;
  for(unsigned i=0; i<ndert; ++i) tmpder[ hasDerivatives[i] ]=0.;
}

void MultiValue::chainRule( const unsigned& ival, const unsigned& iout, const unsigned& stride, const unsigned& off,
                            const double& df, const unsigned& bufstart, std::vector<double>& buffer ) {
  if( !hasDerivatives.updateComplete() ) hasDerivatives.updateActiveMembers();

  plumed_dbg_assert( off<stride );
  unsigned base=nderivatives*ival, ndert=hasDerivatives.getNumberActive();
  unsigned start=bufstart+stride*(nderivatives+1)*iout + stride;
  for(unsigned i=0; i<ndert; ++i) {
    unsigned jder=hasDerivatives[i];
    buffer[start+jder*stride] += df*derivatives[base+jder];
  }
}

void MultiValue::copyValues( MultiValue& outvals ) const {
  plumed_dbg_assert( values.size()<=outvals.getNumberOfValues() );
  for(unsigned i=0; i<values.size(); ++i) outvals.setValue( i, values[i] );

}

void MultiValue::copyDerivatives( MultiValue& outvals ) {
  plumed_dbg_assert( values.size()<=outvals.getNumberOfValues() && nderivatives<=outvals.getNumberOfDerivatives() );
  if( !hasDerivatives.updateComplete() ) hasDerivatives.updateActiveMembers();

  outvals.atLeastOneSet=true; unsigned ndert=hasDerivatives.getNumberActive();
  for(unsigned j=0; j<ndert; ++j) {
    unsigned jder=hasDerivatives[j]; outvals.hasDerivatives.activate(jder);
  }

  unsigned base=0, obase=0;
  for(unsigned i=0; i<values.size(); ++i) {
    for(unsigned j=0; j<ndert; ++j) {
      unsigned jder=hasDerivatives[j];
      outvals.derivatives[obase+jder] += derivatives[base+jder];
    }
    obase+=outvals.nderivatives; base+=nderivatives;
  }
}

void MultiValue::quotientRule( const unsigned& nder, const unsigned& oder ) {
  plumed_dbg_assert( nder<values.size() && oder<values.size() );
  if( !hasDerivatives.updateComplete() ) hasDerivatives.updateActiveMembers();

  unsigned ndert=hasDerivatives.getNumberActive(); double wpref;
  unsigned obase=oder*nderivatives, nbase=nder*nderivatives;

  if( fabs(tmpval)>epsilon ) { wpref=1.0/tmpval; }
  else { wpref=1.0; }

  double pref = values[nder]*wpref*wpref;
  for(unsigned j=0; j<ndert; ++j) {
    unsigned jder=hasDerivatives[j];
    derivatives[obase+jder] = wpref*derivatives[nbase+jder]  - pref*tmpder[jder];
  }
  values[oder] = wpref*values[nder];
}

}
