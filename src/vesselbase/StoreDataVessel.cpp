/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2019 The plumed team
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
#include "StoreDataVessel.h"

namespace PLMD {
namespace vesselbase {

void StoreDataVessel::registerKeywords( Keywords& keys ) {
  Vessel::registerKeywords(keys); keys.remove("LABEL");
}

StoreDataVessel::StoreDataVessel( const VesselOptions& da ):
  Vessel(da),
  max_lowmem_stash(3),
  vecsize(0),
  nspace(0)
{
  ActionWithValue* myval=dynamic_cast<ActionWithValue*>( getAction() );
  if( !myval ) hasderiv=false;
  else hasderiv=!myval->doNotCalculateDerivatives();
}

void StoreDataVessel::addActionThatUses( ActionWithVessel* actionThatUses ) {
  userActions.push_back( actionThatUses );
}

void StoreDataVessel::resize() {
  if( getAction()->lowmem || !getAction()->derivativesAreRequired() ) {
    nspace = 1;
    active_der.resize( max_lowmem_stash * ( 1 + getAction()->getNumberOfDerivatives() ) );
  } else {
    if( getAction()->getNumberOfDerivatives()>getAction()->maxderivatives ) {
      error("not enough memory to store derivatives for action " + getAction()->getLabel() + " use LOWMEM option");
    }
    nspace = 1 + getAction()->maxderivatives;
    active_der.resize( getNumberOfStoredValues() * ( 1 + getAction()->maxderivatives ) );
  }
  vecsize=getAction()->getNumberOfQuantities();
  plumed_dbg_assert( vecsize>0 );
  resizeBuffer( getNumberOfStoredValues()*vecsize*nspace );
  local_buffer.resize( getNumberOfStoredValues()*vecsize*nspace );
}

void StoreDataVessel::storeValues( const unsigned& myelem, MultiValue& myvals, std::vector<double>& buffer ) const {
  plumed_dbg_assert( vecsize>0 );
  unsigned jelem = getAction()->getPositionInCurrentTaskList( myelem ); plumed_dbg_assert( jelem<getNumberOfStoredValues() );
  unsigned ibuf = bufstart + jelem * vecsize * nspace;
  for(unsigned icomp=0; icomp<vecsize; ++icomp) {
    buffer[ibuf] += myvals.get(icomp); ibuf+=nspace;
  }
}

void StoreDataVessel::storeDerivatives( const unsigned& myelem, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const {
  plumed_dbg_assert( vecsize>0 && getAction()->derivativesAreRequired() && myelem<getAction()->getFullNumberOfTasks() );
  unsigned jelem = getAction()->getPositionInCurrentTaskList( myelem );

  if( getAction()->getFullNumberOfTasks()==getNumberOfStoredValues() ) {
    der_list[jelem]=myvals.getNumberActive();
    unsigned kder = getNumberOfStoredValues() + jelem * ( nspace - 1 );
    for(unsigned j=0; j<myvals.getNumberActive(); ++j) { der_list[kder] = myvals.getActiveIndex(j); kder++; }
  } else {
    // This ensures that active indices are gathered correctly if
    // we have multiple tasks contributing to a stored quantity
    unsigned kder = getNumberOfStoredValues() + jelem * ( nspace - 1 );
    for(unsigned j=0; j<myvals.getNumberActive(); ++j) {
      bool found=false; unsigned jder = myvals.getActiveIndex(j);
      for(unsigned k=0; k<der_list[jelem]; ++k) {
        if( der_list[kder+k]==jder ) { found=true; break; }
      }
      if(!found) { der_list[kder+der_list[jelem]]=jder; der_list[jelem]++; }
    }
  }

  // Store the values of the components and the derivatives
  for(unsigned icomp=0; icomp<vecsize; ++icomp) {
    unsigned ibuf = bufstart + jelem * ( vecsize*nspace ) + icomp*nspace + 1;
    for(unsigned j=0; j<myvals.getNumberActive(); ++j) {
      unsigned jder=myvals.getActiveIndex(j);
      buffer[ibuf] += myvals.getDerivative( icomp, jder ); ibuf++;
    }
  }
}

void StoreDataVessel::retrieveSequentialValue( const unsigned& jelem, const bool& normed, std::vector<double>& values ) const {
  plumed_dbg_assert( values.size()==vecsize );
  unsigned ibuf = jelem * vecsize * nspace;
  for(unsigned i=0; i<vecsize; ++i) { values[i]=local_buffer[ibuf]; ibuf+=nspace; }
  if( normed && values.size()>2 ) getAction()->normalizeVector( values );
}

void StoreDataVessel::retrieveValueWithIndex( const unsigned& myelem, const bool& normed, std::vector<double>& values ) const {
  plumed_dbg_assert( values.size()==vecsize );
  unsigned jelem = getStoreIndex( myelem );
  retrieveSequentialValue( jelem, normed, values );
}

double StoreDataVessel::retrieveWeightWithIndex( const unsigned& myelem ) const {
  plumed_dbg_assert( vecsize>0 );
  unsigned jelem = getStoreIndex( myelem ); unsigned ibuf = jelem * vecsize * nspace; return local_buffer[ibuf];
}

void StoreDataVessel::retrieveDerivatives( const unsigned& myelem, const bool& normed, MultiValue& myvals ) {
  plumed_dbg_assert( myvals.getNumberOfValues()==vecsize && myvals.getNumberOfDerivatives()==getAction()->getNumberOfDerivatives() );

  myvals.clearAll();
  if( getAction()->lowmem ) {
    recalculateStoredQuantity( myelem, myvals );
    if( normed ) getAction()->normalizeVectorDerivatives( myvals );
  } else {
    unsigned jelem = getAction()->getPositionInCurrentTaskList( myelem );
    // Retrieve the derivatives for elements 0 and 1 - weight and norm
    for(unsigned icomp=0; icomp<vecsize; ++icomp) {
      unsigned ibuf = jelem * ( vecsize*nspace ) + icomp*nspace + 1;
      unsigned kder = getNumberOfStoredValues() + jelem * ( nspace - 1 );
      for(unsigned j=0; j<active_der[jelem]; ++j) {
        myvals.addDerivative( icomp, active_der[kder], local_buffer[ibuf] );
        kder++; ibuf++;
      }
    }
    if( normed ) getAction()->normalizeVectorDerivatives( myvals );
    // Now ensure appropriate parts of list are activated
    myvals.emptyActiveMembers();
    unsigned kder = getNumberOfStoredValues() + jelem * ( nspace - 1 );
    for(unsigned j=0; j<active_der[jelem]; ++j) { myvals.putIndexInActiveArray( active_der[kder] ); kder++; }
    myvals.sortActiveList();
  }
}

void StoreDataVessel::calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const {

  if( myvals.get(0)>epsilon ) {
    storeValues( current, myvals, buffer );
    if( !(getAction()->lowmem) && getAction()->derivativesAreRequired() ) storeDerivatives( current, myvals, buffer, der_list );
  }

  return;
}

void StoreDataVessel::finish( const std::vector<double>& buffer ) {
  // Store the buffer locally
  for(unsigned i=0; i<local_buffer.size(); ++i) local_buffer[i]=buffer[bufstart+i];
}


void StoreDataVessel::setActiveValsAndDerivatives( const std::vector<unsigned>& der_index ) {
  if( !getAction()->lowmem && getAction()->derivativesAreRequired() ) {
    for(unsigned i=0; i<der_index.size(); ++i) active_der[i]=der_index[i];
  }
}

void StoreDataVessel::resizeTemporyMultiValues( const unsigned& nvals ) {
  for(unsigned i=0; i<nvals; ++i) my_tmp_vals.push_back( MultiValue(0,0) );
}

MultiValue& StoreDataVessel::getTemporyMultiValue( const unsigned& ind ) {
  plumed_dbg_assert( ind<my_tmp_vals.size() ); return my_tmp_vals[ind];
}

}
}
