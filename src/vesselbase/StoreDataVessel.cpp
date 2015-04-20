/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013,2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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

void StoreDataVessel::registerKeywords( Keywords& keys ){
  Vessel::registerKeywords(keys); keys.remove("LABEL");
}

StoreDataVessel::StoreDataVessel( const VesselOptions& da ):
Vessel(da),
max_lowmem_stash(3),
vecsize(0),
hard_cut(false)
{
  ActionWithValue* myval=dynamic_cast<ActionWithValue*>( getAction() );
  if( !myval ) hasderiv=false;
  else hasderiv=!myval->doNotCalculateDerivatives();
}

void StoreDataVessel::setHardCutoffOnWeight( const double& mytol ){
  hard_cut=true; wtol=mytol;
}

bool StoreDataVessel::weightCutoffIsOn() const {
  return hard_cut;
}

void StoreDataVessel::completeSetup( const unsigned& dstart, const unsigned& nvec ){
  if( (dstart+nvec)>getAction()->getNumberOfQuantities() ) error("this vessel can not be used with this action");
  data_start=dstart; vecsize = nvec; fvec.resize( vecsize ); 
}

void StoreDataVessel::resize(){
  plumed_dbg_assert( vecsize>0 );
  if( getAction()->derivativesAreRequired() ) final_derivatives.resize( getAction()->getNumberOfDerivatives() );

  if( getAction()->lowmem || !getAction()->derivativesAreRequired() ){
     nspace = 1;
     active_der.resize( max_lowmem_stash * ( 1 + getAction()->getNumberOfDerivatives() ) );
     local_derivatives.resize( max_lowmem_stash * vecsize * getAction()->getNumberOfDerivatives() );
  } else {
     nspace = 1 + getAction()->maxderivatives;
     active_der.resize( getAction()->getFullNumberOfTasks() * ( 1 + getAction()->maxderivatives ) );
  }
  active_val.resize( getAction()->getFullNumberOfTasks() );
  resizeBuffer( getAction()->getFullNumberOfTasks()*vecsize*nspace );
}

void StoreDataVessel::prepare(){
// Clear bookeeping arrays  
  active_der.assign(active_der.size(),0);
  active_val.assign(active_val.size(),0);
}

void StoreDataVessel::getIndexList( const unsigned& ntotal, const unsigned& jstore, const unsigned& maxder, std::vector<unsigned>& indices ){
  getAction()->getIndexList( ntotal, jstore, maxder, indices );
}

void StoreDataVessel::setTaskToRecompute( const unsigned& ivec ){
 getAction()->current = getAction()->fullTaskList[ivec];
 getAction()->task_index = ivec;
}

void StoreDataVessel::recompute( const unsigned& ivec, const unsigned& jstore ){
  plumed_dbg_assert( getAction()->derivativesAreRequired() && getAction()->lowmem && jstore<max_lowmem_stash );
  // Set the task we want to reperform
  setTaskToRecompute( ivec );
  // Reperform the task
  performTask( jstore );
  // Store the derivatives
  storeDerivativesLowMem( jstore );
 
  // Clear up afterwards
  getAction()->clearAfterTask();
}

void StoreDataVessel::storeValues( const unsigned& myelem ){
  ActionWithVessel* act = getAction(); active_val[myelem]=1; // Keeps track of which values are stashed
  unsigned ibuf = myelem * vecsize * nspace;
  for(unsigned icomp=data_start;icomp<data_start + vecsize;++icomp){
     setBufferElement( ibuf, act->getElementValue( icomp ) ); ibuf+=nspace;  
  }
}

void StoreDataVessel::storeDerivativesHighMem( const unsigned& myelem ){
  plumed_dbg_assert( getAction()->derivativesAreRequired() && myelem<getAction()->getFullNumberOfTasks() );
  ActionWithVessel* act=getAction();
  getIndexList( act->getFullNumberOfTasks(), myelem, nspace-1, active_der );

  // Store the values of the components and the derivatives if 
  unsigned nder = act->getNumberOfDerivatives();
  for(unsigned icomp=data_start;icomp<data_start + vecsize;++icomp){
     unsigned ibuf = myelem * ( vecsize*nspace ) + (icomp-data_start)*nspace + 1;
     unsigned kder = act->getFullNumberOfTasks() + myelem * ( nspace - 1 );
     for(unsigned jder=0;jder<active_der[myelem];++jder){
        addToBufferElement( ibuf, act->getElementDerivative(nder*icomp + active_der[ kder ]) );
        kder++; ibuf++;
     }
  }
}

void StoreDataVessel::storeDerivativesLowMem( const unsigned& jstore ){
  plumed_dbg_assert( getAction()->derivativesAreRequired() );
  // Store the indexes that have derivatives
  ActionWithVessel* act=getAction();
  unsigned nder = act->getNumberOfDerivatives();
  getIndexList( max_lowmem_stash, jstore, nder, active_der );

  // Stash the derivatives
  for(unsigned icomp=data_start;icomp<data_start + vecsize;++icomp){
     unsigned ibuf = jstore * vecsize * nder + (icomp-data_start)*nder;
     unsigned kder = max_lowmem_stash + jstore*nder;
     for(unsigned jder=0;jder<active_der[jstore];++jder){
        local_derivatives[ibuf] = act->getElementDerivative( nder*icomp + active_der[ kder ] );
        ibuf++; kder++;
     }
  }
}

bool StoreDataVessel::calculate(){
  unsigned myelem = getAction()->getCurrentPositionInTaskList();
  // Normalize vector if it is required
  finishTask( myelem );

  // Store the values 
  if( !hard_cut ){
     storeValues( myelem );
     // Store the derivatives if we are not using low memory
     if( !(getAction()->lowmem) && getAction()->derivativesAreRequired() ) storeDerivativesHighMem( myelem );  
  } else {
     if( getAction()->getElementValue(getAction()->getIndexOfWeight())>wtol ){ 
         storeValues( myelem );
         // Store the derivatives if we are not using low memory
         if( !(getAction()->lowmem) && getAction()->derivativesAreRequired() ) storeDerivativesHighMem( myelem );
     }
  }

  return true;
}

void StoreDataVessel::finish(){
  // Update what variables we need
  comm.Sum( &active_val[0], active_val.size() );

  // Get the active derivatives
  if(!getAction()->lowmem && getAction()->derivativesAreRequired() ) comm.Sum( &active_der[0], active_der.size() );
}

double StoreDataVessel::chainRule( const unsigned& ival, const unsigned& ider, const std::vector<double>& df ){
  plumed_dbg_assert( getAction()->derivativesAreRequired() && df.size()==vecsize );
  // Clear final derivatives array
  final_derivatives.assign( final_derivatives.size(), 0.0 );

  double dfout = 0.0;
  if(getAction()->lowmem){
     plumed_dbg_assert( ival<max_lowmem_stash );
     unsigned maxder = getAction()->getNumberOfDerivatives();
     unsigned ibuf=ival*(vecsize*maxder) + ider;
     for(unsigned icomp=0;icomp<vecsize;++icomp){
         dfout+=df[icomp]*local_derivatives[ibuf];
         ibuf+=maxder;
     }  
  } else {
     plumed_dbg_assert( ival<getAction()->getFullNumberOfTasks() );
     unsigned ibuf=ival*(vecsize*nspace) + 1 + ider;
     for(unsigned icomp=0;icomp<vecsize;++icomp){
         dfout+=df[icomp]*getBufferElement(ibuf);
         ibuf+=nspace;
     }
  }
  return dfout;
}

void StoreDataVessel::chainRule( const unsigned& ival, const std::vector<double>& df ){
  plumed_dbg_assert( getAction()->derivativesAreRequired() && df.size()==vecsize );
  // Clear final derivatives array
  final_derivatives.assign( final_derivatives.size(), 0.0 );

  if(getAction()->lowmem){
     plumed_dbg_assert( ival<max_lowmem_stash );
     unsigned maxder = getAction()->getNumberOfDerivatives();
     for(unsigned ider=0;ider<active_der[ival];++ider){
        final_derivatives[ider]=0.0;
        unsigned ibuf=ival*(vecsize*maxder) + ider; 
        for(unsigned jcomp=0;jcomp<vecsize;++jcomp){
            final_derivatives[ider]+=df[jcomp]*local_derivatives[ibuf]; 
            ibuf+=maxder;
        }
     }
  } else {
     plumed_dbg_assert( ival<getAction()->getFullNumberOfTasks() );
     for(unsigned ider=0;ider<active_der[ival];++ider){
         final_derivatives[ider]=0.0;
         unsigned ibuf=ival*(vecsize*nspace) + 1 + ider; 
         for(unsigned jcomp=0;jcomp<vecsize;++jcomp){
             final_derivatives[ider]+=df[jcomp]*getBufferElement(ibuf);
             ibuf+=nspace;
         }
     }
  }
}

void StoreDataVessel::chainRule( const unsigned& ival, const std::vector<double>& df, Value* val ){
  plumed_dbg_assert( getAction()->derivativesAreRequired() && val->getNumberOfDerivatives()==final_derivatives.size() );
  chainRule( ival, df );

  unsigned kder;
  if( getAction()->lowmem ) kder = max_lowmem_stash + ival*getAction()->getNumberOfDerivatives();
  else kder = getAction()->getFullNumberOfTasks() + ival*(nspace-1);

  for(unsigned i=0;i<active_der[ival];++i) val->addDerivative( active_der[kder+i] , final_derivatives[i] );
}

}
}
