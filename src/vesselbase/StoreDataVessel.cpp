/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
  Vessel::registerKeywords(keys);
}

StoreDataVessel::StoreDataVessel( const VesselOptions& da ):
Vessel(da),
max_lowmem_stash(3),
vecsize(0)
{
}

void StoreDataVessel::completeSetup( const unsigned& dstart, const unsigned& nvec ){
  if( (dstart+nvec)>getAction()->getNumberOfQuantities() ) error("this vessel can not be used with this action");
  data_start=dstart; vecsize = nvec; fvec.resize( vecsize ); 
}

void StoreDataVessel::resize(){
  plumed_dbg_assert( vecsize>0 );
  final_derivatives.resize( getAction()->getNumberOfDerivatives() );
  if( getAction()->lowmem ){
     nspace = 1;
     active_der.resize( max_lowmem_stash * ( 1 + getAction()->getNumberOfDerivatives() ) );
     local_derivatives.resize( max_lowmem_stash * vecsize * getAction()->getNumberOfDerivatives() );
  } else {
     nspace = 1 + getAction()->maxderivatives;
     active_der.resize( getAction()->getNumberOfTasks() * ( 1 + getAction()->maxderivatives ) );
  }
  resizeBuffer( getAction()->getNumberOfTasks()*vecsize*nspace );
}

void StoreDataVessel::prepare(){
// Clear bookeeping arrays  
  active_der.assign(active_der.size(),0);
}

void StoreDataVessel::getIndexList( const unsigned& ntotal, const unsigned& jstore, const unsigned& maxder, std::vector<unsigned>& indices ){
  getAction()->getIndexList( ntotal, jstore, maxder, indices );
}

void StoreDataVessel::recompute( const unsigned& ivec, const unsigned& jstore ){
  plumed_dbg_assert( getAction()->lowmem && jstore<max_lowmem_stash );
  // Get the underlying action with value
  ActionWithVessel* act = getAction();

  // Set the task we want to reperform
  act->current = ivec;
  // Reperform the task
  performTask( jstore );
  // Store the derivatives
  storeDerivativesLowMem( jstore );
 
  // Clear up afterwards
  getAction()->clearAfterTask();
}

void StoreDataVessel::storeValues( const unsigned& myelem ){
  ActionWithVessel* act = getAction();
  unsigned ibuf = myelem * vecsize * nspace;
  for(unsigned icomp=data_start;icomp<data_start + vecsize;++icomp){
     setBufferElement( ibuf, act->getElementValue( icomp ) ); ibuf+=nspace;  
  }
}

void StoreDataVessel::storeDerivativesHighMem( const unsigned& myelem ){
  ActionWithVessel* act=getAction();
  getIndexList( act->getNumberOfTasks(), myelem, nspace-1, active_der );

  // Store the values of the components and the derivatives if 
  unsigned nder = act->getNumberOfDerivatives();
  for(unsigned icomp=data_start;icomp<data_start + vecsize;++icomp){
     unsigned ibuf = myelem * ( vecsize*nspace ) + (icomp-data_start)*nspace + 1;
     unsigned kder = act->getNumberOfTasks() + myelem * ( nspace - 1 );
     for(unsigned jder=0;jder<active_der[myelem];++jder){
        addToBufferElement( ibuf, act->getElementDerivative(nder*icomp + active_der[ kder ]) );
        kder++; ibuf++;
     }
  }
}

void StoreDataVessel::storeDerivativesLowMem( const unsigned& jstore ){
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
  unsigned myelem = getAction()->current;
  // Normalize vector if it is required
  finishTask( myelem );

  // Store the values 
  storeValues( myelem );

  // Store the derivatives if we are not using low memory
  if( !(getAction()->lowmem) ) storeDerivativesHighMem( myelem );  

  return true;
}

void StoreDataVessel::finish(){
  if(!getAction()->lowmem) comm.Sum( &active_der[0], active_der.size() );
}

double StoreDataVessel::chainRule( const unsigned& ival, const unsigned& ider, const std::vector<double>& df ){
  plumed_dbg_assert( df.size()==vecsize );
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
     plumed_dbg_assert( ival<getAction()->getNumberOfTasks() );
     unsigned ibuf=ival*(vecsize*nspace) + 1 + ider;
     for(unsigned icomp=0;icomp<vecsize;++icomp){
         dfout+=df[icomp]*getBufferElement(ibuf);
         ibuf+=nspace;
     }
  }
  return dfout;
}

void StoreDataVessel::chainRule( const unsigned& ival, const std::vector<double>& df ){
  plumed_dbg_assert( df.size()==vecsize );
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
     plumed_dbg_assert( ival<getAction()->getNumberOfTasks() );
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

void StoreDataVessel::chainRule( const unsigned& ival, const unsigned& jout, const std::vector<double>& df, ActionWithVessel* act ){
  plumed_dbg_assert( act->getNumberOfDerivatives()==final_derivatives.size() && jout<act->getNumberOfQuantities() );
  plumed_dbg_assert( act->getNumberOfDerivatives()==final_derivatives.size() && jout<act->getNumberOfQuantities() );
  chainRule( ival, df );

  unsigned kder, outstart = jout*getAction()->getNumberOfDerivatives();
  if( getAction()->lowmem ) kder = max_lowmem_stash + ival*getAction()->getNumberOfDerivatives();
  else kder = getAction()->getNumberOfTasks() + ival*(nspace-1);

  for(unsigned i=0;i<active_der[ival];++i){
     act->addElementDerivative( outstart + active_der[kder+i], final_derivatives[i] );
  }
  act->activateIndexes( kder, active_der[ival], active_der );
}

void StoreDataVessel::chainRule( const unsigned& ival, const unsigned& jout, const unsigned& ider, const std::vector<double>& df, ActionWithVessel* act ){
  plumed_dbg_assert( act->getNumberOfDerivatives()==final_derivatives.size() && jout<act->getNumberOfQuantities() );
  plumed_dbg_assert( act->getNumberOfDerivatives()==final_derivatives.size() && jout<act->getNumberOfQuantities() );

  unsigned kder, jder;
  if( getAction()->lowmem ) kder = max_lowmem_stash + ival*getAction()->getNumberOfDerivatives();
  else kder = getAction()->getNumberOfTasks() + ival*(nspace-1);

  jder = jout*act->getNumberOfDerivatives() + active_der[kder+ider];
  act->addElementDerivative( jder, chainRule( ival, ider, df ) );
  act->activateIndex( active_der[kder+ider] );
}

void StoreDataVessel::chainRule( const unsigned& ival, const std::vector<double>& df, Value* val ){
  plumed_dbg_assert( val->getNumberOfDerivatives()==final_derivatives.size() );
  chainRule( ival, df ); 

  unsigned kder;
  if( getAction()->lowmem ) kder = max_lowmem_stash + ival*getAction()->getNumberOfDerivatives();
  else kder = getAction()->getNumberOfTasks() + ival*(nspace-1);

  for(unsigned i=0;i<active_der[ival];++i) val->addDerivative( active_der[kder+i] , final_derivatives[i] );
}

void StoreDataVessel::chainRule( const unsigned& ival, const std::vector<double>& df, std::vector<double>& derout ){
  plumed_dbg_assert( derout.size()==final_derivatives.size() );
  chainRule( ival, df );

  unsigned kder;
  if( getAction()->lowmem ) kder = max_lowmem_stash + ival*getAction()->getNumberOfDerivatives();
  else kder = getAction()->getNumberOfTasks() + ival*(nspace-1);

  for(unsigned i=0;i<active_der[ival];++i) derout[ active_der[kder+i] ] += final_derivatives[i];
}

void StoreDataVessel::transformComponents( const unsigned& jstore, const double& weight, double& wdf, const std::vector<double>& dfvec ){
  plumed_dbg_assert( dfvec.size()==vecsize );
  ActionWithVessel* act = getAction();
  unsigned myelem = act->current;

  unsigned ibuf = myelem * vecsize * nspace;
  for(unsigned icomp=0;icomp<vecsize;++icomp){
     fvec[icomp]=getBufferElement(ibuf);
     setBufferElement( ibuf, weight*getBufferElement(ibuf)  ); ibuf+=nspace;
  }  

  if( !act->lowmem ) {
      for(unsigned ider=0;ider<active_der[myelem];++ider){
          double comp2=0.0; unsigned ibuf = myelem * vecsize * nspace + 1 + ider;
          for(unsigned jcomp=0;jcomp<vecsize;++jcomp){
              comp2  += fvec[jcomp]*dfvec[jcomp]*getBufferElement(ibuf);
              ibuf += nspace;
          }
          ibuf = myelem * vecsize * nspace + 1 + ider;
          for(unsigned jcomp=0;jcomp<vecsize;++jcomp){
             setBufferElement( ibuf, weight*dfvec[jcomp]*getBufferElement(ibuf) + wdf*comp2*fvec[jcomp] );
             ibuf += nspace;
          }
      }
  } else {
      plumed_dbg_assert( jstore<max_lowmem_stash );
      unsigned maxder = act->getNumberOfDerivatives();
      for(unsigned ider=0;ider<active_der[jstore];++ider){
          double comp2=0.0; 
          unsigned ibuf = jstore * vecsize * maxder + ider;
          for(unsigned jcomp=0;jcomp<vecsize;++jcomp){
              comp2 += fvec[jcomp]*dfvec[jcomp]*local_derivatives[ibuf];
              ibuf += maxder; 
          }
          ibuf = jstore * vecsize * maxder + ider; 
          for(unsigned jcomp=0;jcomp<vecsize;++jcomp){
             local_derivatives[ibuf] = weight*dfvec[jcomp]*local_derivatives[ibuf] + wdf*comp2*fvec[jcomp]; 
             ibuf += maxder;
          }
      }
  } 
}

void StoreDataVessel::chainRuleForComponent( const unsigned& ival, const unsigned& jcomp, const unsigned& jout, const double& df, ActionWithVessel* act ){ 
  unsigned outstart=jout*act->getNumberOfDerivatives();

  if(getAction()->lowmem){
     plumed_dbg_assert( ival<max_lowmem_stash );
     unsigned maxder = getAction()->getNumberOfDerivatives();
     unsigned ibuf = ival*(vecsize*maxder) + jcomp*maxder; 
     unsigned kder = max_lowmem_stash + ival*maxder;
     for(unsigned ider=0;ider<active_der[ival];++ider){
        act->addElementDerivative( outstart + active_der[kder+ider], df*local_derivatives[ibuf+ider] ); 
     }
     act->activateIndexes( kder, active_der[ival], active_der );
  } else {
     plumed_dbg_assert( ival<getAction()->getNumberOfTasks() );
     unsigned ibuf = ival*(vecsize*nspace) + jcomp*nspace + 1;
     unsigned kder = getAction()->getNumberOfTasks() + ival*(nspace-1);
     for(unsigned ider=0;ider<active_der[ival];++ider){
        act->addElementDerivative( outstart + active_der[kder+ider], df*getBufferElement(ibuf+ider) );  
     }
     act->activateIndexes( kder, active_der[ival], active_der );
  }
}

}
}
