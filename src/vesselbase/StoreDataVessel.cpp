/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2015 The plumed team
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

void StoreDataVessel::registerKeywords( Keywords& keys ){
  Vessel::registerKeywords(keys); keys.remove("LABEL");
}

StoreDataVessel::StoreDataVessel( const VesselOptions& da ):
Vessel(da),
max_lowmem_stash(3),
vecsize(0),
hard_cut(false),
nspace(0)
{
  ActionWithValue* myval=dynamic_cast<ActionWithValue*>( getAction() );
  if( !myval ) hasderiv=false;
  else hasderiv=!myval->doNotCalculateDerivatives();

  vecsize=getAction()->getNumberOfQuantities();
}

void StoreDataVessel::setHardCutoffOnWeight( const double& mytol ){
  hard_cut=true; wtol=mytol;
}

bool StoreDataVessel::weightCutoffIsOn() const {
  return hard_cut;
}

void StoreDataVessel::resize(){
  plumed_dbg_assert( vecsize>0 );

  if( getAction()->lowmem || !getAction()->derivativesAreRequired() ){
     nspace = 1;
     active_der.resize( max_lowmem_stash * ( 1 + getAction()->getNumberOfDerivatives() ) );
  } else {
     nspace = 1 + getAction()->maxderivatives;
     active_der.resize( getAction()->getFullNumberOfTasks() * ( 1 + getAction()->maxderivatives ) );
  }
  //active_val.resize( getAction()->getFullNumberOfTasks() );
  resizeBuffer( getAction()->getFullNumberOfTasks()*vecsize*nspace );
  local_buffer.resize( getAction()->getFullNumberOfTasks()*vecsize*nspace );
}

// void StoreDataVessel::prepare(){
// // Clear bookeeping arrays  
// //  active_der.assign(active_der.size(),0);
//   //active_val.assign(active_val.size(),0);
// }

void StoreDataVessel::storeValues( const unsigned& myelem, MultiValue& myvals, std::vector<double>& buffer ) const {
  unsigned ibuf = bufstart + myelem * vecsize * nspace; // active_val[myelem]=1;
  for(unsigned icomp=0;icomp<vecsize;++icomp){
      buffer[ibuf] = myvals.get(icomp); ibuf+=nspace;
  }
}

void StoreDataVessel::storeDerivatives( const unsigned& myelem, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const {
  plumed_dbg_assert( getAction()->derivativesAreRequired() && myelem<getAction()->getFullNumberOfTasks() );
  der_list[myelem]=myvals.getNumberActive();

  // Store the values of the components and the derivatives if 
  for(unsigned icomp=0;icomp<vecsize;++icomp){
     unsigned ibuf = bufstart + myelem * ( vecsize*nspace ) + icomp*nspace + 1;
     unsigned kder = getAction()->getFullNumberOfTasks() + myelem * ( nspace - 1 );
     for(unsigned j=0;j<myvals.getNumberActive();++j){
        unsigned jder=myvals.getActiveIndex(j);
        buffer[ibuf] = myvals.getDerivative( icomp, jder ); 
        der_list[kder] = jder; ibuf++; kder++;
     }
  }
}

void StoreDataVessel::retrieveValue( const unsigned& myelem, const bool& normed, std::vector<double>& values ) const {
  plumed_dbg_assert( values.size()==vecsize );
  if( normed && values.size()>2 ){
     unsigned ibuf = myelem * vecsize * nspace;
     values[0]=local_buffer[ibuf]; ibuf+=nspace;
     values[1]=local_buffer[ibuf]; ibuf+=nspace;   // Element 1 contains the norm of the vector
     for(unsigned i=2;i<vecsize;++i){ values[i]=local_buffer[ibuf]/values[1]; ibuf+=nspace; } 
  } else {
     unsigned ibuf = myelem * vecsize * nspace;
     for(unsigned i=0;i<vecsize;++i){ values[i]=local_buffer[ibuf]; ibuf+=nspace; }
  }
}

void StoreDataVessel::retrieveDerivatives( const unsigned& myelem, const bool& normed, MultiValue& myvals ){
  plumed_dbg_assert( myvals.getNumberOfValues()==vecsize && myvals.getNumberOfDerivatives()==getAction()->getNumberOfDerivatives() );

  myvals.clearAll();
  if( getAction()->lowmem ){
      getAction()->performTask( getAction()->getPositionInFullTaskList(myelem), getAction()->getTaskCode(myelem), myvals );
      if( normed ){
          plumed_dbg_assert( myvals.getNumberOfValues()>2 );
          double v = myvals.get(1), weight = 1.0 / v,  wdf = 1.0 / ( v*v*v );
          for(unsigned j=0;j<myvals.getNumberActive();++j){
              double comp2=0.0; unsigned jder=myvals.getActiveIndex(j);
              for(unsigned jcomp=2;jcomp<vecsize;++jcomp) comp2 += myvals.get(jcomp)*myvals.getDerivative( jcomp, jder );
              for(unsigned jcomp=2;jcomp<vecsize;++jcomp) myvals.setDerivative( jcomp, jder, weight*myvals.getDerivative( jcomp, jder ) - wdf*comp2*myvals.get(jcomp) );
          }
      }
  } else {
      // Retrieve the derivatives for elements 0 and 1 - weight and norm
      for(unsigned icomp=0;icomp<2;++icomp){
          unsigned ibuf = myelem * ( vecsize*nspace ) + icomp*nspace + 1;
          unsigned kder = getAction()->getFullNumberOfTasks() + myelem * ( nspace - 1 );
          for(unsigned j=0;j<active_der[myelem];++j){
              myvals.addDerivative( icomp, active_der[kder], local_buffer[ibuf] );
              kder++; ibuf++;
          }
      }

      // Retrieve the derivatives for the vector
      if( vecsize>2 && normed ){
         plumed_dbg_assert( myvals.getNumberOfValues()>2 );
         unsigned kder = getAction()->getFullNumberOfTasks() + myelem * ( nspace - 1 );
         double v = local_buffer[myelem*vecsize*nspace + nspace], weight = 1.0 / v, wdf = 1.0 / ( v*v*v );
         for(unsigned ider=0;ider<active_der[myelem];++ider){
             unsigned ibuf = myelem * vecsize * nspace + 2 * nspace + 1 + ider; double comp2=0.0; 
             for(unsigned jcomp=2;jcomp<vecsize;++jcomp){ comp2 += local_buffer[ibuf-ider-1]*local_buffer[ibuf]; ibuf+=nspace; }  
             ibuf = myelem * vecsize * nspace + 2 * nspace + 1 + ider;
             for(unsigned jcomp=2;jcomp<vecsize;++jcomp){
                 myvals.addDerivative( jcomp, active_der[kder], weight*local_buffer[ibuf] - wdf*comp2*local_buffer[ibuf-ider-1] );
                 ibuf+=nspace;
             }
             kder++;
         }
      } else if( vecsize>2 ){
         for(unsigned icomp=2;icomp<vecsize;++icomp){
             unsigned ibuf = myelem * ( vecsize*nspace ) + icomp*nspace + 1;
             unsigned kder = getAction()->getFullNumberOfTasks() + myelem * ( nspace - 1 );
             for(unsigned j=0;j<active_der[myelem];++j){
                 myvals.addDerivative( icomp, active_der[kder], local_buffer[ibuf] );
                 kder++; ibuf++;
             } 
         } 
      }
      // Now ensure appropriate parts of list are activated
      myvals.emptyActiveMembers();
      unsigned kder = getAction()->getFullNumberOfTasks() + myelem * ( nspace - 1 );
      for(unsigned j=0;j<active_der[myelem];++j){ myvals.putIndexInActiveArray( active_der[kder] ); kder++; }
      myvals.sortActiveList();
  } 
}

bool StoreDataVessel::calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const {

  if( !hard_cut ){
     storeValues( current, myvals, buffer );
     if( !(getAction()->lowmem) && getAction()->derivativesAreRequired() ) storeDerivatives( current, myvals, buffer, der_list );
  } else if( myvals.get(0)>wtol ){
     storeValues( current, myvals, buffer );
     if( !(getAction()->lowmem) && getAction()->derivativesAreRequired() ) storeDerivatives( current, myvals, buffer, der_list );
  } 

  return true;
}

// void StoreDataVessel::buildIndexStores( const unsigned& current, MultiValue& myvals, std::vector<unsigned>& val_index, std::vector<unsigned>& der_index ) const {
//   if( !hard_cut ){
//      val_index[current]=1;
//      if( !(getAction()->lowmem) && getAction()->derivativesAreRequired() ){
//          der_index[current]=myvals.getNumberActive();
//          unsigned kder = getAction()->getFullNumberOfTasks() + current * ( nspace - 1 );
//          for(unsigned j=0;j<myvals.getNumberActive();++j){ der_index[kder] = myvals.getActiveIndex(j); kder++; }
//      }
//   } else if( myvals.get(0)>wtol ){
//      val_index[current]=1;
//      if( !(getAction()->lowmem) && getAction()->derivativesAreRequired() ){
//          der_index[current]=myvals.getNumberActive();
//          unsigned kder = getAction()->getFullNumberOfTasks() + current * ( nspace - 1 );
//          for(unsigned j=0;j<myvals.getNumberActive();++j){ der_index[kder] = myvals.getActiveIndex(j); kder++; }
//      } 
//   }  
// }

void StoreDataVessel::finish( const std::vector<double>& buffer ){
  // Store the buffer locally
  for(unsigned i=0;i<local_buffer.size();++i) local_buffer[i]=buffer[bufstart+i];
}


void StoreDataVessel::setActiveValsAndDerivatives( const std::vector<unsigned>& der_index ){
//   for(unsigned i=0;i<val_index.size();++i){   // Need to think about this bit GAT
//      if( ) active_val[i]    // =val_index[i];
//   }
  if( !getAction()->lowmem && getAction()->derivativesAreRequired() ){
     for(unsigned i=0;i<der_index.size();++i) active_der[i]=der_index[i]; 
  }
}

}
}
