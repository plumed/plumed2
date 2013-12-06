/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
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
#include "vesselbase/VesselRegister.h"
#include "vesselbase/ActionWithVessel.h"
#include "multicolvar/MultiColvarFunction.h"
#include "StoreVectorsVessel.h"
#include "VectorMultiColvar.h"

namespace PLMD {
namespace crystallization{

void StoreVectorsVessel::registerKeywords( Keywords& keys ){
  vesselbase::StoreDataVessel::registerKeywords(keys);
}

StoreVectorsVessel::StoreVectorsVessel( const vesselbase::VesselOptions& da ):
StoreDataVessel(da),
store_director(false)
{
  vecs=dynamic_cast<VectorMultiColvar*>( getAction() );
  plumed_assert( vecs );
  if( vecs->complexvec ) ncomponents=2*vecs->ncomponents;  
  else ncomponents = vecs->ncomponents;   

  completeSetup( 5, ncomponents ); myfvec.resize( ncomponents );
}

void StoreVectorsVessel::usedInFunction( const bool& store ){
  store_director=store; resize();
}

void StoreVectorsVessel::recompute( const unsigned& ivec, const unsigned& jstore ){
  plumed_dbg_assert( usingLowMem() ); 
  // Set the task we want to reperform
  setTaskToRecompute( ivec );
  // Reperform the task
  vecs->performTask();
  // Store the derivatives
  storeDerivativesLowMem( jstore );
  // Normalize the vector if it is required
  if( store_director ) normalizeVector( jstore ); 
  // Clear up afterwards
  vecs->clearAfterTask();

}

bool StoreVectorsVessel::calculate(){
  storeValues( vecs->getCurrentPositionInTaskList() );  // Store the values of the components of the vector

  if(!store_director) return true;
  if( !usingLowMem() ) normalizeVector( vecs->getCurrentPositionInTaskList() );
  else normalizeVector( -1 );  // Ensures vector components are normalized 
  return true;
}

void StoreVectorsVessel::normalizeVector( const int& jstore ){
  unsigned myelem = vecs->getCurrentPositionInTaskList();
  bool lowmemory = usingLowMem(); double norm2=0.0, norm;
  
  if( (lowmemory && jstore<0) || !lowmemory ){
     for(unsigned icomp=0;icomp<ncomponents;++icomp) norm2 += getComponent( myelem, icomp ) * getComponent( myelem, icomp );
     norm=sqrt( norm2 ); 
     for(unsigned icomp=0;icomp<ncomponents;++icomp){
        myfvec[icomp]=getComponent( myelem, icomp );      
        setComponent( myelem, icomp, getComponent( myelem, icomp ) / norm );
     }
  } else {
     for(unsigned icomp=0;icomp<ncomponents;++icomp){
        myfvec[icomp] = vecs->getElementValue( 5+icomp ); norm2 += myfvec[icomp] * myfvec[icomp];
     }
     norm=sqrt( norm2 );
  }
  double norm3=norm2*norm, weight = 1.0 / norm, wdf = -1.0 / norm3;

  if( !lowmemory ) {
      plumed_dbg_assert( jstore<getAction()->getFullNumberOfTasks() );
      for(unsigned ider=0;ider<getNumberOfDerivatives(myelem);++ider){
          double comp2=0.0; unsigned ibuf = myelem * ncomponents * getNumberOfDerivativeSpacesPerComponent() + 1 + ider;
          for(unsigned jcomp=0;jcomp<ncomponents;++jcomp){
              comp2  += myfvec[jcomp]*getBufferElement(ibuf);
              ibuf += getNumberOfDerivativeSpacesPerComponent();
          }
          ibuf = myelem * ncomponents * getNumberOfDerivativeSpacesPerComponent() + 1 + ider;
          for(unsigned jcomp=0;jcomp<ncomponents;++jcomp){
             setBufferElement( ibuf, weight*getBufferElement(ibuf) + wdf*comp2*myfvec[jcomp] );
             ibuf += getNumberOfDerivativeSpacesPerComponent();
          }
      }
  } else if( jstore>0 ) {
      unsigned maxder = vecs->getNumberOfDerivatives();
      for(unsigned ider=0;ider<getNumberOfDerivatives(jstore);++ider){
          double comp2=0.0; unsigned ibuf = jstore * ncomponents * maxder + ider;
          for(unsigned jcomp=0;jcomp<ncomponents;++jcomp){
              comp2 += myfvec[jcomp]*getLocalDerivative( ibuf );    
              ibuf += maxder;
          }
          ibuf = jstore * ncomponents * maxder + ider;
          for(unsigned jcomp=0;jcomp<ncomponents;++jcomp){
             setLocalDerivative( ibuf,  weight*getLocalDerivative( ibuf ) + wdf*comp2*myfvec[jcomp] ); 
             ibuf += maxder;
          }
      }
  }
}

void StoreVectorsVessel::chainRuleForComponent( const unsigned& icolv, const unsigned& jin, const unsigned& jout, const unsigned& base_cv_no, 
                                                const double& weight, multicolvar::MultiColvarFunction* funcout ){
  if( usingLowMem() ){
     unsigned ibuf = ( icolv*ncomponents + jin ) * getAction()->getNumberOfDerivatives();
     for(unsigned ider=0;ider<getNumberOfDerivatives(icolv);++ider){
         funcout->addStoredDerivative( jout, base_cv_no, getStoredIndex( icolv, ider ), weight*getLocalDerivative(ibuf+ider) );
     }                                         
  } else {
     unsigned ibuf = (icolv*ncomponents + jin ) * getNumberOfDerivativeSpacesPerComponent() + 1;
     for(unsigned ider=0;ider<getNumberOfDerivatives(icolv);++ider){
         funcout->addStoredDerivative( jout, base_cv_no, getStoredIndex( icolv, ider ), weight*getBufferElement(ibuf+ider) );
     }   
  }  
}

void StoreVectorsVessel::chainRuleForVector( const unsigned& icolv, const unsigned& jout, const unsigned& base_cv_no, 
                                             const std::vector<double>& df, multicolvar::MultiColvarFunction* funcout ){ 
   chainRule( icolv, df );
   for(unsigned ider=0;ider<getNumberOfDerivatives(icolv);++ider){
       funcout->addStoredDerivative( jout, base_cv_no, getStoredIndex( icolv, ider ), getFinalDerivative(ider) );
   }
}

}
}
