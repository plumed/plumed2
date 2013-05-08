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
#include "BridgeVessel.h"
#include "core/ActionWithArguments.h"

namespace PLMD {
namespace vesselbase {

BridgeVessel::BridgeVessel( const VesselOptions& da ):
Vessel(da),
inum(0)
{
  resizeBuffer(0);
}

void BridgeVessel::resize(){
  forces.resize( myOutputAction->getNumberOfDerivatives() );
  myOutputAction->resizeFunctions();
  if( myOutputAction->checkNumericalDerivatives() ){
      mynumerical_values.resize( myOutputAction->getNumberOfDerivatives()*myOutputValues->getNumberOfComponents() );
      inum=0;
  }
}

void BridgeVessel::setOutputAction( ActionWithVessel* myact ){
  myOutputAction=myact;
  myOutputValues=dynamic_cast<ActionWithValue*>( myact );
  plumed_assert( myOutputValues );
}

std::string BridgeVessel::description(){
  plumed_merror("I shouldn't end up here");
}

void BridgeVessel::prepare(){
  myOutputAction->doJobsRequiredBeforeTaskList();
}

bool BridgeVessel::calculate(){
  unsigned j;
  myOutputAction->performTask(j);
  if( myOutputAction->thisval[1]<myOutputAction->getTolerance() ){
      myOutputAction->clearAfterTask();
      return ( !myOutputAction->areContributorsUnlocked() || myOutputAction->thisval[1]>=myOutputAction->getNLTolerance() );
  }
  bool keep=myOutputAction->calculateAllVessels();
  return ( !myOutputAction->areContributorsUnlocked() || keep );
}

void BridgeVessel::finish(){
  myOutputAction->finishComputations();
  if( myOutputAction->checkNumericalDerivatives() ){
     if ( inum<mynumerical_values.size() ){
         for(unsigned i=0;i<myOutputValues->getNumberOfComponents();++i){
             mynumerical_values[inum]=myOutputValues->getOutputQuantity(i);
             inum++;
         }
     }
  } 
}

void BridgeVessel::completeNumericalDerivatives(){
  inum=0; const double delta=sqrt(epsilon);
  ActionAtomistic* aa=dynamic_cast<ActionAtomistic*>( getAction() );
  unsigned nder=myOutputAction->getNumberOfDerivatives();
  for(unsigned j=0;j<myOutputValues->getNumberOfComponents();++j) ( myOutputValues->copyOutput(j) )->clearDerivatives();   

  if( aa ){
      ActionWithArguments* aarg=dynamic_cast<ActionWithArguments*>( getAction() );
      plumed_assert( !aarg ); Tensor box=aa->getBox();
      for(unsigned j=0;j<myOutputValues->getNumberOfComponents();++j){
          double ref=( myOutputValues->copyOutput(j) )->get();
          if( ( myOutputValues->copyOutput(j) )->hasDerivatives() ){
              for(unsigned i=0;i<nder-9;++i){
                  double d=( mynumerical_values[j*nder+i] - ref)/delta;
                  ( myOutputValues->copyOutput(j) )->addDerivative(i,d);
              }
              Tensor virial;
              for(int i=0;i<3;i++) for(int k=0;k<3;k++) virial(i,k)=( mynumerical_values[nder-9+3*k+i]-ref)/delta;
              virial=-1.0*matmul(box.transpose(),virial);
              for(int i=0;i<3;i++) for(int k=0;k<3;k++) ( myOutputValues->copyOutput(j) )->addDerivative(nder-9+3*k+i,virial(i,k));
          }
      }
  } else {
      for(unsigned j=0;j<myOutputValues->getNumberOfComponents();++j){
          double ref=( myOutputValues->copyOutput(j) )->get();
          if( ( myOutputValues->copyOutput(j) )->hasDerivatives() ){
              for(unsigned i=0;i<nder;++i){
                  double d=( mynumerical_values[j*nder+i] - ref)/delta;
                  ( myOutputValues->copyOutput(j) )->addDerivative(i,d);
              }
          }
      }
  }
}

bool BridgeVessel::applyForce( std::vector<double>& outforces ){
  bool hasforce=false; outforces.assign(outforces.size(),0.0);
  for(unsigned i=0;i<myOutputAction->getNumberOfVessels();++i){
     if( ( myOutputAction->getPntrToVessel(i) )->applyForce( forces ) ){
         hasforce=true;
         for(unsigned j=0;j<outforces.size();++j) outforces[j]+=forces[j];
     }
  }
  return hasforce;
}

}
}

