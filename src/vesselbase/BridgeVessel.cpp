/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
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
#include "BridgeVessel.h"
#include "tools/Matrix.h"
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
      mynumerical_values.resize( getAction()->getNumberOfDerivatives()*myOutputValues->getNumberOfComponents() );
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
  myOutputAction->performTask();
  if( myOutputAction->thisval[1]<myOutputAction->getTolerance() ){
      myOutputAction->clearAfterTask();
      return ( !myOutputAction->contributorsAreUnlocked || myOutputAction->thisval[1]>=myOutputAction->getNLTolerance() );
  }
  bool keep=myOutputAction->calculateAllVessels();
  return ( !myOutputAction->contributorsAreUnlocked || keep );
}

void BridgeVessel::finish(){
  myOutputAction->finishComputations();
  if( myOutputAction->checkNumericalDerivatives() ){
     if ( inum<mynumerical_values.size() ){
         for(int i=0;i<myOutputValues->getNumberOfComponents();++i){
             mynumerical_values[inum]=myOutputValues->getOutputQuantity(i);
             inum++;
         }
         plumed_dbg_assert( inum<=mynumerical_values.size() );
     } else {
         plumed_assert( inum==mynumerical_values.size() );
     }
  } 
}

void BridgeVessel::completeNumericalDerivatives(){
  unsigned nextra = myOutputAction->getNumberOfDerivatives() - getAction()->getNumberOfDerivatives();
  Matrix<double> tmpder( myOutputValues->getNumberOfComponents(), nextra );
  ActionWithVessel* vval=dynamic_cast<ActionWithVessel*>( myOutputAction );
  for(unsigned i=0;i<nextra;++i){
      vval->bridgeVariable=i; getAction()->calculate();
      for(int j=0;j<myOutputValues->getNumberOfComponents();++j) tmpder(j,i) = myOutputValues->getOutputQuantity(j);
  }
  vval->bridgeVariable=nextra; getAction()->calculate(); 
  inum=0;  // Reset inum now that we have finished calling calculate
  std::vector<double> base( myOutputValues->getNumberOfComponents() );
  for(int j=0;j<myOutputValues->getNumberOfComponents();++j) base[j] = myOutputValues->getOutputQuantity(j);

  const double delta=sqrt(epsilon);
  ActionAtomistic* aa=dynamic_cast<ActionAtomistic*>( getAction() );
  unsigned nvals=myOutputValues->getNumberOfComponents();
  for(unsigned j=0;j<nvals;++j) ( myOutputValues->copyOutput(j) )->clearDerivatives();   

  if( aa ){
      ActionWithArguments* aarg=dynamic_cast<ActionWithArguments*>( getAction() );
      plumed_assert( !aarg ); Tensor box=aa->getBox(); 
      unsigned natoms=aa->getNumberOfAtoms();
      for(unsigned j=0;j<nvals;++j){
          double ref=( myOutputValues->copyOutput(j) )->get();
          if( ( myOutputValues->copyOutput(j) )->hasDerivatives() ){
              for(unsigned i=0;i<3*natoms;++i){
                  double d=( mynumerical_values[i*nvals+j] - ref)/delta;
                  ( myOutputValues->copyOutput(j) )->addDerivative(i,d);
              }
              Tensor virial;
              for(int i=0;i<3;i++) for(int k=0;k<3;k++){
                 virial(i,k)=( mynumerical_values[ nvals*(3*natoms + 3*i + k) + j ]-ref)/delta;
              }
              virial=-matmul(box.transpose(),virial);
              for(int i=0;i<3;i++) for(int k=0;k<3;k++) ( myOutputValues->copyOutput(j) )->addDerivative(3*natoms+3*k+i,virial(k,i));
          }
      }
  } else {
      plumed_merror("not implemented or tested yet");
//      unsigned nder=myOutputAction->getNumberOfDerivatives();
//      for(unsigned j=0;j<nvals;++j){
//          double ref=( myOutputValues->copyOutput(j) )->get();
//          if( ( myOutputValues->copyOutput(j) )->hasDerivatives() ){
//              for(unsigned i=0;i<nder;++i){
//                  double d=( mynumerical_values[i*nvals+j] - ref)/delta;
//                  ( myOutputValues->copyOutput(j) )->addDerivative(i,d);
//              }
//          }
//      }
  }
  // Add the derivatives wrt to the local quantities we are working with
  for(unsigned j=0;j<nvals;++j){
     unsigned k=0;
     for(unsigned i=getAction()->getNumberOfDerivatives();i<myOutputAction->getNumberOfDerivatives();++i){
        ( myOutputValues->copyOutput(j) )->addDerivative( i, (tmpder(j,k)-base[j])/sqrt(epsilon) ); k++;
     }
  }
}

bool BridgeVessel::applyForce( std::vector<double>& outforces ){
  bool hasforce=false; outforces.assign(outforces.size(),0.0);
  unsigned nextra = myOutputAction->getNumberOfDerivatives() - getAction()->getNumberOfDerivatives();
  std::vector<double> eforces( nextra, 0.0 );
  for(unsigned i=0;i<myOutputAction->getNumberOfVessels();++i){
     if( ( myOutputAction->getPntrToVessel(i) )->applyForce( forces ) ){
         hasforce=true;
         for(unsigned j=0;j<outforces.size();++j) outforces[j]+=forces[j];
         for(unsigned j=0;j<nextra;++j) eforces[j]+=forces[ outforces.size()+j ];
     }
  }
  if(hasforce) myOutputAction->applyBridgeForces( eforces );
  return hasforce;
}

}
}

