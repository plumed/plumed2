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
#include "BridgeVessel.h"
#include "tools/Matrix.h"
#include "core/ActionWithArguments.h"
#include "StoreDataVessel.h"

namespace PLMD {
namespace vesselbase {

BridgeVessel::BridgeVessel( const VesselOptions& da ):
  Vessel(da),
  inum(0),
// in_normal_calculate(false)
  myOutputAction(NULL),
  myOutputValues(NULL),
  my_tmp_val(0,0)
{
}

void BridgeVessel::resize() {
  if( myOutputAction->checkNumericalDerivatives() ) {
    mynumerical_values.resize( getAction()->getNumberOfDerivatives()*myOutputValues->getNumberOfComponents() );
    inum=0;
  }
  // This bit ensures that we can store data in a bridge function if needs be
  // Notice we need to resize der_list in the underlying action as this is called
  // from a bridge
  if( myOutputAction->mydata ) {
    unsigned dsize=(myOutputAction->mydata)->getSizeOfDerivativeList();
    if( getAction()->der_list.size()!=dsize ) getAction()->der_list.resize( dsize );
  }
  unsigned tmp=0; resizeBuffer( myOutputAction->getSizeOfBuffer( tmp ) );
}

void BridgeVessel::setOutputAction( ActionWithVessel* myact ) {
  ActionWithValue* checkme=dynamic_cast<ActionWithValue*>( getAction() );
  plumed_massert( checkme, "vessel in bridge must inherit from ActionWithValue");

  myOutputAction=myact;
  myOutputValues=dynamic_cast<ActionWithValue*>( myact );
  plumed_massert( myOutputValues, "bridging vessel must inherit from ActionWithValue");
}

std::string BridgeVessel::description() {
  plumed_merror("I shouldn't end up here");
}

void BridgeVessel::prepare() {
  myOutputAction->doJobsRequiredBeforeTaskList();
}

void BridgeVessel::setBufferStart( unsigned& start ) {
  myOutputAction->getSizeOfBuffer( start );
}

MultiValue& BridgeVessel::transformDerivatives( const unsigned& current, MultiValue& invals, MultiValue& outvals ) {
  if( outvals.getNumberOfValues()!=myOutputAction->getNumberOfQuantities() ||
      outvals.getNumberOfDerivatives()!=myOutputAction->getNumberOfDerivatives() ) {
    outvals.resize( myOutputAction->getNumberOfQuantities(), myOutputAction->getNumberOfDerivatives() );
  }
  myOutputAction->transformBridgedDerivatives( current, invals, outvals );
  return outvals;
}

void BridgeVessel::calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const {
  // in_normal_calculate=true;
  if( myvals.get(0)<myOutputAction->getTolerance() ) return;
  myOutputAction->calculateAllVessels( current, myvals, myvals, buffer, der_list );
  return;
}

void BridgeVessel::finish( const std::vector<double>& buffer ) {
  myOutputAction->finishComputations( buffer );
  if( myOutputAction->checkNumericalDerivatives() ) {
    if ( inum<mynumerical_values.size() ) {
      for(int i=0; i<myOutputValues->getNumberOfComponents(); ++i) {
        mynumerical_values[inum]=myOutputValues->getOutputQuantity(i);
        inum++;
      }
      plumed_dbg_assert( inum<=mynumerical_values.size() );
    } else {
      plumed_assert( inum==mynumerical_values.size() );
    }
  }
}

void BridgeVessel::completeNumericalDerivatives() {
  unsigned nextra = myOutputAction->getNumberOfDerivatives() - getAction()->getNumberOfDerivatives();
  Matrix<double> tmpder( myOutputValues->getNumberOfComponents(), nextra );
  ActionWithVessel* vval=dynamic_cast<ActionWithVessel*>( myOutputAction );
  for(unsigned i=0; i<nextra; ++i) {
    vval->bridgeVariable=i; getAction()->calculate();
    for(int j=0; j<myOutputValues->getNumberOfComponents(); ++j) tmpder(j,i) = myOutputValues->getOutputQuantity(j);
  }
  vval->bridgeVariable=nextra; getAction()->calculate();
  plumed_assert( inum==mynumerical_values.size() ); inum=0;  // Reset inum now that we have finished calling calculate
  std::vector<double> base( myOutputValues->getNumberOfComponents() );
  for(int j=0; j<myOutputValues->getNumberOfComponents(); ++j) base[j] = myOutputValues->getOutputQuantity(j);

  const double delta=sqrt(epsilon);
  ActionAtomistic* aa=dynamic_cast<ActionAtomistic*>( getAction() );
  unsigned nvals=myOutputValues->getNumberOfComponents();
  for(unsigned j=0; j<nvals; ++j) ( myOutputValues->copyOutput(j) )->clearDerivatives();

  if( aa ) {
    ActionWithArguments* aarg=dynamic_cast<ActionWithArguments*>( getAction() );
    plumed_assert( !aarg ); Tensor box=aa->getBox();
    unsigned natoms=aa->getNumberOfAtoms();
    for(unsigned j=0; j<nvals; ++j) {
      double ref=( myOutputValues->copyOutput(j) )->get();
      if( ( myOutputValues->copyOutput(j) )->getNumberOfDerivatives()>0 ) {
        for(unsigned i=0; i<3*natoms; ++i) {
          double d=( mynumerical_values[i*nvals+j] - ref)/delta;
          ( myOutputValues->copyOutput(j) )->addDerivative(i,d);
        }
        Tensor virial;
        for(int i=0; i<3; i++) for(int k=0; k<3; k++) {
            virial(i,k)=( mynumerical_values[ nvals*(3*natoms + 3*i + k) + j ]-ref)/delta;
          }
        virial=-matmul(box.transpose(),virial);
        for(int i=0; i<3; i++) for(int k=0; k<3; k++) ( myOutputValues->copyOutput(j) )->addDerivative(3*natoms+3*k+i,virial(k,i));
      }
    }
  } else {
    plumed_merror("not implemented or tested yet");
//      unsigned nder=myOutputAction->getNumberOfDerivatives();
//      for(unsigned j=0;j<nvals;++j){
//          double ref=( myOutputValues->copyOutput(j) )->get();
//              for(unsigned i=0;i<nder;++i){
//                  double d=( mynumerical_values[i*nvals+j] - ref)/delta;
//                  ( myOutputValues->copyOutput(j) )->addDerivative(i,d);
//              }
//          }
//      }
  }
  // Add the derivatives wrt to the local quantities we are working with
  for(unsigned j=0; j<nvals; ++j) {
    unsigned k=0;
    for(unsigned i=getAction()->getNumberOfDerivatives(); i<myOutputAction->getNumberOfDerivatives(); ++i) {
      ( myOutputValues->copyOutput(j) )->addDerivative( i, (tmpder(j,k)-base[j])/sqrt(epsilon) ); k++;
    }
  }
}

bool BridgeVessel::applyForce( std::vector<double>& outforces ) {
  bool hasforce=false; outforces.assign(outforces.size(),0.0);
  unsigned ndertot = myOutputAction->getNumberOfDerivatives();
  unsigned nextra = ndertot - getAction()->getNumberOfDerivatives();
  std::vector<double> forces( ndertot ), eforces( nextra, 0.0 );
  for(unsigned i=0; i<myOutputAction->getNumberOfVessels(); ++i) {
    if( ( myOutputAction->getPntrToVessel(i) )->applyForce( forces ) ) {
      hasforce=true;
      for(unsigned j=0; j<outforces.size(); ++j) outforces[j]+=forces[j];
      for(unsigned j=0; j<nextra; ++j) eforces[j]+=forces[ outforces.size()+j ];
    }
  }
  if(hasforce) myOutputAction->applyBridgeForces( eforces );
  return hasforce;
}

void BridgeVessel::copyTaskFlags() {
  myOutputAction->deactivateAllTasks();
  for(unsigned i=0; i<getAction()->nactive_tasks; ++i) myOutputAction->taskFlags[ getAction()->indexOfTaskInFullList[i] ] = 1;
  myOutputAction->lockContributors();
}

MultiValue& BridgeVessel::getTemporyMultiValue() {
  return my_tmp_val;
}

}
}

