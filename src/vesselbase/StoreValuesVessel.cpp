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
#include "ActionWithVessel.h"
#include "StoreValuesVessel.h"

#define MAXDERIVATIVES 300

namespace PLMD {
namespace vesselbase{

void StoreValuesVessel::registerKeywords( Keywords& keys ){
  Vessel::registerKeywords( keys );
  plumed_assert( keys.size()==0 );
}

StoreValuesVessel::StoreValuesVessel( const VesselOptions& da ):
Vessel(da)
{
  diffweight=getAction()->weightHasDerivatives;
}

void StoreValuesVessel::resize(){
  ActionWithVessel* aa=getAction();
  unsigned nfunc=aa->getNumberOfFunctionsInAction();
  bufsize=0; start.resize( nfunc +1 );
  for(unsigned i=0;i<nfunc;++i){
      start[i] = bufsize;
      bufsize += 1 + MAXDERIVATIVES; 
  }
  start[nfunc]=bufsize;
  resizeBuffer( 2*bufsize ); local_resizing();
  // Build the arrays of indexes
  getAction()->buildDerivativeIndexArrays( active_der );
}

bool StoreValuesVessel::calculate(){
  ActionWithVessel* aa=getAction();
  unsigned ibuf=start[aa->current]; setBufferElement( ibuf, aa->getElementValue(0) ); ibuf++;
  for(unsigned ider=getAction()->getFirstDerivativeToMerge();ider<getAction()->getNumberOfDerivatives();ider=getAction()->getNextDerivativeToMerge(ider)){ 
     active_der[aa->current].activate(ider);
     setBufferElement( ibuf, aa->getElementDerivative(ider) ); ibuf++; 
  } 
  ibuf=bufsize+start[aa->current]; setBufferElement( ibuf, aa->getElementValue(1) ); ibuf++;
  if(diffweight){
    unsigned nder=aa->getNumberOfDerivatives();
    for(unsigned ider=getAction()->getFirstDerivativeToMerge();ider<getAction()->getNumberOfDerivatives();ider=getAction()->getNextDerivativeToMerge(ider)){ 
       active_der[aa->current].activate(ider);
       setBufferElement( ibuf, aa->getElementDerivative(nder+ider) ); ibuf++; 
    }
  }
  return true;
}

void StoreValuesVessel::finish(){
  mpi_gatherActiveMembers( comm, active_der );
  performCalculationUsingAllValues();
}

void StoreValuesVessel::addDerivatives( const unsigned& ival, double& pref, Value* value_out ){
  unsigned jbuf=start[ival]+1;
  for(unsigned i=0;i<active_der[ival].getNumberActive();++i){
      value_out->addDerivative( active_der[ival][i], pref*getBufferElement(jbuf) ); jbuf++;
  }
}

void StoreValuesVessel::addWeightDerivatives( const unsigned& ival, double& pref, Value* value_out ){
  if(!diffweight) return;

  unsigned jbuf=bufsize+start[ival]+1;
  for(unsigned i=0;i<active_der[ival].getNumberActive();++i){
      value_out->addDerivative( active_der[ival][i], pref*getBufferElement(jbuf) ); jbuf++;
  }
}


}
}

