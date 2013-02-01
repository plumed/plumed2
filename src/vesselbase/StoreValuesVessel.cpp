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
      bufsize += 1 + aa->getNumberOfDerivatives(i); 
  }
  start[nfunc]=bufsize;
  resizeBuffer( 2*bufsize ); local_resizing();
}

bool StoreValuesVessel::calculate(){
  ActionWithVessel* aa=getAction();
  unsigned ibuf=start[aa->current]; setBufferElement( ibuf, aa->getElementValue(0) ); ibuf++;
  for(unsigned ider=0;ider<getAction()->nderivatives;++ider){ setBufferElement( ibuf, aa->getElementDerivative(ider) ); ibuf++; } 
  plumed_dbg_assert( ibuf==start[aa->current+1] );
  ibuf=bufsize+start[aa->current]; setBufferElement( ibuf, aa->getElementValue(1) ); ibuf++;
  if(diffweight){
    unsigned nder=aa->getNumberOfDerivatives();
    for(unsigned ider=0;ider<getAction()->nderivatives;++ider){ setBufferElement( ibuf, aa->getElementDerivative(nder+ider) ); ibuf++; }
    plumed_dbg_assert( ibuf==(bufsize+start[aa->current+1]) );
  }
  return true;
}

void StoreValuesVessel::addDerivatives( const unsigned& ival, double& pref, Value* value_out ){
  unsigned j=0; ActionWithVessel* aa=getAction();
  for(unsigned i=start[ival]+1;i<start[ival+1];++i){
      value_out->addDerivative( aa->getOutputDerivativeIndex( ival, j ), pref*getBufferElement(i) ); j++;
  }
}

void StoreValuesVessel::addWeightDerivatives( const unsigned& ival, double& pref, Value* value_out ){
  if(!diffweight) return;

  unsigned j=0; ActionWithVessel* aa=getAction();
  for(unsigned i=bufsize+start[ival]+1;i<bufsize+start[ival+1];++i){
      value_out->addDerivative( aa->getOutputDerivativeIndex( ival, j ), pref*getBufferElement(i) ); j++;
  }
}


}
}

