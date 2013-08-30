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
#include "vesselbase/ActionWithVessel.h"
#include "StoreColvarVessel.h"
#include "MultiColvarBase.h"
#include "tools/DynamicList.h"

#define MAXATOMS 300

namespace PLMD {
namespace multicolvar{

void StoreColvarVessel::registerKeywords( Keywords& keys ){
  Vessel::registerKeywords( keys );
  plumed_assert( keys.size()==0 );
}

StoreColvarVessel::StoreColvarVessel( const vesselbase::VesselOptions& da ):
Vessel(da)
{
  mycolv=dynamic_cast<MultiColvarBase*>( getAction() );
  plumed_assert( mycolv );
  diffweight=mycolv->weightHasDerivatives;
  // Resize the active derivative lists
  active_atoms.resize( mycolv->colvar_atoms.size() + mycolv->colvar_atoms.size()*MAXATOMS );
}

void StoreColvarVessel::resize(){
  unsigned nfunc=mycolv->colvar_atoms.size();
  bufsize=0; start.resize( nfunc +1 );
  for(unsigned i=0;i<nfunc;++i){
      start[i] = bufsize;
      bufsize += 1 + 3*MAXATOMS + 9; 
  }
  start[nfunc]=bufsize;
  resizeBuffer( 2*bufsize ); local_resizing();
}

void StoreColvarVessel::prepare(){
  active_atoms.assign( active_atoms.size(),0 );
}

bool StoreColvarVessel::calculate(){
  plumed_massert( mycolv->atoms_with_derivatives.getNumberActive()<=MAXATOMS,
        "Error increase MAXATOMS in StoreColvarVessel");

  unsigned ibuf=start[mycolv->current]; 
  setBufferElement( ibuf, mycolv->getElementValue(0) ); ibuf++;
  unsigned atom_der_index=0; unsigned atom_der_index_start = mycolv->colvar_atoms.size() + mycolv->current*MAXATOMS;
  for(unsigned j=0;j<mycolv->atoms_with_derivatives.getNumberActive();++j){
     unsigned iatom=mycolv->atoms_with_derivatives[j]; 
     active_atoms[ atom_der_index_start + atom_der_index ] = iatom; atom_der_index++;
     unsigned ider=3*iatom;
     setBufferElement( ibuf, mycolv->getElementDerivative(ider) ); ider++; ibuf++; 
     setBufferElement( ibuf, mycolv->getElementDerivative(ider) ); ider++; ibuf++;
     setBufferElement( ibuf, mycolv->getElementDerivative(ider) ); ibuf++;
  } 
  active_atoms[mycolv->current] = atom_der_index;
  unsigned ivir=3*mycolv->getNumberOfAtoms();
  for(unsigned j=0;j<9;++j){
     setBufferElement( ibuf, mycolv->getElementDerivative(ivir) ); ivir++; ibuf++;
  }
  ibuf=bufsize+start[mycolv->current]; 
  setBufferElement( ibuf, mycolv->getElementValue(1) ); ibuf++;
  if(diffweight){
    unsigned nder=mycolv->getNumberOfDerivatives();
    for(unsigned j=0;j<mycolv->atoms_with_derivatives.getNumberActive();++j){ 
       unsigned ider=nder+3*mycolv->atoms_with_derivatives[j];
       setBufferElement( ibuf, mycolv->getElementDerivative(ider) ); ider++; ibuf++; 
       setBufferElement( ibuf, mycolv->getElementDerivative(ider) ); ider++; ibuf++;
       setBufferElement( ibuf, mycolv->getElementDerivative(ider) ); ibuf++;
    }
    unsigned ivir=nder + 3*mycolv->getNumberOfAtoms();
    for(unsigned j=0;j<9;++j){
       setBufferElement( ibuf, mycolv->getElementDerivative(ivir) ); ivir++; ibuf++;
    }
  }
  return true;
}

void StoreColvarVessel::finish(){
  comm.Sum(active_atoms);
  performCalculationUsingAllValues();
}

void StoreColvarVessel::addDerivatives( const unsigned& ival, double& pref, Value* value_out ){
  unsigned atom_der_index_start = mycolv->colvar_atoms.size() + ival*MAXATOMS;
  for(unsigned i=0;i<active_atoms[ival];++i){
      unsigned jbuf=start[ival] + 1 + 3*i; unsigned aind=3*active_atoms[ atom_der_index_start + i ];
      value_out->addDerivative( aind, pref*getBufferElement(jbuf) ); jbuf++; aind++;
      value_out->addDerivative( aind, pref*getBufferElement(jbuf) ); jbuf++; aind++;
      value_out->addDerivative( aind, pref*getBufferElement(jbuf) ); 
  }
  unsigned jbuf=start[ival] + 1 + 3*active_atoms[ival];
  unsigned nder=3*mycolv->getNumberOfAtoms();
  for(unsigned i=0;i<9;++i){
     value_out->addDerivative( nder, pref*getBufferElement(jbuf) ); nder++; jbuf++;
  }
}

void StoreColvarVessel::addWeightDerivatives( const unsigned& ival, double& pref, Value* value_out ){
  if(!diffweight) return;

  unsigned atom_der_index_start = mycolv->colvar_atoms.size() + ival*MAXATOMS;
  for(unsigned i=0;i<active_atoms[ival];++i){
      unsigned jbuf=bufsize+start[ival] + 1 +3*i; unsigned aind=3*active_atoms[ atom_der_index_start + i ];
      value_out->addDerivative( aind, pref*getBufferElement(jbuf) ); jbuf++; aind++;
      value_out->addDerivative( aind, pref*getBufferElement(jbuf) ); jbuf++; aind++;
      value_out->addDerivative( aind, pref*getBufferElement(jbuf) ); 
  }
  unsigned jbuf=bufsize+start[ival] + 1 + 3*active_atoms[ival];
  unsigned nder=3*mycolv->getNumberOfAtoms();
  for(unsigned i=0;i<9;++i){ 
     value_out->addDerivative( nder, pref*getBufferElement(jbuf) ); nder++; jbuf++;
  }
}


}
}

