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
#include "vesselbase/VesselRegister.h"
#include "vesselbase/ActionWithVessel.h"
#include "tools/DynamicList.h"
#include "MultiColvar.h"
#include "MultiColvarFunction.h"
#include "StoreCentralAtomsVessel.h"

namespace PLMD {
namespace multicolvar {

StoreCentralAtomsVessel::StoreCentralAtomsVessel( const vesselbase::VesselOptions& da ):
Vessel(da)
{
  mycolv=dynamic_cast<MultiColvar*>( getAction() );
  plumed_assert( mycolv );
  // Get the number of atoms with derivatives for each central atom position
  nder_atoms = mycolv->getNumberOfAtomsInCentralAtomDerivatives();
  // Resize the active derivative lists
  active_atoms.resize( mycolv->colvar_atoms.size() + mycolv->colvar_atoms.size()*nder_atoms );
  // The number of scalars stored for each central atom (1 value + natoms*3 derivatives )
  nspace = 1 + 3*nder_atoms;
}

void StoreCentralAtomsVessel::resize(){
  unsigned nfunc=mycolv->colvar_atoms.size(); 
  resizeBuffer( nfunc*3*nspace ); 
}

void StoreCentralAtomsVessel::prepare(){
  active_atoms.assign( active_atoms.size(),0 );
}

bool StoreCentralAtomsVessel::calculate(){
  Vector catom_pos=mycolv->retrieveCentralAtomPos();

  // Store the value
  unsigned ibuf=3*nspace*mycolv->current;     // start[mycolv->current]; 
  for(unsigned i=0;i<3;++i){ addToBufferElement( ibuf, catom_pos[i] ); ibuf+=nspace; }
  plumed_dbg_assert( ibuf==3*(mycolv->current+1)*nspace );

  // Chect there is sufficent space for the derivatives
  plumed_dbg_assert( 3*mycolv->atomsWithCatomDer.getNumberActive()==3*nder_atoms ); 

  // Store the derivatives
  unsigned atom_der_index=0; unsigned atom_der_index_start = mycolv->colvar_atoms.size() + mycolv->current*nder_atoms;
  for(unsigned j=0;j<mycolv->atomsWithCatomDer.getNumberActive();++j){
      unsigned n=mycolv->atomsWithCatomDer[j]; 
      active_atoms[ atom_der_index_start + atom_der_index ] = n; atom_der_index++;
      ibuf = 3*nspace*mycolv->current + 1 + 3*j;       
      for(unsigned i=0;i<3;++i){
          for(unsigned k=0;k<3;++k) addToBufferElement( ibuf+k, mycolv->central_derivs[n](i,k) );
          ibuf += nspace;
      }
  }
  active_atoms[mycolv->current] = atom_der_index;

  return true;
}

void StoreCentralAtomsVessel::finish(){
  comm.Sum( &active_atoms[0], active_atoms.size() );
}

Vector StoreCentralAtomsVessel::getPosition( const unsigned& ivec ) const {
  plumed_dbg_assert( ivec<mycolv->colvar_atoms.size() );
  unsigned pos=3*nspace*ivec; Vector mypos;
  for(unsigned i=0;i<3;++i){
      mypos[i] = getBufferElement( pos ); 
      pos+=nspace;      
  }
  plumed_dbg_assert( pos==3*(ivec+1)*nspace );
  return mypos;
}

void StoreCentralAtomsVessel::addAtomsDerivatives( const unsigned& iatom, const Vector& df, MultiColvarFunction* funcout ) const {
  plumed_dbg_assert( iatom<mycolv->colvar_atoms.size() );

  Vector thisder; unsigned atom_der_index_start = mycolv->colvar_atoms.size() + iatom*nder_atoms; 
  for(unsigned ider=0;ider<active_atoms[iatom];++ider){
      unsigned ibuf=3*iatom*nspace + 1 + 3*ider; thisder.zero();
      for(unsigned jcomp=0;jcomp<3;++jcomp){
          thisder[0]+=df[jcomp]*getBufferElement(ibuf);
          thisder[1]+=df[jcomp]*getBufferElement(ibuf+1);
          thisder[2]+=df[jcomp]*getBufferElement(ibuf+2);
          ibuf+=nspace;
      }
      unsigned jatom = active_atoms[ atom_der_index_start + ider ];
      funcout->addAtomsDerivatives( jatom, thisder );
  }
}

void StoreCentralAtomsVessel::addAtomsDerivativeOfWeight( const unsigned& iatom, const Vector& df, MultiColvarFunction* funcout ) const {
  plumed_dbg_assert( iatom<mycolv->colvar_atoms.size() );
 
  Vector thisder; unsigned atom_der_index_start = mycolv->colvar_atoms.size() + iatom*nder_atoms; 
  for(unsigned ider=0;ider<active_atoms[iatom];++ider){
      unsigned ibuf=3*iatom*nspace + 1 + 3*ider; thisder.zero();
      for(unsigned jcomp=0;jcomp<3;++jcomp){
          thisder[0]+=df[jcomp]*getBufferElement(ibuf);
          thisder[1]+=df[jcomp]*getBufferElement(ibuf+1);
          thisder[2]+=df[jcomp]*getBufferElement(ibuf+2);
          ibuf+=nspace;
      }   
      unsigned jatom = active_atoms[ atom_der_index_start + ider];
      funcout->addAtomsDerivativeOfWeight( jatom, thisder );
  }   
}

void StoreCentralAtomsVessel::addDerivativeOfCentralAtomPos( const unsigned& iatom, const Tensor& df, MultiColvarFunction* funcout ) const {
  plumed_dbg_assert( iatom<mycolv->colvar_atoms.size() );

  Tensor thisder; unsigned atom_der_index_start = mycolv->colvar_atoms.size() + iatom*nder_atoms;
  for(unsigned ider=0;ider<active_atoms[iatom];++ider){
     unsigned ibuf=3*iatom*nspace + 1 + 3*ider; thisder.zero();
     for(unsigned jcomp=0;jcomp<3;++jcomp){
         for(unsigned kcomp=0;kcomp<3;++kcomp){
             for(unsigned k=0;k<3;++k) thisder(jcomp,kcomp)+=df(jcomp,k)*getBufferElement(ibuf+k*nspace+kcomp);
         }
     }
     unsigned jatom = active_atoms[ atom_der_index_start + ider];
     funcout->addCentralAtomDerivatives( jatom, thisder );
  }
}

unsigned StoreCentralAtomsVessel::getDerivativeIndex( const unsigned& iatom, const unsigned& jder ) const {
  plumed_dbg_assert( iatom<mycolv->colvar_atoms.size() && jder<nder_atoms );
  unsigned atom_der_index_start = mycolv->colvar_atoms.size() + iatom*nder_atoms;
  return  active_atoms[ atom_der_index_start + jder ];
}

Tensor StoreCentralAtomsVessel::getDerivatives( const unsigned& iatom, const unsigned& jder ) const {
  plumed_dbg_assert( iatom<mycolv->colvar_atoms.size() && jder<nder_atoms );
  Tensor der; unsigned ibuf = 3*iatom*nspace + 1 + 3*jder;
  for(unsigned i=0;i<3;++i){
      for(unsigned k=0;k<3;++k) der(i,k) = getBufferElement( ibuf + k );
      ibuf += nspace;
  }
  return der;
}     

}
}
