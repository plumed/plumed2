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
#include "MultiColvar.h"
#include "MultiColvarFunction.h"
#include "StoreCentralAtomsVessel.h"

#define CATOMS_MAXATOMS 2

namespace PLMD {
namespace multicolvar {

StoreCentralAtomsVessel::StoreCentralAtomsVessel( const vesselbase::VesselOptions& da ):
Vessel(da),
nspace(1 + 3*CATOMS_MAXATOMS)
{
  mycolv=dynamic_cast<MultiColvar*>( getAction() );
  plumed_assert( mycolv );
  // Resize the active derivative lists
  active_atoms.resize( mycolv->colvar_atoms.size() );
}

void StoreCentralAtomsVessel::resize(){
  unsigned nfunc=mycolv->colvar_atoms.size(); 
  unsigned bsize=0; start.resize( nfunc +1 );
  for(unsigned i=0;i<nfunc;++i){
      start[i] = bsize;
      bsize += 3*( 1 + 3*CATOMS_MAXATOMS ); 
  }
  start[nfunc]=bsize;
  resizeBuffer( bsize ); 
  // Update the active_derivative lists
  for(unsigned i=0;i<active_atoms.size();++i){
     active_atoms[i].clear();
     for(unsigned j=0;j<mycolv->getNumberOfAtoms();++j){
         active_atoms[i].addIndexToList(j);
     }
  }
}

bool StoreCentralAtomsVessel::calculate(){
  Vector catom_pos=mycolv->retrieveCentralAtomPos( false );

  // Store the value
  unsigned ibuf=start[mycolv->current]; 
  for(unsigned i=0;i<3;++i){ addToBufferElement( ibuf, catom_pos[i] ); ibuf+=nspace; }
  plumed_dbg_assert( ibuf==start[mycolv->current+1] );

  // Chect there is sufficent space for the derivatives
  plumed_massert( 3*mycolv->atomsWithCatomDer.getNumberActive()<=3*CATOMS_MAXATOMS, 
    "increase CATOMS_MAXATOMS in StoreCentralAtomsVessel");

  // Store the derivatives
  active_atoms[mycolv->current].deactivateAll();
  for(unsigned j=0;j<mycolv->atomsWithCatomDer.getNumberActive();++j){
      unsigned n=mycolv->atomsWithCatomDer[j]; 
      active_atoms[mycolv->current].activate(n); 
      ibuf = start[mycolv->current] + 1 + 3*j;       
      for(unsigned i=0;i<3;++i){
          for(unsigned k=0;k<3;++k) addToBufferElement( ibuf+k, mycolv->central_derivs[n](i,k) );
          ibuf += nspace;
      }
  }

  return true;
}

void StoreCentralAtomsVessel::finish(){
  mpi_gatherActiveMembers( comm, active_atoms );
}

Vector StoreCentralAtomsVessel::getPosition( const unsigned& ivec ) const {
  plumed_dbg_assert( ivec<mycolv->colvar_atoms.size() );
  unsigned pos=start[ivec]; Vector mypos;
  for(unsigned i=0;i<3;++i){
      mypos[i] = getBufferElement( pos ); 
      pos+=nspace;      
  }
  plumed_dbg_assert( pos==start[ivec+1] );
  return mypos;
}

void StoreCentralAtomsVessel::addAtomsDerivatives( const unsigned& iatom, const Vector& df, MultiColvarFunction* funcout ) const {
  plumed_dbg_assert( iatom<mycolv->colvar_atoms.size() );

  Vector thisder;
  for(unsigned ider=0;ider<active_atoms[iatom].getNumberActive();++ider){
      unsigned ibuf=start[iatom] + 1 + 3*ider; thisder.zero();
      for(unsigned jcomp=0;jcomp<3;++jcomp){
          thisder[0]+=df[jcomp]*getBufferElement(ibuf);
          thisder[1]+=df[jcomp]*getBufferElement(ibuf+1);
          thisder[2]+=df[jcomp]*getBufferElement(ibuf+2);
          ibuf+=nspace;
      }
      unsigned jatom = active_atoms[iatom][ider];
      funcout->addAtomsDerivatives( jatom, thisder );
  }
}

void StoreCentralAtomsVessel::addAtomsDerivativeOfWeight( const unsigned& iatom, const Vector& df, MultiColvarFunction* funcout ) const {
  plumed_dbg_assert( iatom<mycolv->colvar_atoms.size() );
 
  Vector thisder; 
  for(unsigned ider=0;ider<active_atoms[iatom].getNumberActive();++ider){
      unsigned ibuf=start[iatom] + 1 + 3*ider; thisder.zero();
      for(unsigned jcomp=0;jcomp<3;++jcomp){
          thisder[0]+=df[jcomp]*getBufferElement(ibuf);
          thisder[1]+=df[jcomp]*getBufferElement(ibuf+1);
          thisder[2]+=df[jcomp]*getBufferElement(ibuf+2);
          ibuf+=nspace;
      }   
      unsigned jatom = active_atoms[iatom][ider];
      funcout->addAtomsDerivativeOfWeight( jatom, thisder );
  }   
}

void StoreCentralAtomsVessel::addDerivativeOfCentralAtomPos( const unsigned& iatom, const Tensor& df, MultiColvarFunction* funcout ) const {
  plumed_dbg_assert( iatom<mycolv->colvar_atoms.size() );

  Tensor thisder;
  for(unsigned ider=0;ider<active_atoms[iatom].getNumberActive();++ider){
     unsigned ibuf=start[iatom] + 1 + 3*ider; thisder.zero();
     for(unsigned jcomp=0;jcomp<3;++jcomp){
         for(unsigned kcomp=0;kcomp<3;++kcomp){
             for(unsigned k=0;k<3;++k) thisder(jcomp,kcomp)+=df(jcomp,k)*getBufferElement(ibuf+k*nspace+kcomp);
         }
     }
     unsigned jatom = active_atoms[iatom][ider];
     funcout->addCentralAtomDerivatives( jatom, thisder );
  }
}     


}
}
