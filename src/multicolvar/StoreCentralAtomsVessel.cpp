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
#include "MultiColvar.h"
#include "StoreCentralAtomsVessel.h"

#define CATOMS_MAXATOMS 2

namespace PLMD {
namespace multicolvar {

StoreCentralAtomsVessel::StoreCentralAtomsVessel( const vesselbase::VesselOptions& da ):
Vessel(da),
wasforced(false),
nspace(1 + 3*CATOMS_MAXATOMS)
{
  mycolv=dynamic_cast<MultiColvar*>( getAction() );
  plumed_assert( mycolv );
}

void StoreCentralAtomsVessel::resize(){
  unsigned nfunc=mycolv->getNumberOfFunctionsInAction();
  unsigned bsize=0; start.resize( nfunc +1 );
  for(unsigned i=0;i<nfunc;++i){
      start[i] = bsize;
      bsize += 3*( 1 + 3*CATOMS_MAXATOMS ); 
  }
  start[nfunc]=bsize;
  resizeBuffer( bsize ); 
  forces.resize( mycolv->getNumberOfDerivatives() );
  // Update the active_derivative lists
  mycolv->buildDerivativeIndexArrays( active_der );
}

void StoreCentralAtomsVessel::prepare(){
  wasforced=false;
  forces.assign( forces.size(), 0.0 );
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
  active_der[mycolv->current].deactivateAll();
  for(unsigned j=0;j<mycolv->atomsWithCatomDer.getNumberActive();++j){
      unsigned n=mycolv->atomsWithCatomDer[j]; 
      ibuf = start[mycolv->current] + 1 + 3*j;       
      for(unsigned i=0;i<3;++i){
          active_der[mycolv->current].activate(3*n+i); 
          for(unsigned k=0;k<3;++k) addToBufferElement( ibuf+k, mycolv->central_derivs[n](i,k) );
          ibuf += nspace;
      }
  }

  return true;
}

void StoreCentralAtomsVessel::finish(){
  mpi_gatherActiveMembers( comm, active_der );
}

Vector StoreCentralAtomsVessel::getPosition( const unsigned& ivec ) const {
  plumed_dbg_assert( ivec<mycolv->getNumberOfFunctionsInAction() );
  unsigned pos=start[ivec]; Vector mypos;
  for(unsigned i=0;i<3;++i){
      mypos[i] = getBufferElement( pos ); 
      pos+=nspace;      
  }
  plumed_dbg_assert( pos==start[ivec+1] );
  return mypos;
}

void StoreCentralAtomsVessel::addForces( const std::vector<double>& ff){
  plumed_dbg_assert( ff.size()==forces.size() );
  wasforced=true;
  for(unsigned i=0;i<forces.size();++i) forces[i]+=ff[i];
}

bool StoreCentralAtomsVessel::applyForce(std::vector<double>& ff){
  plumed_dbg_assert( ff.size()==forces.size() );
  if(wasforced){
    for(unsigned i=0;i<forces.size();++i) ff[i]=forces[i];
  }
  return wasforced;
}

void StoreCentralAtomsVessel::chainRuleForCentralAtom( const unsigned& iatom, const unsigned& iderno, const Vector& df, vesselbase::ActionWithVessel* act) const {
  plumed_dbg_assert( iatom<mycolv->getNumberOfFunctionsInAction() );

  unsigned nder=mycolv->getNumberOfDerivatives();
  for(unsigned ider=0;ider<active_der[iatom].getNumberActive();++ider){
      unsigned ibuf=start[iatom] + 1 + ider;   
      double thisder=0.0;
      for(unsigned jcomp=0;jcomp<3;++jcomp){
          thisder+=df[jcomp]*getBufferElement(ibuf);
          ibuf+=nspace;
      }
      unsigned jder = iderno*nder + active_der[iatom][ider];
      act->addElementDerivative( jder, thisder );
  }
}

}
}
