/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2014 The plumed team
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
#include "AtomValuePack.h"
#include "CatomPack.h"
#include "tools/LinkCells.h"

namespace PLMD {
namespace multicolvar {

AtomValuePack::AtomValuePack( MultiValue& vals, MultiColvarBase const * mcolv ):
myvals(vals),
mycolv(mcolv),
indices( vals.getIndices() ),
sort_vector( vals.getSortIndices() ),
myatoms( vals.getAtomVector() )
{
  if( indices.size()!=mcolv->getNumberOfAtoms() ){ 
    indices.resize( mcolv->getNumberOfAtoms() ); 
    sort_vector.resize( mcolv->getNumberOfAtoms() ); 
    myatoms.resize( mcolv->getNumberOfAtoms() );
  }
}

unsigned AtomValuePack::setupAtomsFromLinkCells( const unsigned& cind, const Vector& cpos, const LinkCells& linkcells ){
  indices[0]=cind; natoms=1; linkcells.retrieveNeighboringAtoms( cpos, natoms, indices );
  myatoms[0]=Vector(0.0,0.0,0.0);
  for(unsigned i=1;i<natoms;++i) myatoms[i]=mycolv->getPositionOfAtomForLinkCells( indices[i] ) - cpos;
  if( mycolv->usesPbc() ) mycolv->applyPbc( myatoms, natoms );
  return natoms;
}

void AtomValuePack::updateUsingIndices(){
  if( myvals.updateComplete() ) return;

  unsigned jactive=0;
  for(unsigned i=0;i<natoms;++i){
      unsigned base=3*indices[i];
      if( myvals.isActive( base ) ){ sort_vector[jactive]=indices[i]; jactive++; } 
  }
  std::sort( sort_vector.begin(), sort_vector.begin()+jactive );

  myvals.emptyActiveMembers();
  for(unsigned i=0;i<jactive;++i){
     unsigned base=3*sort_vector[i]; // indices[i]; 
     myvals.putIndexInActiveArray( base );
     myvals.putIndexInActiveArray( base + 1 );
     myvals.putIndexInActiveArray( base + 2 ); 
  }
  unsigned nvir=3*mycolv->getNumberOfAtoms();
  if( myvals.isActive( nvir ) ){
      for(unsigned i=0;i<9;++i) myvals.putIndexInActiveArray( nvir + i );
  }
  myvals.completeUpdate();
}

void AtomValuePack::addComDerivatives( const unsigned& ind, const Vector& der, CatomPack& catom_der ){
  for(unsigned ider=0;ider<catom_der.getNumberOfAtomsWithDerivatives();++ider){
      unsigned jder=3*catom_der.getIndex(ider);
      myvals.addDerivative( ind, jder+0, catom_der.getDerivative(ider,0,der) );
      myvals.addDerivative( ind, jder+1, catom_der.getDerivative(ider,1,der) );
      myvals.addDerivative( ind, jder+2, catom_der.getDerivative(ider,2,der) );
  }
}

}
}
