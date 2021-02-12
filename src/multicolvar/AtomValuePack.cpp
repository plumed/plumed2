/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2020 The plumed team
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
#include "AtomValuePack.h"
#include "CatomPack.h"
#include "tools/LinkCells.h"

namespace PLMD {
namespace multicolvar {

AtomValuePack::AtomValuePack( MultiValue& vals, MultiColvarBase const * mcolv ):
  myvals(vals),
  mycolv(mcolv),
  natoms(0),
  indices( vals.getIndices() ),
  sort_vector( vals.getSortIndices() ),
  myatoms( vals.getAtomVector() )
{
  if( indices.size()!=mcolv->getNumberOfAtoms() ) {
    indices.resize( mcolv->getNumberOfAtoms() );
    sort_vector.resize( mcolv->getNumberOfAtoms() );
    myatoms.resize( mcolv->getNumberOfAtoms() );
  }
}

unsigned AtomValuePack::setupAtomsFromLinkCells( const std::vector<unsigned>& cind, const Vector& cpos, const LinkCells& linkcells ) {
  if( cells_required.size()!=linkcells.getNumberOfCells() ) cells_required.resize( linkcells.getNumberOfCells() );
  // Build the list of cells that we need
  unsigned ncells_required=0; linkcells.addRequiredCells( linkcells.findMyCell( cpos ), ncells_required, cells_required );
  // Now build the list of atoms we need
  natoms=cind.size(); for(unsigned i=0; i<natoms; ++i) indices[i]=cind[i];
  linkcells.retrieveAtomsInCells( ncells_required, cells_required, natoms, indices );
//  linkcells.retrieveNeighboringAtoms( cpos, natoms, indices );
  for(unsigned i=0; i<natoms; ++i) myatoms[i]=mycolv->getPositionOfAtomForLinkCells( indices[i] ) - cpos;
  if( mycolv->usesPbc() ) mycolv->applyPbc( myatoms, natoms );
  return natoms;
}

void AtomValuePack::updateUsingIndices() {
  if( myvals.updateComplete() ) return;

  unsigned jactive=0;
  for(unsigned i=0; i<natoms; ++i) {
    unsigned base=3*indices[i];
    if( myvals.isActive( base ) ) { sort_vector[jactive]=indices[i]; jactive++; }
  }
  std::sort( sort_vector.begin(), sort_vector.begin()+jactive );

  myvals.emptyActiveMembers();
  for(unsigned i=0; i<jactive; ++i) {
    unsigned base=3*sort_vector[i]; // indices[i];
    myvals.putIndexInActiveArray( base );
    myvals.putIndexInActiveArray( base + 1 );
    myvals.putIndexInActiveArray( base + 2 );
  }
  unsigned nvir=3*mycolv->getNumberOfAtoms();
  if( myvals.isActive( nvir ) ) {
    for(unsigned i=0; i<9; ++i) myvals.putIndexInActiveArray( nvir + i );
  }
  myvals.completeUpdate();
}

void AtomValuePack::addComDerivatives( const int& ind, const Vector& der, const CatomPack& catom_der ) {
  if( ind<0 ) {
    for(unsigned ider=0; ider<catom_der.getNumberOfAtomsWithDerivatives(); ++ider) {
      unsigned jder=3*catom_der.getIndex(ider);
      myvals.addTemporyDerivative( jder+0, catom_der.getDerivative(ider,0,der) );
      myvals.addTemporyDerivative( jder+1, catom_der.getDerivative(ider,1,der) );
      myvals.addTemporyDerivative( jder+2, catom_der.getDerivative(ider,2,der) );
    }
  } else {
    for(unsigned ider=0; ider<catom_der.getNumberOfAtomsWithDerivatives(); ++ider) {
      unsigned jder=3*catom_der.getIndex(ider);
      myvals.addDerivative( ind, jder+0, catom_der.getDerivative(ider,0,der) );
      myvals.addDerivative( ind, jder+1, catom_der.getDerivative(ider,1,der) );
      myvals.addDerivative( ind, jder+2, catom_der.getDerivative(ider,2,der) );
    }
  }
}

}
}
