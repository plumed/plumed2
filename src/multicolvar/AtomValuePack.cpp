/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2017 The plumed team
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

namespace PLMD {
namespace multicolvar {

AtomValuePack::AtomValuePack( MultiValue& vals, MultiColvarBase const * mcolv ):
  myvals(vals),
  mycolv(mcolv),
  natoms(0),
  indices( vals.getIndices() ),
  sort_vector( vals.getSortIndices() ),
  myatoms( vals.getFirstAtomVector() )
{
  if( indices.size()!=mcolv->getNumberOfAtoms() ) {
    indices.resize( mcolv->getNumberOfAtoms() );
    sort_vector.resize( mcolv->getNumberOfAtoms() );
    myatoms.resize( mcolv->getNumberOfAtoms() );
  }
}

void AtomValuePack::makeWhole(){
  for(unsigned j=0; j<myatoms.size()-1; ++j) {
    const Vector & first (myatoms[j]);
    Vector & second (myatoms[j+1]);
    second=first+mycolv->pbcDistance(first,second);
  }
}

void AtomValuePack::updateUsingIndices() {
#ifdef DNDEBUG
  for(unsigned j=0;j<mycolv->getNumberOfComponents();++j) plumed_dbg_assert( myvals.getNumberActive((mycolv->getPntrToOutput(j))->getPositionInStream())==0 );
#endif  

  for(unsigned i=0; i<natoms; ++i) {
    unsigned base=3*indices[i]; // indices[i];
    for(unsigned j=0;j<mycolv->getNumberOfComponents();++j){
         myvals.updateIndex( (mycolv->getPntrToOutput(j))->getPositionInStream(), base );
         myvals.updateIndex( (mycolv->getPntrToOutput(j))->getPositionInStream(), base + 1 );
         myvals.updateIndex( (mycolv->getPntrToOutput(j))->getPositionInStream(), base + 2 );
    }
  }
  unsigned nvir=3*mycolv->getNumberOfAtoms();
  for(unsigned j=0;j<mycolv->getNumberOfComponents();++j){ 
      for(unsigned i=0; i<9; ++i) myvals.updateIndex( (mycolv->getPntrToOutput(j))->getPositionInStream(), nvir + i );
  }
}

void AtomValuePack::setBoxDerivativesNoPbc(const unsigned& ival) {
  if( mycolv->doNotCalculateDerivatives() ) return;
  Tensor virial; 
  for(unsigned i=0; i<natoms; i++) virial-=Tensor(getPosition(i), Vector(myvals.getDerivative(ival,3*indices[i]+0),
                                                                         myvals.getDerivative(ival,3*indices[i]+1), 
                                                                         myvals.getDerivative(ival,3*indices[i]+2)));
  addBoxDerivatives(ival,virial);
}

}
}
