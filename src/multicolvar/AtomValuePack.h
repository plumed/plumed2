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
#ifndef __PLUMED_multicolvar_AtomValuePack_h
#define __PLUMED_multicolvar_AtomValuePack_h

#include "tools/MultiValue.h"
#include "MultiColvarBase.h"

namespace PLMD {
namespace multicolvar {

class AtomValuePack {
  friend class MultiColvarBase;
private:
/// Copy of the values that we are adding to
  MultiValue& myvals;
/// Copy of the underlying multicolvar
  MultiColvarBase const * mycolv;
/// Number of atoms at the moment
  unsigned natoms;
/// Atom indices
  std::vector<unsigned>& indices;
/// This holds atom positions
  std::vector<Vector>& myatoms;
public:
  AtomValuePack( MultiValue& vals, MultiColvarBase const * mcolv );
/// Set the number of atoms
  void setNumberOfAtoms( const unsigned& );
///
  void addAtomsDerivatives( const unsigned&, const unsigned&, const Vector& );
///
  void setAtomIndex( const unsigned& j, const unsigned& ind );
///
  void setAtom( const unsigned& j, const unsigned& ind );
///
  void makeWhole();
///
  unsigned getNumberOfAtoms() const ;
///
  unsigned getNumberOfDerivatives() const ;
/// Get the position of the ith atom
  Vector& getPosition( const unsigned& );
///
  void setValue( const unsigned&, const double& );
///
  void addValue( const unsigned& ival, const double& vv );
///
  double getValue( const unsigned& ) const ;
///
  void setBoxDerivativesNoPbc(const unsigned& ival);
///
  void addBoxDerivatives( const unsigned&, const Tensor& );
///
  void updateUsingIndices();
};

inline
void AtomValuePack::setNumberOfAtoms( const unsigned& nat ) {
  natoms=nat;
}

inline
unsigned AtomValuePack::getNumberOfAtoms() const {
  return natoms;
}

inline
unsigned AtomValuePack::getNumberOfDerivatives() const {
  return myvals.getNumberOfDerivatives();
}

inline
void AtomValuePack::setAtomIndex( const unsigned& j, const unsigned& ind ) {
  plumed_dbg_assert( j<natoms ); indices[j]=ind;
}

inline
void AtomValuePack::setAtom( const unsigned& j, const unsigned& ind ) {
  setAtomIndex( j, ind ); myatoms[j]=mycolv->getPosition( ind );
}

inline
Vector& AtomValuePack::getPosition( const unsigned& iatom ) {
  plumed_dbg_assert( iatom<natoms );
  return myatoms[iatom];
}

inline
void AtomValuePack::setValue( const unsigned& ival, const double& vv ) {
  myvals.setValue( (mycolv->getPntrToOutput(ival))->getPositionInStream(), vv );
}

inline
void AtomValuePack::addValue( const unsigned& ival, const double& vv ) {
  myvals.addValue( (mycolv->getPntrToOutput(ival))->getPositionInStream(), vv );
}

inline
double AtomValuePack::getValue( const unsigned& ival ) const {
  return myvals.get( (mycolv->getPntrToOutput(ival))->getPositionInStream() );
}

inline
void AtomValuePack::addAtomsDerivatives( const unsigned& ival, const unsigned& jder, const Vector& der ) {
  if( mycolv->doNotCalculateDerivatives() ) return;
  plumed_dbg_assert( jder<natoms );
  unsigned jval=(mycolv->getPntrToOutput(ival))->getPositionInStream();
  myvals.addDerivative( jval, 3*indices[jder] + 0, der[0] );
  myvals.addDerivative( jval, 3*indices[jder] + 1, der[1] );
  myvals.addDerivative( jval, 3*indices[jder] + 2, der[2] );
}

inline
void AtomValuePack::addBoxDerivatives( const unsigned& ival, const Tensor& vir ) {
  if( mycolv->doNotCalculateDerivatives() ) return;
  unsigned nvir=3*mycolv->getNumberOfAtoms(), jval=(mycolv->getPntrToOutput(ival))->getPositionInStream();
  for(unsigned i=0; i<3; ++i) for(unsigned j=0; j<3; ++j) myvals.addDerivative( jval, nvir + 3*i+j, vir(i,j) );
}

}
}
#endif

