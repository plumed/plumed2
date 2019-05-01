/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2019 The plumed team
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

class LinkCells;

namespace multicolvar {

class CatomPack;

class AtomValuePack {
  friend class MultiColvarBase;
  friend class LocalAverage;
private:
/// Copy of the values that we are adding to
  MultiValue& myvals;
/// Copy of the underlying multicolvar
  MultiColvarBase const * mycolv;
/// Number of atoms at the moment
  unsigned natoms;
/// Atom indices
  std::vector<unsigned>& indices;
/// This is used to sort the atom indices
  std::vector<unsigned>& sort_vector;
/// This holds atom positions
  std::vector<Vector>& myatoms;
/// This is stuff for link cells
  std::vector<unsigned> cells_required;
///
  void addAtomsDerivatives( const unsigned&, const unsigned&, const Vector& );
///
  void addTemporyAtomsDerivatives( const unsigned& jder, const Vector& der );
public:
  AtomValuePack( MultiValue& vals, MultiColvarBase const * mcolv );
/// Set the number of atoms
  void setNumberOfAtoms( const unsigned& );
/// Set the index for one of the atoms
  void setIndex( const unsigned&, const unsigned& );
///
  void setAtomIndex( const unsigned& j, const unsigned& ind );
///
  void setAtom( const unsigned& j, const unsigned& ind );
///
  unsigned setupAtomsFromLinkCells( const std::vector<unsigned>& cind, const Vector& cpos, const LinkCells& linkcells );
///
  unsigned getIndex( const unsigned& j ) const ;
///
  unsigned getNumberOfAtoms() const ;
///
  unsigned getNumberOfDerivatives() const ;
/// Get the position of the ith atom
  Vector& getPosition( const unsigned& );
/// Get the absolute index of the ith atom in the list
  AtomNumber getAbsoluteIndex( const unsigned& j ) const ;
///
  void setValue( const unsigned&, const double& );
///
  void addValue( const unsigned& ival, const double& vv );
///
  double getValue( const unsigned& ) const ;
///
  void addBoxDerivatives( const unsigned&, const Tensor& );
///
  void addTemporyBoxDerivatives( const Tensor& vir );
///
  void updateUsingIndices();
///
  void updateDynamicList();
///
  void addComDerivatives( const int&, const Vector&, CatomPack& );
///
  MultiValue& getUnderlyingMultiValue();
///
  void addDerivative( const unsigned&, const unsigned&, const double& );
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
void AtomValuePack::setIndex( const unsigned& j, const unsigned& ind ) {
  plumed_dbg_assert( j<natoms ); indices[j]=ind;
}

inline
void AtomValuePack::setAtomIndex( const unsigned& j, const unsigned& ind ) {
  plumed_dbg_assert( j<natoms ); indices[j]=ind;
}

inline
void AtomValuePack::setAtom( const unsigned& j, const unsigned& ind ) {
  setAtomIndex( j, ind ); myatoms[j]=mycolv->getPositionOfAtomForLinkCells( ind );
}

inline
unsigned AtomValuePack::getIndex( const unsigned& j ) const {
  plumed_dbg_assert( j<natoms ); return indices[j];
}

inline
AtomNumber AtomValuePack::getAbsoluteIndex( const unsigned& j ) const {
  plumed_dbg_assert( j<natoms ); unsigned jatom=indices[j];
  if( mycolv->atom_lab[jatom].first>0 ) {
    unsigned mmc=mycolv->atom_lab[jatom].first - 1;
    return (mycolv->mybasemulticolvars[mmc])->getAbsoluteIndexOfCentralAtom( mycolv->atom_lab[jatom].second );
  }
  return mycolv->getAbsoluteIndex( mycolv->atom_lab[jatom].second );
}

inline
Vector& AtomValuePack::getPosition( const unsigned& iatom ) {
  plumed_dbg_assert( iatom<natoms );
  return myatoms[iatom];
}

inline
void AtomValuePack::setValue( const unsigned& ival, const double& vv ) {
  myvals.setValue( ival, vv );
}

inline
void AtomValuePack::addValue( const unsigned& ival, const double& vv ) {
  myvals.addValue( ival, vv );
}

inline
double AtomValuePack::getValue( const unsigned& ival ) const {
  return myvals.get( ival );
}

inline
void AtomValuePack::addDerivative( const unsigned& ival, const unsigned& jder, const double& der ) {
  myvals.addDerivative( ival, jder, der );
}

inline
void AtomValuePack::addAtomsDerivatives( const unsigned& ival, const unsigned& jder, const Vector& der ) {
  plumed_dbg_assert( jder<natoms );
  myvals.addDerivative( ival, 3*indices[jder] + 0, der[0] );
  myvals.addDerivative( ival, 3*indices[jder] + 1, der[1] );
  myvals.addDerivative( ival, 3*indices[jder] + 2, der[2] );
}

inline
void AtomValuePack::addTemporyAtomsDerivatives( const unsigned& jder, const Vector& der ) {
  plumed_dbg_assert( jder<natoms );
  myvals.addTemporyDerivative( 3*indices[jder] + 0, der[0] );
  myvals.addTemporyDerivative( 3*indices[jder] + 1, der[1] );
  myvals.addTemporyDerivative( 3*indices[jder] + 2, der[2] );
}

inline
void AtomValuePack::addTemporyBoxDerivatives( const Tensor& vir ) {
  unsigned nvir=3*mycolv->getNumberOfAtoms();
  for(unsigned i=0; i<3; ++i) for(unsigned j=0; j<3; ++j) myvals.addTemporyDerivative( nvir + 3*i+j, vir(i,j) );
}

inline
void AtomValuePack::addBoxDerivatives( const unsigned& ival, const Tensor& vir ) {
  unsigned nvir=3*mycolv->getNumberOfAtoms();
  for(unsigned i=0; i<3; ++i) for(unsigned j=0; j<3; ++j) myvals.addDerivative( ival, nvir + 3*i+j, vir(i,j) );
}

inline
void AtomValuePack::updateDynamicList() {
  if( myvals.updateComplete() ) return;
  myvals.updateDynamicList();
}

inline
MultiValue& AtomValuePack::getUnderlyingMultiValue() {
  return myvals;
}

}
}
#endif

