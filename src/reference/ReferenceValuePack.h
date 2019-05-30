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
#ifndef __PLUMED_reference_ReferenceValuePack_h
#define __PLUMED_reference_ReferenceValuePack_h

#include "tools/MultiValue.h"

namespace PLMD {

class ReferenceValuePack {
  friend class MultiDomainRMSD;
  friend class OptimalRMSD;
private:
/// Was the virial set
  bool boxWasSet;
///
  unsigned numberOfArgs;
///
  bool oind_set;
  unsigned oind;
/// Copy of the values that we are adding to
  MultiValue& myvals;
/// Ths list of atom indices
  std::vector<unsigned>& atom_indices;
/// Are we using pca
  bool pca;
/// A vector of vectors to save us some overhead for vector resizes
  std::vector<Vector> centeredpos;
///
  std::vector<Vector> displacement;
///
  std::vector<Tensor> rot;
///
  Matrix< std::vector<Vector> >  DRotDPos;
public:
  ReferenceValuePack( const unsigned& nargs, const unsigned& natoms, MultiValue& vals );
///
  void resize( const unsigned& nargs, const unsigned& natoms );
///
  void clear();
///
  unsigned getNumberOfDerivatives() const ;
///
  unsigned getNumberOfArguments() const ;
///
  unsigned getNumberOfAtoms() const ;
///
  void setAtomIndex( const unsigned& iatom, const unsigned& jindex );
///
  unsigned getAtomIndex( const unsigned& iatom ) const ;
///
  void copyScaledDerivatives( const unsigned& from, const double& scalef, const MultiValue& tvals );
///
  void addArgumentDerivatives( const unsigned& iarg, const double& der );
///
  void setArgumentDerivatives( const unsigned& iarg, const double& der );
///
  void setAtomDerivatives( const unsigned& jder, const Vector& der );
///
  void addAtomDerivatives( const unsigned& iatom, const Vector& der );
///
  void addBoxDerivatives( const Tensor& vir );
///
  bool updateComplete() const ;
///
  void updateDynamicLists();
///
  void scaleAllDerivatives( const double& scalef );
///
  void setValIndex( const unsigned& ind );
///
  void moveDerivatives( const unsigned& from, const unsigned& to );
///
  bool virialWasSet() const ;
///
  Vector getAtomDerivative( const unsigned& iatom ) const ;
///
  double getArgumentDerivative( const unsigned& ival ) const ;
///
  Tensor getBoxDerivatives() const ;
///
  bool calcUsingPCAOption() const ;
///
  void switchOnPCAOption();
///
  std::vector<Vector>& getAtomVector();
///
  std::vector<Vector>& getAtomsDisplacementVector();
};

inline
unsigned ReferenceValuePack::getNumberOfDerivatives() const {
  return myvals.getNumberOfDerivatives();
}

inline
unsigned ReferenceValuePack::getNumberOfArguments() const {
  return numberOfArgs;
}

inline
unsigned ReferenceValuePack::getNumberOfAtoms() const {
  return atom_indices.size();
}

inline
void ReferenceValuePack::setAtomIndex( const unsigned& iatom, const unsigned& jindex ) {
  plumed_dbg_assert( iatom<atom_indices.size() ); atom_indices[iatom]=jindex;
}

inline
unsigned ReferenceValuePack::getAtomIndex( const unsigned& iatom ) const {
  plumed_dbg_assert( iatom<atom_indices.size() );
  return atom_indices[iatom];
}

inline
void ReferenceValuePack::addArgumentDerivatives( const unsigned& iarg, const double& der ) {
  plumed_dbg_assert( iarg<numberOfArgs && oind_set ); myvals.addDerivative( oind, iarg, der );
}

inline
void ReferenceValuePack::setArgumentDerivatives( const unsigned& iarg, const double& der ) {
  plumed_dbg_assert( iarg<numberOfArgs && oind_set ); myvals.setDerivative( oind, iarg, der );
}

inline
bool ReferenceValuePack::updateComplete() const {
  return myvals.updateComplete();
}

inline
void ReferenceValuePack::setAtomDerivatives( const unsigned& jder, const Vector& der ) {
  plumed_dbg_assert( oind_set && jder<atom_indices.size() );
  myvals.setDerivative( oind, numberOfArgs + 3*atom_indices[jder] + 0, der[0] );
  myvals.setDerivative( oind, numberOfArgs + 3*atom_indices[jder] + 1, der[1] );
  myvals.setDerivative( oind, numberOfArgs + 3*atom_indices[jder] + 2, der[2] );
}

inline
void ReferenceValuePack::addAtomDerivatives( const unsigned& jder, const Vector& der ) {
  plumed_dbg_assert( oind_set && jder<atom_indices.size() );
  myvals.addDerivative( oind, numberOfArgs + 3*atom_indices[jder] + 0, der[0] );
  myvals.addDerivative( oind, numberOfArgs + 3*atom_indices[jder] + 1, der[1] );
  myvals.addDerivative( oind, numberOfArgs + 3*atom_indices[jder] + 2, der[2] );
}

inline
void ReferenceValuePack::addBoxDerivatives( const Tensor& vir ) {
  plumed_dbg_assert( oind_set && atom_indices.size()>0 );
  boxWasSet=true; unsigned nbase = myvals.getNumberOfDerivatives() - 9;
  for(unsigned i=0; i<3; ++i) for(unsigned j=0; j<3; ++j) myvals.addDerivative( oind, nbase + 3*i + j, vir(i,j) );
}

inline
void ReferenceValuePack::setValIndex( const unsigned& ind ) {
  oind=ind; oind_set=true;
}

inline
bool ReferenceValuePack::virialWasSet() const {
  return boxWasSet;
}

inline
Vector ReferenceValuePack::getAtomDerivative( const unsigned& iatom ) const {
  Vector tmp; plumed_dbg_assert( oind_set && iatom<atom_indices.size() );
  tmp[0]=myvals.getDerivative( oind, numberOfArgs + 3*atom_indices[iatom] + 0 );
  tmp[1]=myvals.getDerivative( oind, numberOfArgs + 3*atom_indices[iatom] + 1 );
  tmp[2]=myvals.getDerivative( oind, numberOfArgs + 3*atom_indices[iatom] + 2 );
  return tmp;
}

inline
double ReferenceValuePack::getArgumentDerivative( const unsigned& ival ) const {
  plumed_dbg_assert( oind_set && ival<numberOfArgs );
  return myvals.getDerivative( oind, ival );
}

inline
Tensor ReferenceValuePack::getBoxDerivatives() const {
  plumed_dbg_assert( oind_set && boxWasSet ); Tensor tvir; unsigned nbase = myvals.getNumberOfDerivatives() - 9;
  for(unsigned i=0; i<3; ++i) for(unsigned j=0; j<3; ++j) tvir(i,j)=myvals.getDerivative( oind, nbase + 3*i + j );
  return tvir;
}

inline
bool ReferenceValuePack::calcUsingPCAOption() const {
  return pca;
}

inline
void ReferenceValuePack::switchOnPCAOption() {
  pca=true;
}

inline
std::vector<Vector>& ReferenceValuePack::getAtomVector() {
  return myvals.getAtomVector();
}

inline
std::vector<Vector>& ReferenceValuePack::getAtomsDisplacementVector() {
  return displacement;
}

}

#endif
