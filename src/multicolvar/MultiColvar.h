/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
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
#ifndef __PLUMED_multicolvar_MultiColvar_h
#define __PLUMED_multicolvar_MultiColvar_h

#include "MultiColvarBase.h"
#include "tools/SwitchingFunction.h"
#include <vector>

#define PLUMED_MULTICOLVAR_INIT(ao) Action(ao),MultiColvar(ao)

namespace PLMD {
namespace multicolvar {

/**
\ingroup INHERIT
This is the abstract base class to use for creating distributions of colvars and functions
thereof, whtin it there is \ref AddingAMultiColvar "information" as to how to go implementing these types of actions.
*/

class MultiColvar : public MultiColvarBase {
private:
/// Do we want lots of details in the output
  bool verbose_output;
/// Read in the various GROUP keywords
  void readGroupsKeyword( int& natoms );
/// Read in the various SPECIES keywords
  void readSpeciesKeyword( int& natoms );
protected:
/// Read in all the keywords that can be used to define atoms
  void readAtoms( int& natoms );
/// Read in ATOMS keyword
  void readAtomsLikeKeyword( const std::string & key, int& natoms );
/// Read two group keywords
  void readTwoGroups( const std::string& key1, const std::string& key2 );
/// Read three groups
  void readThreeGroups( const std::string& key1, const std::string& key2, const std::string& key3, const bool& allow2 );
/// This is used to make neighbor list update fast when three atoms are involved in the colvar (e.g. ANGLES, WATERBRIDGE)
  void threeBodyNeighborList( const SwitchingFunction& sf );
/// Add a collective variable
  void addColvar( const std::vector<unsigned>& newatoms );
/// Add some derivatives for an atom 
  void addAtomsDerivatives(const int&,const Vector&);
/// Set the derivative of the weight (used in MEAN and HISTOGRAM)
  void addAtomsDerivativeOfWeight( const unsigned& i, const Vector& wder );
/// Add derivatives to the central atom position
  void addCentralAtomDerivatives( const unsigned& iatom, const Tensor& der );
/// Get the index of an atom
   unsigned getAtomIndex( const unsigned& iatom ) const ;
public:
  MultiColvar(const ActionOptions&);
  ~MultiColvar(){}
  static void registerKeywords( Keywords& keys );
/// Resize all the dynamic arrays (used at neighbor list update time and during setup)
//  virtual void resizeDynamicArrays();
/// Get the position of atom iatom
  const Vector & getPosition(unsigned) const;
/// This is used in VectorMultiColvar.  In there we use a MultiColvar as if it is 
/// a MultiColvarFunction so we can calculate functions of vectors that output functions
  virtual void calculationsRequiredBeforeNumericalDerivatives(){}
/// Finish the update of the task list
  void finishTaskListUpdate();
/// Calculate the multicolvar
  virtual void calculate();
/// Update the atoms that have derivatives
  void updateActiveAtoms();
/// Calculate the position of the central atom
  Vector calculateCentralAtomPosition();
/// Get the position of the central atom
  virtual Vector getCentralAtom()=0;
/// Get the mass of atom iatom
  double getMass(unsigned) const ;
/// Get the charge of atom iatom
  double getCharge(unsigned) const ;
/// Get the absolute index of atom iatom
  AtomNumber getAbsoluteIndex(unsigned) const ;
/// Get base quantity index
  unsigned getBaseQuantityIndex( const unsigned& code );
};

inline
unsigned MultiColvar::getBaseQuantityIndex( const unsigned& code ){
  return all_atoms.linkIndex( code );
}

inline
unsigned MultiColvar::getAtomIndex( const unsigned& iatom ) const {
  plumed_dbg_assert( iatom<natomsper );
  return current_atoms[iatom];
//  return all_atoms.linkIndex( colvar_atoms[current][iatom] );
}

inline
const Vector & MultiColvar::getPosition( unsigned iatom ) const {
  return ActionAtomistic::getPosition( getAtomIndex(iatom) );
}

inline
double MultiColvar::getMass(unsigned iatom ) const {
  return ActionAtomistic::getMass( getAtomIndex(iatom) );
}

inline
double MultiColvar::getCharge(unsigned iatom ) const {
  return ActionAtomistic::getCharge( getAtomIndex(iatom) );
}

inline
AtomNumber MultiColvar::getAbsoluteIndex(unsigned iatom) const {
  return ActionAtomistic::getAbsoluteIndex( getAtomIndex(iatom) );
}

inline
void MultiColvar::addAtomsDerivatives(const int& iatom, const Vector& der){
  MultiColvarBase::addAtomsDerivatives( 0, getAtomIndex(iatom), der );
}

inline
void MultiColvar::addAtomsDerivativeOfWeight( const unsigned& iatom, const Vector& wder ){
  MultiColvarBase::addAtomsDerivatives( 1, getAtomIndex(iatom), wder );
}

inline
void MultiColvar::addCentralAtomDerivatives( const unsigned& iatom, const Tensor& der ){
  MultiColvarBase::addCentralAtomDerivatives( getAtomIndex(iatom), der );
}

}
}

#endif
