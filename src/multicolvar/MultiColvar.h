/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2014 The plumed team
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
  void readGroupsKeyword( int& natoms, std::vector<AtomNumber>& all_atoms );
/// Read in the various SPECIES keywords
  void readSpeciesKeyword( int& natoms, std::vector<AtomNumber>& all_atoms );
protected:
/// Read in all the keywords that can be used to define atoms
  void readAtoms( int& natoms );
/// Read in ATOMS keyword
  void readAtomsLikeKeyword( const std::string & key, int& natoms, std::vector<AtomNumber>& all_atoms );
/// Read two group keywords
  void readTwoGroups( const std::string& key1, const std::string& key2, std::vector<AtomNumber>& all_atoms );
/// Read three groups
  void readThreeGroups( const std::string& key1, const std::string& key2, const std::string& key3, const bool& allow2, std::vector<AtomNumber>& all_atoms );
/// Add a collective variable
  void addColvar( const std::vector<unsigned>& newatoms );
/// Add some derivatives for an atom 
  void addAtomsDerivatives(const int&,const Vector&);
/// Set the derivative of the weight (used in MEAN and HISTOGRAM)
  void addAtomsDerivativeOfWeight( const unsigned& i, const Vector& wder );
/// Add derivatives to the central atom position
  void addCentralAtomDerivatives( const unsigned& iatom, const Tensor& der );
public:
  MultiColvar(const ActionOptions&);
  ~MultiColvar(){}
  static void registerKeywords( Keywords& keys );
/// Get the position of atom iatom
  const Vector & getPosition(unsigned) const;
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
/// This is used in MultiColvarBase only - it is used to setup the link cells
  Vector getPositionOfAtomForLinkCells( const unsigned& iatom );
/// Atoms are always active
  bool isCurrentlyActive( const unsigned& code ){ return true; }
};

inline
Vector MultiColvar::getPositionOfAtomForLinkCells( const unsigned& iatom ){
  return ActionAtomistic::getPosition( iatom );
}

inline
const Vector & MultiColvar::getPosition( unsigned iatom ) const {
  return ActionAtomistic::getPosition( current_atoms[iatom] );
}

inline
double MultiColvar::getMass(unsigned iatom ) const {
  return ActionAtomistic::getMass( current_atoms[iatom] );
}

inline
double MultiColvar::getCharge(unsigned iatom ) const {
  return ActionAtomistic::getCharge( current_atoms[iatom] );
}

inline
AtomNumber MultiColvar::getAbsoluteIndex(unsigned iatom) const {
  return ActionAtomistic::getAbsoluteIndex( current_atoms[iatom] );
}

inline
void MultiColvar::addAtomsDerivatives(const int& iatom, const Vector& der){
  MultiColvarBase::addAtomsDerivatives( 0, current_atoms[iatom], der );
}

inline
void MultiColvar::addAtomsDerivativeOfWeight( const unsigned& iatom, const Vector& wder ){
  MultiColvarBase::addAtomsDerivatives( 1, current_atoms[iatom], wder );
}

inline
void MultiColvar::addCentralAtomDerivatives( const unsigned& iatom, const Tensor& der ){
  MultiColvarBase::addCentralAtomDerivatives( current_atoms[iatom], der );
}

}
}

#endif
