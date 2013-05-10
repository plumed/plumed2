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
#ifndef __PLUMED_multicolvar_MultiColvar_h
#define __PLUMED_multicolvar_MultiColvar_h

#include "MultiColvarBase.h"
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
/// Have atoms been read in
  bool readatoms;
  bool verbose_output;
/// Read in the various GROUP keywords
  void readGroupsKeyword( int& natoms );
/// Read in the various SPECIES keywords
  void readSpeciesKeyword( int& natoms );
/// Update the atoms request
  void requestAtoms();
protected:
/// Read in all the keywords that can be used to define atoms
  void readAtoms( int& natoms );
/// Read in ATOMS keyword
  void readAtomsLikeKeyword( const std::string key, int& natoms );
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
/// Resize all the dynamic arrays (used at neighbor list update time and during setup)
  void resizeDynamicArrays();
/// Get the position of atom iatom
  const Vector & getPosition(unsigned) const;
/// Calculate the multicolvar
  virtual void calculate();
/// Do the calculation
  double doCalculation( const unsigned& j );
/// Actually compute the colvar
  virtual double compute( const unsigned& j )=0;
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
};

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
  MultiColvarBase::addAtomsDerivatives( getAtomIndex(iatom), der );
}

inline
void MultiColvar::addAtomsDerivativeOfWeight( const unsigned& iatom, const Vector& wder ){
  MultiColvarBase::addAtomsDerivativeOfWeight( getAtomIndex(iatom), wder );
}

inline
void MultiColvar::addCentralAtomDerivatives( const unsigned& iatom, const Tensor& der ){
  MultiColvarBase::addCentralAtomDerivatives( getAtomIndex(iatom), der );
}

}
}

#endif
