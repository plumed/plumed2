/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#ifndef __PLUMED_core_ActionAtomistic_h
#define __PLUMED_core_ActionAtomistic_h

#include "Action.h"
#include "tools/Tensor.h"
#include "Atoms.h"
#include "tools/Pbc.h"
#include "tools/ForwardDecl.h"
#include <vector>
#include <set>
#include <map>

namespace PLMD {

class Pbc;
class PDB;
class ActionWithVirtualAtom;

/// \ingroup MULTIINHERIT
/// Action used to create objects that access the positions of the atoms from the MD code
class ActionAtomistic :
  virtual public Action
{
  std::vector<AtomNumber> indexes;         // the set of needed atoms
/// unique should be an ordered set since we later create a vector containing the corresponding indexes
  std::set<AtomNumber>  unique;
  std::vector<Vector>   positions;       // positions of the needed atoms
  double                energy;
  Value*                boxValue;
  ForwardDecl<Pbc>      pbc_fwd;
  Pbc&                  pbc=*pbc_fwd;
  Tensor                virial;
  std::vector<double>   masses;
  std::vector<double>   charges;


  std::string           extraCV;

  bool                  lockRequestAtoms; // forbid changes to request atoms

  bool                  donotretrieve;
  bool                  donotforce;

protected:
  bool                  chargesWereSet;
  Atoms&                atoms;

public:
/// Request an array of atoms.
/// This method is used to ask for a list of atoms. Atoms
/// should be asked for by number. If this routine is called
/// during the simulation, atoms will be available at the next step
/// MAYBE WE HAVE TO FIND SOMETHING MORE CLEAR FOR DYNAMIC
/// LISTS OF ATOMS
  void requestAtoms(const std::vector<AtomNumber> & a, const bool clearDep=true);
/// Get position of i-th atom (access by relative index)
  const Vector & getPosition(int)const;
/// Get position of i-th atom (access by absolute AtomNumber).
/// With direct access to the global atom array.
/// \warning Should be only used by actions that need to read the shared position array.
///          This array is insensitive to local changes such as makeWhole(), numerical derivatives, etc.
  Vector getGlobalPosition(AtomNumber)const;
/// Modify position of i-th atom (access by absolute AtomNumber).
/// \warning Should be only used by actions that need to modify the shared position array.
///          This array is insensitive to local changes such as makeWhole(), numerical derivatives, etc.
  void setGlobalPosition(AtomNumber, const Vector& pos);
/// Get total number of atoms, including virtual ones.
/// Can be used to make a loop on modifyGlobalPosition or getGlobalPosition.
  unsigned getTotAtoms()const;
/// Get modifiable force of i-th atom (access by absolute AtomNumber).
/// \warning  Should be used by action that need to modify the stored atomic forces.
///           This should be used with great care since the plumed core does
///           not usually keep all these forces up to date. In particular,
///           if an action require this, one should during constructor
///           call allowToAccessGlobalForces().
///           Notice that for efficiency reason plumed does not check if this is done!
//  Vector & modifyGlobalForce(AtomNumber);
/// Get modifiable virial
  Tensor & modifyVirial();
/// Get the value of the virial
  const Tensor & getVirial() const ;
/// Should be used by action that need to modify the stored virial
//  Tensor & modifyGlobalVirial();
/// Get modifiable PBC
/// Get box shape
  const Tensor & getBox()const;
/// Get the array of all positions
  const std::vector<Vector> & getPositions()const;
/// Get mass of i-th atom
  double getMass(int i)const;
/// Get charge of i-th atom
  double getCharge(int i)const;
/// Add force to an atom
  void addForce( const unsigned& i, const Vector& f );
/// Get a reference to virial array
  void addVirial( const Tensor& v );
/// Get number of available atoms
  unsigned getNumberOfAtoms()const {return indexes.size();}
/// Compute the pbc distance between two positions
  Vector pbcDistance(const Vector&,const Vector&)const;
/// Applies  PBCs to a seriens of positions or distances
  void pbcApply(std::vector<Vector>& dlist, unsigned max_index=0) const;
/// Get the vector of absolute indexes
  virtual const std::vector<AtomNumber> & getAbsoluteIndexes()const;
/// Get the absolute index of an atom
  AtomNumber getAbsoluteIndex(int i)const;
/// Parse a list of atoms without a numbered keyword
  void parseAtomList(const std::string&key,std::vector<AtomNumber> &t);
/// Parse an list of atom with a numbred keyword
  void parseAtomList(const std::string&key,const int num, std::vector<AtomNumber> &t);
/// Convert a set of read in strings into an atom list (this is used in parseAtomList)
  void interpretAtomList( std::vector<std::string>& strings, std::vector<AtomNumber> &t);
/// Change the box shape
  void changeBox( const Tensor& newbox );
/// Get reference to Pbc
  const Pbc & getPbc() const;
/// Add the forces to the atoms
  void setForcesOnAtoms( const std::vector<double>& forcesToApply, unsigned& ind );
/// Skip atom retrieval - use with care.
/// If this function is called during initialization, then atoms are
/// not going to be retrieved. Can be used for optimization. Notice that
/// calling getPosition(int) in an Action where DoNotRetrieve() was called might
/// lead to undefined behavior.
  void doNotRetrieve() {donotretrieve=true;}
/// Skip atom forces - use with care.
/// If this function is called during initialization, then forces are
/// not going to be propagated. Can be used for optimization.
  void doNotForce() {donotforce=true;}
/// Make atoms whole, assuming they are in the proper order
/// If start and end are set equal to 0 then do all atoms.
/// If start and end are set to particular indices then do a subset
  void makeWhole( const unsigned start=0, const unsigned end=0 );
/// Allow calls to modifyGlobalForce()
//   void allowToAccessGlobalForces() {atoms.zeroallforces=true;}
public:

// virtual functions:

  explicit ActionAtomistic(const ActionOptions&ao);
  ~ActionAtomistic();

  static void registerKeywords( Keywords& keys );

/// N.B. only pass an ActionWithValue to this routine if you know exactly what you
/// are doing.  The default will be correct for the vast majority of cases
  void   calculateNumericalDerivatives( ActionWithValue* a=NULL ) override;
/// Numerical derivative routine to use when using Actions that inherit from BOTH
/// ActionWithArguments and ActionAtomistic
  void calculateAtomicNumericalDerivatives( ActionWithValue* a, const unsigned& startnum );

  virtual void retrieveAtoms();
  void lockRequests() override;
  void unlockRequests() override;
  const std::set<AtomNumber> & getUnique()const;
/// Read in an input file containing atom positions and calculate the action for the atomic
/// configuration therin
  void readAtomsFromPDB( const PDB& pdb );
};

inline
const Vector & ActionAtomistic::getPosition(int i)const {
  return positions[i];
}


// inline
// Vector & ActionAtomistic::modifyGlobalForce(AtomNumber i) {
//   return atoms.forces[i.index()];
// }

// inline
// Tensor & ActionAtomistic::modifyGlobalVirial() {
//   return atoms.virial;
// }

inline
double ActionAtomistic::getMass(int i)const {
  plumed_assert( i<masses.size() );
  return masses[i];
}

inline
double ActionAtomistic::getCharge(int i) const {
  if( !chargesWereSet ) error("charges were not passed to plumed");
  return charges[i];
}

inline
const std::vector<AtomNumber> & ActionAtomistic::getAbsoluteIndexes()const {
  return indexes;
}

inline
AtomNumber ActionAtomistic::getAbsoluteIndex(int i)const {
  return indexes[i];
}

inline
const std::vector<Vector> & ActionAtomistic::getPositions()const {
  return positions;
}

inline
const Tensor & ActionAtomistic::getBox()const {
  return pbc.getBox();
}

inline
const Tensor & ActionAtomistic::getVirial() const {
  return virial;
}

inline
const Pbc & ActionAtomistic::getPbc() const {
  return pbc;
}

inline
void ActionAtomistic::lockRequests() {
  lockRequestAtoms=true;
}

inline
void ActionAtomistic::unlockRequests() {
  lockRequestAtoms=false;
}

inline
const std::set<AtomNumber> & ActionAtomistic::getUnique()const {
  return unique;
}

}

#endif
