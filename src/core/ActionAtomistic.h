/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "tools/Pbc.h"
#include "tools/ForwardDecl.h"
#include "Value.h"
#include <vector>
#include <map>

namespace PLMD {

class Pbc;
class PDB;
class GenericMolInfo;
class Tree;

namespace colvar {
class SelectMassCharge;
}

/// \ingroup MULTIINHERIT
/// Action used to create objects that access the positions of the atoms from the MD code
class ActionAtomistic :
  virtual public Action {
  friend class Group;
  friend class DomainDecomposition;
  friend class colvar::SelectMassCharge;
  friend class ActionWithVirtualAtom;

  std::vector<AtomNumber> indexes;         // the set of needed atoms
  std::vector<std::size_t>   value_depends;   // The list of values that are being used
  std::vector<std::pair<std::size_t, std::size_t > > atom_value_ind;  // The list of values and indices for the atoms that are being used
  std::vector<std::pair<std::size_t,std::vector<std::size_t>>> atom_value_ind_grouped;
/// unique should be an ordered set since we later create a vector containing the corresponding indexes
  std::vector<AtomNumber>  unique;
/// unique_local should be an ordered set since we later create a vector containing the corresponding indexes
  bool unique_local_needs_update;
  std::vector<AtomNumber>  unique_local;
  std::vector<Vector>   actionPositions;       // positions of the needed atoms
  double                energy;
  Value*                boxValue;
  ForwardDecl<Pbc>      pbc_fwd;
  Pbc&                  actionPbc=*pbc_fwd;
  std::vector<double>   masses;
  std::vector<double>   charges;

  std::vector<Vector>   forces;          // forces on the needed atoms
  double                forceOnEnergy;

  double                forceOnExtraCV;

  bool                  lockRequestAtoms; // forbid changes to request atoms

  bool                  donotretrieve;
  bool                  donotforce;

  // EMST
  GenericMolInfo* actionMoldat{nullptr};
  std::unique_ptr<Tree> tree;

/// Values that hold information about atom positions and charges
  std::vector<Value*>   xpos, ypos, zpos, masv, chargev;
  void updateUniqueLocal( const bool& useunique, const std::vector<int>& g2l );
protected:
  bool                  massesWereSet;
  bool                  chargesWereSet;
  void setExtraCV(const std::string &name);
/// Used to interpret whether this index is a virtual atom or a real atom
  std::pair<std::size_t, std::size_t> getValueIndices( const AtomNumber& i ) const ;
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
  Vector getGlobalPosition(const std::pair<std::size_t,std::size_t>& ) const ;
/// Modify position of i-th atom (access by absolute AtomNumber).
/// \warning Should be only used by actions that need to modify the shared position array.
///          This array is insensitive to local changes such as makeWhole(), numerical derivatives, etc.
  void setGlobalPosition(const std::pair<std::size_t,std::size_t>&, const Vector& pos);
/// Get total number of atoms, including virtual ones.
/// Can be used to make a loop on modifyGlobalPosition or getGlobalPosition.
  unsigned getTotAtoms()const;
/// Get box shape
  const Tensor & getBox()const;
/// Get the array of all positions
  const std::vector<Vector> & getPositions()const;
/// Get the array of all masses
  const std::vector<double>& getMasses()const;
/// Get the array of all charges
  const std::vector<double>& getCharges( const bool allowempty=false )const;
/// Get the virial that is acting
  Tensor getVirial() const ;
/// Get energy
  const double & getEnergy()const;
/// Get mass of i-th atom
  double getMass(int i)const;
/// Get charge of i-th atom
  double getCharge(int i)const;
/// Get the force acting on a particular atom
  Vector getForce( const std::pair<std::size_t, std::size_t>& a ) const ;
/// Add force to an atom
  void addForce( const std::pair<std::size_t, std::size_t>& a, const Vector& f );
/// Get a reference to force on energy
  double & modifyForceOnEnergy();
/// Get number of available atoms
  unsigned getNumberOfAtoms()const {
    return indexes.size();
  }
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
/// Interpret the atom selection.  Just a wrapper to the static function with four arguments called interpretAtomList that passes xpos and this.
  void interpretAtomList( std::vector<std::string>& strings, std::vector<AtomNumber> &t);
/// Convert a set of read in strings into an atom list (this is used in parseAtomList)
  static void interpretAtomList( std::vector<std::string>& strings, const std::vector<Value*>& xpos, Action* action, std::vector<AtomNumber> &t);
/// This gets std::vector that contain the PLMD::Value objects that contain xpositions, ypositions, zpositions, masses and charges
  static void getAtomValuesFromPlumedObject( const PlumedMain& plumed, std::vector<Value*>& xpos, std::vector<Value*>& ypos, std::vector<Value*>& zpos, std::vector<Value*>& masv, std::vector<Value*>& chargev );
/// Change the box shape
  void changeBox( const Tensor& newbox );
/// Get reference to Pbc
  const Pbc & getPbc() const;
/// Add the forces to the atoms
  void setForcesOnAtoms( const std::vector<double>& forcesToApply, unsigned& ind );
/// Add the virial forces
  void setForcesOnCell(const std::vector<double>& forcesToApply, unsigned& ind);
/// Add the virial forces (span-like syntax)
  void setForcesOnCell(const double* forcesToApply, std::size_t size, unsigned& ind);
/// Skip atom retrieval - use with care.
/// If this function is called during initialization, then atoms are
/// not going to be retrieved. Can be used for optimization. Notice that
/// calling getPosition(int) in an Action where DoNotRetrieve() was called might
/// lead to undefined behavior.
  void doNotRetrieve() {
    donotretrieve=true;
  }
/// Skip atom forces - use with care.
/// If this function is called during initialization, then forces are
/// not going to be propagated. Can be used for optimization.
  void doNotForce() {
    donotforce=true;
  }
/// Make atoms whole, assuming they are in the proper order
  void makeWhole();
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

  virtual void retrieveAtoms( const bool& force=false );
  void lockRequests() override;
  void unlockRequests() override;
  const std::vector<AtomNumber> & getUnique()const;
  const std::vector<AtomNumber> & getUniqueLocal()const;
/// Read in an input file containing atom positions and calculate the action for the atomic
/// configuration therin
  void readAtomsFromPDB( const PDB& pdb ) override;
/// Transfer the gradients
  void getGradient( const unsigned& ind, Vector& deriv, std::map<AtomNumber,Vector>& gradients ) const ;
  ActionAtomistic* castToActionAtomistic() noexcept final {
    return this;
  }
  virtual bool actionHasForces();
};

inline
const Vector & ActionAtomistic::getPosition(int i)const {
  return actionPositions[i];
}

inline
double ActionAtomistic::getMass(int i)const {
  if( !massesWereSet ) {
    log.printf("WARNING: masses were not passed to plumed\n");
  }
  return masses[i];
}

inline
double ActionAtomistic::getCharge(int i) const {
  if( !chargesWereSet ) {
    error("charges were not passed to plumed");
  }
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
  return actionPositions;
}

inline
const std::vector<double> & ActionAtomistic::getMasses()const {
  return masses;
}

inline
const std::vector<double> & ActionAtomistic::getCharges( const bool allowempty )const {
  if( !allowempty && !chargesWereSet ) {
    error("charges were not passed to plumed");
  }
  return charges;
}

inline
const double & ActionAtomistic::getEnergy()const {
  return energy;
}

inline
const Tensor & ActionAtomistic::getBox()const {
  return actionPbc.getBox();
}

inline
double & ActionAtomistic::modifyForceOnEnergy() {
  return forceOnEnergy;
}

inline
const Pbc & ActionAtomistic::getPbc() const {
  return actionPbc;
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
const std::vector<AtomNumber> & ActionAtomistic::getUnique()const {
  return unique;
}

inline
const std::vector<AtomNumber> & ActionAtomistic::getUniqueLocal() const {
  return unique_local;
}

inline
Vector ActionAtomistic::getGlobalPosition(const std::pair<std::size_t,std::size_t>& a) const {
  Vector pos;
  pos[0]=xpos[a.first]->data[a.second];
  pos[1]=ypos[a.first]->data[a.second];
  pos[2]=zpos[a.first]->data[a.second];
  return pos;
}

inline
void ActionAtomistic::setGlobalPosition(const std::pair<std::size_t, std::size_t>& a, const Vector& pos ) {
  xpos[a.first]->data[a.second]=pos[0];
  ypos[a.first]->data[a.second]=pos[1];
  zpos[a.first]->data[a.second]=pos[2];
}

inline
Vector ActionAtomistic::getForce( const std::pair<std::size_t, std::size_t>& a ) const {
  Vector f;
  f[0]=xpos[a.first]->getForce(a.second);
  f[1]=ypos[a.first]->getForce(a.second);
  f[2]=zpos[a.first]->getForce(a.second);
  return f;
}

inline
void ActionAtomistic::addForce( const std::pair<std::size_t, std::size_t>& a, const Vector& f ) {
  xpos[a.first]->addForce( a.second, f[0] );
  ypos[a.first]->addForce( a.second, f[1] );
  zpos[a.first]->addForce( a.second, f[2] );
}

inline
Vector ActionAtomistic::pbcDistance(const Vector &v1,const Vector &v2)const {
  return actionPbc.distance(v1,v2);
}

}
#endif
