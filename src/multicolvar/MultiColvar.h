/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
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

#include "core/ActionAtomistic.h"
#include "core/ActionWithValue.h"
#include "tools/DynamicList.h"
#include "vesselbase/ActionWithVessel.h"
#include "StoreCentralAtomsVessel.h"
#include <vector>

#define PLUMED_MULTICOLVAR_INIT(ao) Action(ao),MultiColvar(ao)

namespace PLMD {
namespace multicolvar {

/**
\ingroup INHERIT
This is the abstract base class to use for creating distributions of colvars and functions
thereof, whtin it there is \ref AddingAMultiColvar "information" as to how to go implementing these types of actions.
*/

class MultiColvar :
  public ActionAtomistic,
  public ActionWithValue,
  public vesselbase::ActionWithVessel
  {
friend class Region;
friend class StoreCentralAtomsVessel;
private:
  bool usepbc;
  bool readatoms;
  bool verbose_output;
/// Everything for controlling the updating of neighbor lists
  int updateFreq;
  unsigned lastUpdate;
  bool reduceAtNextStep;
/// A flag that tells us the position has been set already
  bool posHasBeenSet;
/// The list of all the atoms involved in the colvar
  DynamicList<AtomNumber> all_atoms;
/// A DynamicList containing the numbers of 1 - ncolvars
  DynamicList<unsigned> colvar_list; 
/// The lists of the atoms involved in each of the individual colvars
/// note these refer to the atoms in all_atoms
  std::vector< DynamicList<unsigned> > colvar_atoms;
/// These are used to store the values of CVs etc so they can be retrieved by distribution
/// functions
  std::vector<Vector> pos;
  Tensor ibox;
  bool centralAtomDerivativesAreInFractional;
  std::vector<Tensor> central_derivs;
/// Read in ATOMS keyword
  void readAtomsKeyword( int& natoms );
/// Read in the various GROUP keywords
  void readGroupsKeyword( int& natoms );
/// Read in the various SPECIES keywords
  void readSpeciesKeyword( int& natoms );
/// Update the atoms request
  void requestAtoms();
protected:
/// Are we on an update step
  bool updatetime;
/// Read in all the keywords that can be used to define atoms
  void readAtoms( int& natoms );
/// Read in the atoms that form the backbone of a polymeric chain
  void readBackboneAtoms( const std::vector<std::string>& backnames, std::vector<unsigned>& chain_lengths );
/// Add a colvar to the set of colvars we are calculating (in practise just a list of atoms)
  void addColvar( const std::vector<unsigned>& newatoms );
/// Get the separation between a pair of vectors
  Vector getSeparation( const Vector& vec1, const Vector& vec2 ) const ;
/// Find out if it is time to do neighbor list update
  bool isTimeForNeighborListUpdate() const ;
/// Update the list of atoms after the neighbor list step
  void removeAtomRequest( const unsigned& aa );
/// Do we use pbc to calculate this quantity
  bool usesPbc() const ;
/// Add some derivatives for an atom 
  void addAtomsDerivatives(const int&,const Vector&);
/// Add some derivatives to the virial
  void addBoxDerivatives(const Tensor&);
/// Add derivative of central atom position wrt to position of iatom'th atom
  void addCentralAtomDerivatives( const unsigned& iatom, const Tensor& der );
/// Retrieve derivative of central atom position wrt jcomp'th component of position of iatom'th atom
  double getCentralAtomDerivative( const unsigned& iatom, const unsigned jcomp, const Vector& df ) const ;
/// Set a weight for this colvar (used in MEAN and HISTOGRAM)
  void setWeight( const double& weight );
/// Set the derivative of the weight (used in MEAN and HISTOGRAM)
  void addAtomsDerivativeOfWeight( const unsigned& i, const Vector& wder );
  void addBoxDerivativesOfWeight( const Tensor& vir );
/// Get the number of atoms in this particular colvar
  unsigned getNAtoms() const;
/// Get all the positions
  const std::vector<Vector> & getPositions();
/// This can be used to get rid of erroenous effects that might happen
/// because molecules are split by the pbcs.
  void setAlignedPositions( const std::vector<Vector>& );
public:
  MultiColvar(const ActionOptions&);
  ~MultiColvar(){};
  static void registerKeywords( Keywords& keys );
/// Calculate the multicolvar
  void calculate();
/// Prepare for the calculation
  virtual void prepare();
/// Apply the forces on the values
  void apply();
/// Get the number of derivatives for this action
  unsigned getNumberOfDerivatives();  // N.B. This is replacing the virtual function in ActionWithValue
/// Return the number of Colvars this is calculating
  unsigned getNumberOfFunctionsInAction();  
/// Return the number of derivatives for a given colvar
  unsigned getNumberOfDerivatives( const unsigned& j );
/// Retrieve the position of the central atom
  Vector retrieveCentralAtomPos( const bool& frac );
/// Make sure we calculate the position of the central atom
  void useCentralAtom();
/// Merge the derivatives 
  void chainRuleForElementDerivatives( const unsigned& , const unsigned& , const unsigned& , const unsigned& , const double& , vesselbase::Vessel* );
/// Also used for derivative merging
  unsigned getOutputDerivativeIndex( const unsigned& ival, const unsigned& i );
/// Can we skip the calculations of quantities
  virtual bool isPossibleToSkip();
/// Turn of atom requests when this colvar is deactivated cos its small
  void deactivate_task();
/// Turn on atom requests when the colvar is activated
  void activateValue( const unsigned j );
/// Calcualte the colvar
  bool performTask( const unsigned& j );
/// You can use this to screen contributions that are very small so we can avoid expensive (and pointless) calculations
  virtual bool contributionIsSmall(){ plumed_dbg_assert( !isPossibleToSkip() ); return false; }
/// And a virtual function which actually computes the colvar
  virtual double compute( const unsigned& j )=0;  
/// A virtual routine to get the position of the central atom - used for things like cv gradient
  virtual Vector getCentralAtom(); 
/// Is this a density?
  virtual bool isDensity(){ return false; }
/// Get the position of atom iatom
  const Vector & getPosition(unsigned) const;
/// Get the mass of atom iatom
  double getMass(unsigned) const ;
/// Get the charge of atom iatom
  double getCharge(unsigned) const ;
/// Get the absolute index of atom iatom
  AtomNumber getAbsoluteIndex(unsigned) const ;
/// Return a pointer to the vessel that stores the positions of 
/// all the central atoms
  StoreCentralAtomsVessel* getCentralAtoms();
};

inline
bool MultiColvar::isPossibleToSkip(){
  return (colvar_atoms[current].getNumberActive()==0); 
}

inline
unsigned MultiColvar::getNumberOfDerivatives(){
  return 3*getNumberOfAtoms()+9;
} 

inline
unsigned MultiColvar::getNumberOfFunctionsInAction(){
  return colvar_atoms.size();
}

inline
void MultiColvar::deactivate_task(){
  if( !reduceAtNextStep ) return;          // Deactivating tasks only possible during neighbor list update
  colvar_list.deactivate(current);         // Deactivate the colvar from the list
  colvar_atoms[current].deactivateAll();   // Deactivate all atom requests for this colvar
}

inline
void MultiColvar::activateValue( const unsigned j ){
  colvar_atoms[j].activateAll(); 
  colvar_atoms[j].updateActiveMembers();
}

inline
unsigned MultiColvar::getNumberOfDerivatives( const unsigned& j ){
  return 3*colvar_atoms[j].getNumberActive() + 9;
}

inline
bool MultiColvar::isTimeForNeighborListUpdate() const {
  return reduceAtNextStep;
}

inline
void MultiColvar::removeAtomRequest( const unsigned& i ){
  plumed_massert(reduceAtNextStep,"found removeAtomRequest but not during neighbor list step");
  colvar_atoms[current].deactivate(i); 
}

inline
bool MultiColvar::usesPbc() const {
  return usepbc;
}

inline
unsigned MultiColvar::getNAtoms() const {
  return colvar_atoms[current].getNumberActive();
}

inline
const Vector & MultiColvar::getPosition( unsigned iatom ) const {
  plumed_dbg_assert( iatom<colvar_atoms[current].getNumberActive() );
  return ActionAtomistic::getPosition( colvar_atoms[current][iatom] );
}

inline
double MultiColvar::getMass(unsigned iatom ) const {
  plumed_dbg_assert( iatom<colvar_atoms[current].getNumberActive() );
  return ActionAtomistic::getMass( colvar_atoms[current][iatom] );
}

inline
double MultiColvar::getCharge(unsigned iatom ) const {
  plumed_dbg_assert( iatom<colvar_atoms[current].getNumberActive() );
  return ActionAtomistic::getCharge( colvar_atoms[current][iatom] );
}

inline
AtomNumber MultiColvar::getAbsoluteIndex(unsigned iatom) const {
  plumed_dbg_assert( iatom<colvar_atoms[current].getNumberActive() );
  return ActionAtomistic::getAbsoluteIndex( colvar_atoms[current][iatom] );
}

inline
void MultiColvar::addAtomsDerivatives(const int& iatom, const Vector& der){
  plumed_dbg_assert( iatom<colvar_atoms[current].getNumberActive() );
  addElementDerivative( 3*iatom+0, der[0] );
  addElementDerivative( 3*iatom+1, der[1] );
  addElementDerivative( 3*iatom+2, der[2] );
}

inline
void MultiColvar::addBoxDerivatives(const Tensor& vir){
  int natoms=colvar_atoms[current].getNumberActive();
  addElementDerivative( 3*natoms+0, vir(0,0) );
  addElementDerivative( 3*natoms+1, vir(0,1) );
  addElementDerivative( 3*natoms+2, vir(0,2) );
  addElementDerivative( 3*natoms+3, vir(1,0) );
  addElementDerivative( 3*natoms+4, vir(1,1) );
  addElementDerivative( 3*natoms+5, vir(1,2) );
  addElementDerivative( 3*natoms+6, vir(2,0) );
  addElementDerivative( 3*natoms+7, vir(2,1) );
  addElementDerivative( 3*natoms+8, vir(2,2) );
}

inline
unsigned MultiColvar::getOutputDerivativeIndex( const unsigned& ival, const unsigned& i ){
  if( i<3*getNumberOfAtoms() ){
      unsigned inat=std::floor( i/3 );
      unsigned thisatom=linkIndex( inat, colvar_atoms[ival], all_atoms );
      return 3*thisatom + i%3; 
  } 
  return 3*getNumberOfAtoms() + i - 3*colvar_atoms[ival].getNumberActive(); 
}

inline
void MultiColvar::setWeight( const double& weight ){
  setElementValue( 1, weight );
}

inline
void MultiColvar::addAtomsDerivativeOfWeight( const unsigned& iatom, const Vector& wder ){
  plumed_dbg_assert( iatom<colvar_atoms[current].getNumberActive() );
  int nstart  = 3*getNumberOfAtoms() + 9 + 3*iatom;
  addElementDerivative( nstart + 0, wder[0] );
  addElementDerivative( nstart + 1, wder[1] );
  addElementDerivative( nstart + 2, wder[2] );
}

inline
void MultiColvar::addBoxDerivativesOfWeight( const Tensor& vir ){
  int nstart = 3*getNumberOfAtoms() + 9 + 3*colvar_atoms[current].getNumberActive();
  addElementDerivative( nstart+0, vir(0,0) );
  addElementDerivative( nstart+1, vir(0,1) );
  addElementDerivative( nstart+2, vir(0,2) );
  addElementDerivative( nstart+3, vir(1,0) );
  addElementDerivative( nstart+4, vir(1,1) );
  addElementDerivative( nstart+5, vir(1,2) );
  addElementDerivative( nstart+6, vir(2,0) );
  addElementDerivative( nstart+7, vir(2,1) );
  addElementDerivative( nstart+8, vir(2,2) );
}


}
}

#endif
