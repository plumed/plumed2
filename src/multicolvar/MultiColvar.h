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
friend class ActionVolume;
friend class MultiColvarFunction;
friend class StoreCentralAtomsVessel;
private:
  bool usepbc;
/// Have atoms been read in
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
/// A dynamic list containing those atoms with derivatives
  DynamicList<unsigned> atoms_with_derivatives;
/// The lists of the atoms involved in each of the individual colvars
/// note these refer to the atoms in all_atoms
  std::vector< DynamicList<unsigned> > colvar_atoms;
/// These are used to store the values of CVs etc so they can be retrieved by distribution
/// functions
  std::vector<Vector> pos;
  Tensor ibox;
  bool centralAtomDerivativesAreInFractional;
  DynamicList<unsigned> atomsWithCatomDer;
  std::vector<Tensor> central_derivs;
/// Stuff used to merge derivatives
  unsigned imerge_deriv, imerge_natoms;
/// The forces we are going to apply to things
  std::vector<double> forcesToApply;
/// Read in the various GROUP keywords
  void readGroupsKeyword( int& natoms );
/// Read in the various SPECIES keywords
  void readSpeciesKeyword( int& natoms );
/// Update the atoms request
  void requestAtoms();
/// Resize all the dynamic arrays (used at neighbor list update time and during setup)
  void resizeDynamicArrays();
protected:
/// Are we on an update step
  bool updatetime;
/// Read in all the keywords that can be used to define atoms
  void readAtoms( int& natoms );
/// Read in ATOMS keyword
  void readAtomsLikeKeyword( const std::string key, int& natoms );
/// Read in the atoms that form the backbone of a polymeric chain
  void readBackboneAtoms( const std::vector<std::string>& backnames, std::vector<unsigned>& chain_lengths );
/// Add a colvar to the set of colvars we are calculating (in practise just a list of atoms)
  void addColvar( const std::vector<unsigned>& newatoms );
/// Return the index of an atom
  unsigned getAtomIndex( const unsigned& ) const ;
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
/// Retrieve the position of the central atom
  Vector retrieveCentralAtomPos( const bool& frac );
/// Turn of atom requests when this colvar is deactivated cos its small
  void deactivate_task();
/// Turn on atom requests when the colvar is activated
  void activateValue( const unsigned j );
/// Calcualte the colvar
  bool performTask( const unsigned& j );
/// You can use this to screen contributions that are very small so we can avoid expensive (and pointless) calculations
  virtual void calculateWeight();
/// And a virtual function which actually computes the colvar
  virtual double compute( const unsigned& j )=0;  
/// These are replacing the virtual methods in ActionWithVessel
  void mergeDerivatives( const unsigned& ider, const double& df );
  unsigned getFirstDerivativeToMerge();
  unsigned getNextDerivativeToMerge( const unsigned& );
/// Clear tempory data that is calculated for each task
  void clearDerivativesAfterTask( const unsigned& );
/// Build lists of indexes of the derivatives from the colvar atoms arrays
  void buildDerivativeIndexArrays( std::vector< DynamicList<unsigned> >& active_der );
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
unsigned MultiColvar::getFirstDerivativeToMerge(){
  imerge_deriv=0; imerge_natoms=atoms_with_derivatives.getNumberActive();
  return 3*getAtomIndex( atoms_with_derivatives[imerge_deriv] ); 
}

inline
unsigned MultiColvar::getNextDerivativeToMerge( const unsigned& j){
  imerge_deriv++; 
  if( imerge_deriv>=3*imerge_natoms ) return 3*getNumberOfAtoms() - 3*imerge_natoms + imerge_deriv;
  unsigned imerge_atom=std::floor( imerge_deriv / 3 );
  return 3*getAtomIndex( imerge_atom ) + imerge_deriv%3;
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
  deactivateCurrentTask();                 // Deactivate the colvar from the list
  colvar_atoms[current].deactivateAll();   // Deactivate all atom requests for this colvar
}

inline
bool MultiColvar::isTimeForNeighborListUpdate() const {
  return reduceAtNextStep;
}

inline
void MultiColvar::removeAtomRequest( const unsigned& i ){
  plumed_massert(reduceAtNextStep,"found removeAtomRequest but not during neighbor list step");
  colvar_atoms[current].deactivate( getAtomIndex(i) ); 
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
unsigned MultiColvar::getAtomIndex( const unsigned& iatom ) const {
  plumed_dbg_assert( iatom<colvar_atoms[current].getNumberActive() );
  return colvar_atoms[current][iatom];
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
  unsigned jatom=getAtomIndex(iatom); 
  atoms_with_derivatives.activate(iatom);
  addElementDerivative( 3*jatom+0, der[0] );
  addElementDerivative( 3*jatom+1, der[1] );
  addElementDerivative( 3*jatom+2, der[2] );
}

inline
void MultiColvar::addBoxDerivatives(const Tensor& vir){
  unsigned nstart=3*getNumberOfAtoms(); 
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

inline
void MultiColvar::calculateWeight(){
  setElementValue( 1, 1.0 );
}

inline
void MultiColvar::setWeight( const double& weight ){
  setElementValue( 1, weight );
}

inline
void MultiColvar::addAtomsDerivativeOfWeight( const unsigned& iatom, const Vector& wder ){
  unsigned nstart  = 3*getNumberOfAtoms() + 9 + 3*getAtomIndex(iatom);
  atoms_with_derivatives.activate(iatom); 
  addElementDerivative( nstart + 0, wder[0] );
  addElementDerivative( nstart + 1, wder[1] );
  addElementDerivative( nstart + 2, wder[2] );
}

inline
void MultiColvar::addBoxDerivativesOfWeight( const Tensor& vir ){
  int nstart = 6*getNumberOfAtoms() + 9;
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
