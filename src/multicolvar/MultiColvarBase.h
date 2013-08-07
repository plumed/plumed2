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
#ifndef __PLUMED_multicolvar_MultiColvarBase_h
#define __PLUMED_multicolvar_MultiColvarBase_h

#include "core/ActionAtomistic.h"
#include "core/ActionWithValue.h"
#include "tools/DynamicList.h"
#include "vesselbase/ActionWithVessel.h"
#include "StoreCentralAtomsVessel.h"
#include <vector>

namespace PLMD {
namespace multicolvar {

class MultiColvarBase :
  public ActionAtomistic,
  public ActionWithValue,
  public vesselbase::ActionWithVessel
  {
friend class ActionVolume;
friend class VolumeSubcell;
friend class StoreColvarVessel;
friend class StoreCentralAtomsVessel;
friend class MultiColvarFunction;
friend class MultiColvar;
private:
/// Use periodic boundary conditions
  bool usepbc;
/// Everything for controlling the updating of neighbor lists
  int updateFreq;
  unsigned lastUpdate;
/// The list of all the atoms involved in the colvar
  DynamicList<AtomNumber> all_atoms;
/// Variables used for central atoms
  Tensor ibox;
  DynamicList<unsigned> atomsWithCatomDer;
/// The forces we are going to apply to things
  std::vector<double> forcesToApply;
/// This resizes the local arrays after neighbor list updates and during initialization
  void resizeLocalArrays();
protected:
/// A dynamic list containing those atoms with derivatives
  DynamicList<unsigned> atoms_with_derivatives;
/// The lists of the atoms involved in each of the individual colvars
  std::vector< DynamicList<unsigned> > colvar_atoms;
/// Add a colvar to the set of colvars we are calculating (in practise just a list of atoms)
  void addColvar( const std::vector<unsigned>& newatoms );
/// Finish setting up the multicolvar base
  void setupMultiColvarBase();
/// Request the atoms from action atomistic
  void requestAtoms();
/// Get the separation between a pair of vectors
  Vector getSeparation( const Vector& vec1, const Vector& vec2 ) const ;
/// Do we use pbc to calculate this quantity
  bool usesPbc() const ;
/// Add some derivatives for an atom
  void addAtomsDerivatives(const unsigned&, const unsigned&, const Vector& );
/// Add some derivatives for a box
  void addBoxDerivatives(const unsigned&, const Tensor& );
/// Add some derivatives of the value to the virial
  void addBoxDerivatives(const Tensor&);
/// Retrieve derivative of central atom position wrt jcomp'th component of position of iatom'th atom
  double getCentralAtomDerivative( const unsigned& iatom, const unsigned& jcomp, const Vector& df );
/// Set a weight for this colvar (used in MEAN and HISTOGRAM)
  void setWeight( const double& weight );
/// Set the derivative of the weight (used in MEAN and HISTOGRAM)
  void addBoxDerivativesOfWeight( const Tensor& vir );
/// Get the number of atoms in this particular colvar
  unsigned getNAtoms() const;
/// Update the list of atoms after the neighbor list step
  void removeAtomRequest( const unsigned& aa, const double& weight );
/// Add derivative of central atom position wrt to position of iatom'th atom
  void addCentralAtomDerivatives( const unsigned& iatom, const Tensor& der );
/// Get the indices for the central atom
  void getCentralAtomIndexList( const unsigned& ntotal, const unsigned& jstore, const unsigned& maxder, std::vector<unsigned>& indices ) const ;
/// Calculate and store getElementValue(uder)/getElementValue(vder) and its derivatives in getElementValue(iout)
  void quotientRule( const unsigned& uder, const unsigned& vder, const unsigned& iout );
public:
  MultiColvarBase(const ActionOptions&);
  ~MultiColvarBase(){}
  static void registerKeywords( Keywords& keys );
/// Prepare for the calculation
  void prepare();
  virtual void resizeDynamicArrays()=0;
/// Perform one of the tasks
  void performTask( const unsigned& j );
/// And a virtual function which actually computes the colvar
  virtual double doCalculation( const unsigned& j );  
/// Update the atoms that have derivatives
  virtual void updateActiveAtoms()=0;
/// This is replaced once we have a function to calculate the cv
  virtual double compute( const unsigned& j )=0;
/// These replace the functions in ActionWithVessel to make the code faster
  void mergeDerivatives( const unsigned& ider, const double& df );
  void clearDerivativesAfterTask( const unsigned& ider );
/// Apply the forces from this action
  void apply();
/// Deactivate one of the tasks
  void deactivate_task();
/// Get the number of derivatives for this action
  unsigned getNumberOfDerivatives();  // N.B. This is replacing the virtual function in ActionWithValue
/// Get the number of quantities that are calculated each time
  virtual unsigned getNumberOfQuantities();
/// Retrieve the position of the central atom
  Vector retrieveCentralAtomPos();
/// You can use this to screen contributions that are very small so we can avoid expensive (and pointless) calculations
  virtual void calculateWeight();
/// A virtual routine to get the position of the central atom - used for things like cv gradient
  virtual Vector calculateCentralAtomPosition()=0; 
/// Get the list of indices that have derivatives
 void getIndexList( const unsigned& ntotal, const unsigned& jstore, const unsigned& maxder, std::vector<unsigned>& indices );
/// Is this a density?
  virtual bool isDensity(){ return false; }
/// Return a pointer to the vessel that stores the positions of 
/// all the central atoms
  StoreCentralAtomsVessel* getCentralAtoms();
/// Copy the list of atoms involved to a second MultiColvarBase (used by functions)
  void copyAtomListToFunction( MultiColvarBase* myfunction );
/// Return the number of the colvar in which iatom is the first atom
  unsigned getInternalIndex( const AtomNumber& iatom ) const ;
/// Make sure the same list of atoms is active in a function
  void copyActiveAtomsToFunction( MultiColvarBase* myfunction );
/// Activate the atoms that have derivatives from a storeDataVessel
  void activateIndexes( const unsigned& istart, const unsigned& number, const std::vector<unsigned>& indexes ); 
  void activateIndex( const unsigned& );
};

inline
unsigned MultiColvarBase::getNumberOfDerivatives(){
  return 3*getNumberOfAtoms()+9;
}

inline
void MultiColvarBase::removeAtomRequest( const unsigned& i, const double& weight ){
  if( !contributorsAreUnlocked ) return;
  plumed_dbg_assert( weight<getTolerance() );
  if( weight<getNLTolerance() ) colvar_atoms[current].deactivate( i );
}

inline
void MultiColvarBase::deactivate_task(){
  if( !contributorsAreUnlocked ) return;      // Deactivating tasks only possible during neighbor list update
  colvar_atoms[current].deactivateAll();      // Deactivate all atom requests for this colvar
  ActionWithVessel::deactivate_task();        // Deactivate the colvar from the list
}

inline
bool MultiColvarBase::usesPbc() const {
  return usepbc;
}

inline
unsigned MultiColvarBase::getNumberOfQuantities(){
  return 5;
}

inline
unsigned MultiColvarBase::getNAtoms() const {
  return colvar_atoms[current].getNumberActive();
}

inline
void MultiColvarBase::addAtomsDerivatives(const unsigned& ielem, const unsigned& iatom, const Vector& der ){
  atoms_with_derivatives.activate(iatom);
  unsigned ibase=ielem*getNumberOfDerivatives() + 3*iatom;
  for(unsigned i=0;i<3;++i) addElementDerivative( ibase + i, der[i] );
}

inline 
void MultiColvarBase::addBoxDerivatives(const unsigned& ielem, const Tensor& vir ){
  unsigned ibase=ielem*getNumberOfDerivatives() + 3*getNumberOfAtoms();
  for(unsigned i=0;i<3;++i) for(unsigned j=0;j<3;++j) addElementDerivative( ibase+3*i+j, vir(i,j) );
}

inline
void MultiColvarBase::addBoxDerivatives(const Tensor& vir){
  addBoxDerivatives( 0, vir );
}

inline
void MultiColvarBase::calculateWeight(){
  setElementValue( 1, 1.0 );
}

inline
void MultiColvarBase::setWeight( const double& weight ){
  setElementValue( 1, weight );
}

inline
void MultiColvarBase::addBoxDerivativesOfWeight( const Tensor& vir ){
  addBoxDerivatives( 1, vir );
}

inline
void MultiColvarBase::activateIndex( const unsigned& ider ){
  unsigned iatom = std::floor( ider / 3 );
  atoms_with_derivatives.activate( iatom );
}

}
}

#endif
