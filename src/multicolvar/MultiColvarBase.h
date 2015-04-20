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
#ifndef __PLUMED_multicolvar_MultiColvarBase_h
#define __PLUMED_multicolvar_MultiColvarBase_h

#include "core/ActionAtomistic.h"
#include "core/ActionWithValue.h"
#include "tools/DynamicList.h"
#include "tools/LinkCells.h"
#include "vesselbase/ActionWithVessel.h"
#include "StoreColvarVessel.h"
#include "StoreCentralAtomsVessel.h"
#include <vector>

namespace PLMD {
namespace multicolvar {

class BridgedMultiColvarFunction;

class MultiColvarBase :
  public ActionAtomistic,
  public ActionWithValue,
  public vesselbase::ActionWithVessel
  {
friend class StoreCentralAtomsVessel;
friend class MultiColvarFunction;
friend class BridgedMultiColvarFunction;
friend class VolumeGradientBase;
friend class MultiColvarFilter;
friend class MultiColvar;
private:
/// Use periodic boundary conditions
  bool usepbc;
/// Variables used for central atoms
  Tensor ibox;
  DynamicList<unsigned> atomsWithCatomDer;
/// The forces we are going to apply to things
  std::vector<double> forcesToApply;
/// Stuff for link cells - this is used to make coordination number like variables faster
  LinkCells linkcells;
/// This remembers where the boundaries are for the tasks. It makes link cells work fast
  Matrix<std::pair<unsigned,unsigned> > bookeeping;
/// A copy of the vessel containing the catoms
  StoreCentralAtomsVessel* mycatoms;
/// A copy of the vessel containg the values of each colvar
  StoreColvarVessel* myvalues;
/// This resizes the local arrays after neighbor list updates and during initialization
  void resizeLocalArrays();
/// This resizes the arrays that are used for link cell update
  void resizeBookeepingArray( const unsigned& num1, const unsigned& num2 );
protected:
/// A dynamic list containing those atoms with derivatives
  DynamicList<unsigned> atoms_with_derivatives;
/// Using the species keyword to read in atoms
  bool usespecies;
/// Number of atoms in each block
  unsigned nblock;
/// This is used when turning cvcodes into atom numbers
  std::vector<unsigned> decoder;
/// Blocks of atom numbers
  std::vector< std::vector<unsigned> > ablocks;
/// Number of atoms in the cv - set at start of calculation
  unsigned natomsper;  
/// Vector containing the indices of the current atoms
  std::vector<unsigned> current_atoms;
/// Add a task to the list of tasks
  void addTaskToList( const unsigned& taskCode );
/// Finish setting up the multicolvar base
  void setupMultiColvarBase();
/// Set the value of the cutoff for the link cells
  void setLinkCellCutoff( const double& lcut );
/// Setup link cells in order to make this calculation faster
  void setupLinkCells();
/// Get the separation between a pair of vectors
  Vector getSeparation( const Vector& vec1, const Vector& vec2 ) const ;
/// Apply PBCs over a set of distance vectors
  void applyPbc(std::vector<Vector>& dlist, unsigned max_index=0) const;
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
/// Add derivative of central atom position wrt to position of iatom'th atom
  void addCentralAtomDerivatives( const unsigned& iatom, const Tensor& der );
/// Get the indices for the central atom
  void getCentralAtomIndexList( const unsigned& ntotal, const unsigned& jstore, const unsigned& maxder, std::vector<unsigned>& indices ) const ;
/// This sets up the list of atoms that are involved in this colvar
  bool setupCurrentAtomList( const unsigned& taskCode );
public:
  MultiColvarBase(const ActionOptions&);
  ~MultiColvarBase(){}
  static void registerKeywords( Keywords& keys );
/// Used in setupCurrentAtomList to get atom numbers 
/// Base quantities are different in MultiColvar and MultiColvarFunction
/// Turn on the derivatives 
  virtual void turnOnDerivatives();
/// Prepare for the calculation
/// Perform one of the tasks
  virtual void performTask();
/// This gets the position of an atom for the link cell setup
  virtual Vector getPositionOfAtomForLinkCells( const unsigned& iatom )=0;
/// And a virtual function which actually computes the colvar
  virtual double doCalculation();  
/// Update the atoms that have derivatives
  virtual void updateActiveAtoms()=0;
/// This is replaced once we have a function to calculate the cv
  virtual double compute()=0;
/// These replace the functions in ActionWithVessel to make the code faster
  virtual void mergeDerivatives( const unsigned& ider, const double& df );
  virtual void clearDerivativesAfterTask( const unsigned& ider );
/// Apply the forces from this action
  virtual void apply();
/// Get the number of derivatives for this action
  virtual unsigned getNumberOfDerivatives();  // N.B. This is replacing the virtual function in ActionWithValue
/// Get number size of atoms with derivatives array
  virtual unsigned getSizeOfAtomsWithDerivatives();
/// Checks if an task is being performed at the present time
  virtual bool isCurrentlyActive( const unsigned& code )=0;
/// Get the number of quantities that are calculated each time
  virtual unsigned getNumberOfQuantities();
/// Get the index where the central atom is stored
  virtual unsigned getCentralAtomElementIndex();
/// Retrieve the position of the central atom
  virtual Vector retrieveCentralAtomPos();
/// You can use this to screen contributions that are very small so we can avoid expensive (and pointless) calculations
  virtual void calculateWeight();
/// A virtual routine to get the position of the central atom - used for things like cv gradient
  virtual Vector calculateCentralAtomPosition()=0; 
/// Get the list of indices that have derivatives
 virtual void getIndexList( const unsigned& ntotal, const unsigned& jstore, const unsigned& maxder, std::vector<unsigned>& indices );
/// Is this a density?
  virtual bool isDensity(){ return false; }
/// Store central atoms so that this can be used in a function
  virtual vesselbase::StoreDataVessel* buildDataStashes( const bool& allow_wcutoff, const double& wtol );
/// Calculate and store getElementValue(uder)/getElementValue(vder) and its derivatives in getElementValue(iout)
  void quotientRule( const unsigned& uder, const unsigned& vder, const unsigned& iout );
/// Activate the atoms that have derivatives from a storeDataVessel
  void activateIndexes( const unsigned& istart, const unsigned& number, const std::vector<unsigned>& indexes ); 
/// Get the position of the iatom th central atom (used in multicolvarfunction)
  Vector getCentralAtomPosition( const unsigned& iatom ) const ;
/// Add central atom derivatives to a multicolvar function
  void addCentralAtomDerivativeToFunction( const unsigned& iatom, const unsigned& jout, const unsigned& base_cv_no, const Vector& der, MultiColvarFunction* func ); 
/// Get the value for this task
  virtual void getValueForTask( const unsigned& iatom, std::vector<double>& vals ); 
//// Used in ActionVolume and Gradient
  virtual void copyElementsToBridgedColvar( BridgedMultiColvarFunction* );
/// Used to accumulate values
  virtual void addWeightedValueDerivatives( const unsigned& iatom, const unsigned& base_cv_no, const double& weight, MultiColvarFunction* func );
/// Used for calculating weighted averages
  virtual void finishWeightedAverageCalculation( MultiColvarFunction* func );
/// Add derivatives to the orientations
  virtual void addOrientationDerivativesToBase( const unsigned& iatom, const unsigned& jstore, const unsigned& base_cv_no, 
                                                const std::vector<double>& weight, MultiColvarFunction* func );
/// Is the iatom'th stored value currently active
  bool storedValueIsActive( const unsigned& iatom );
/// This is true if multicolvar is calculating a vector or if the multicolvar is the density
  virtual bool hasDifferentiableOrientation() const { return false; }
/// This makes sure we are not calculating the director when we do LocalAverage
  virtual void doNotCalculateDirector(){}
};

inline
unsigned MultiColvarBase::getNumberOfDerivatives(){
  return 3*getNumberOfAtoms()+9;
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
unsigned MultiColvarBase::getCentralAtomElementIndex(){
  return 2;
}

inline
unsigned MultiColvarBase::getNAtoms() const {
  return natomsper;   // colvar_atoms[current].getNumberActive();
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
unsigned MultiColvarBase::getSizeOfAtomsWithDerivatives(){
  return getNumberOfAtoms();
}

}
}

#endif
