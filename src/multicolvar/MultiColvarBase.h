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
/// A dynamic list containing those atoms with derivatives
  DynamicList<unsigned> atoms_with_derivatives;
/// The lists of the atoms involved in each of the individual colvars
/// note these refer to the atoms in all_atoms
  std::vector< DynamicList<unsigned> > colvar_atoms;
/// Variables used for central atoms
  Tensor ibox;
  bool centralAtomDerivativesAreInFractional;
  DynamicList<unsigned> atomsWithCatomDer;
  std::vector<Tensor> central_derivs;
/// The forces we are going to apply to things
  std::vector<double> forcesToApply;
/// This resizes the local arrays after neighbor list updates and during initialization
  void resizeLocalArrays();
protected:
/// Add a colvar to the set of colvars we are calculating (in practise just a list of atoms)
  void addColvar( const std::vector<unsigned>& newatoms );
/// Finish setting up the multicolvar base
  void setupMultiColvarBase();
/// Get the separation between a pair of vectors
  Vector getSeparation( const Vector& vec1, const Vector& vec2 ) const ;
/// Do we use pbc to calculate this quantity
  bool usesPbc() const ;
/// Return the index of an atom
  unsigned getAtomIndex( const unsigned& ) const ;
/// Add some derivatives for an atom 
  void addAtomsDerivatives(const int&,const Vector&);
/// Add some derivatives to the virial
  void addBoxDerivatives(const Tensor&);
/// Retrieve derivative of central atom position wrt jcomp'th component of position of iatom'th atom
  double getCentralAtomDerivative( const unsigned& iatom, const unsigned jcomp, const Vector& df ) const ;
/// Set a weight for this colvar (used in MEAN and HISTOGRAM)
  void setWeight( const double& weight );
/// Set the derivative of the weight (used in MEAN and HISTOGRAM)
  void addAtomsDerivativeOfWeight( const unsigned& i, const Vector& wder );
  void addBoxDerivativesOfWeight( const Tensor& vir );
/// Get the number of atoms in this particular colvar
  unsigned getNAtoms() const;
/// Update the list of atoms after the neighbor list step
  void removeAtomRequest( const unsigned& aa, const double& weight );
/// Add derivative of central atom position wrt to position of iatom'th atom
  void addCentralAtomDerivatives( const unsigned& iatom, const Tensor& der );
public:
  MultiColvarBase(const ActionOptions&);
  ~MultiColvarBase(){};
  static void registerKeywords( Keywords& keys );
/// Prepare for the calculation
  void prepare();
  virtual void resizeDynamicArrays()=0;
/// Return the size of the colvar_atoms array
  unsigned getNumberOfColvars() const ;
/// Perform one of the tasks
  void performTask( const unsigned& j );
/// And a virtual function which actually computes the colvar
  virtual double doCalculation( const unsigned& j )=0;  
/// These replace the functions in ActionWithVessel to make the code faster
  void mergeDerivatives( const unsigned& ider, const double& df );
  void clearDerivativesAfterTask( const unsigned& ider );
/// Apply the forces from this action
  void apply();
/// Deactivate one of the tasks
  void deactivate_task();
/// Get the number of derivatives for this action
  unsigned getNumberOfDerivatives();  // N.B. This is replacing the virtual function in ActionWithValue
/// Retrieve the position of the central atom
  Vector retrieveCentralAtomPos( const bool& frac );
/// You can use this to screen contributions that are very small so we can avoid expensive (and pointless) calculations
  virtual void calculateWeight();
/// A virtual routine to get the position of the central atom - used for things like cv gradient
  virtual Vector calculateCentralAtomPosition()=0; 
/// Is this a density?
  virtual bool isDensity(){ return false; }
/// Return a pointer to the vessel that stores the positions of 
/// all the central atoms
  StoreCentralAtomsVessel* getCentralAtoms();
};

inline
unsigned MultiColvarBase::getNumberOfDerivatives(){
  return 3*getNumberOfAtoms()+9;
}

inline
unsigned MultiColvarBase::getAtomIndex( const unsigned& iatom ) const {
  plumed_dbg_assert( iatom<colvar_atoms[current].getNumberActive() );
  return all_atoms.linkIndex( colvar_atoms[current][iatom] );
}

inline
void MultiColvarBase::removeAtomRequest( const unsigned& i, const double& weight ){
  if( !contributorsAreUnlocked ) return;
  plumed_dbg_assert( weight<getTolerance() );
  if( weight<getNLTolerance() ) colvar_atoms[current].deactivate( i );
}

inline
void MultiColvarBase::deactivate_task(){
  if( !contributorsAreUnlocked ) return;   // Deactivating tasks only possible during neighbor list update
  colvar_atoms[current].deactivateAll();   // Deactivate all atom requests for this colvar
  ActionWithVessel::deactivate_task();     // Deactivate the colvar from the list
}

inline
bool MultiColvarBase::usesPbc() const {
  return usepbc;
}

inline
unsigned MultiColvarBase::getNumberOfColvars() const {
  return colvar_atoms.size();
}

inline
unsigned MultiColvarBase::getNAtoms() const {
  return colvar_atoms[current].getNumberActive();
}

inline
void MultiColvarBase::addAtomsDerivatives(const int& iatom, const Vector& der){
  atoms_with_derivatives.activate(iatom);
  addElementDerivative( 3*iatom+0, der[0] );
  addElementDerivative( 3*iatom+1, der[1] );
  addElementDerivative( 3*iatom+2, der[2] );
} 

inline
void MultiColvarBase::addBoxDerivatives(const Tensor& vir){
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
void MultiColvarBase::calculateWeight(){
  setElementValue( 1, 1.0 );
}

inline
void MultiColvarBase::setWeight( const double& weight ){
  setElementValue( 1, weight );
}

inline
void MultiColvarBase::addAtomsDerivativeOfWeight( const unsigned& iatom, const Vector& wder ){
  unsigned nstart  = 3*getNumberOfAtoms() + 9 + 3*iatom;   
  atoms_with_derivatives.activate(iatom);
  addElementDerivative( nstart + 0, wder[0] );
  addElementDerivative( nstart + 1, wder[1] );
  addElementDerivative( nstart + 2, wder[2] );
}

inline
void MultiColvarBase::addBoxDerivativesOfWeight( const Tensor& vir ){
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
