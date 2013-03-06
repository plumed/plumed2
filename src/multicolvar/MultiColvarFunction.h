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
#ifndef __PLUMED_multicolvar_MultiColvarFunction_h
#define __PLUMED_multicolvar_MultiColvarFunction_h

#include "core/ActionWithValue.h"
#include "vesselbase/ActionWithVessel.h"
#include "MultiColvar.h"
#include "StoreCentralAtomsVessel.h"

namespace PLMD {
namespace multicolvar {

class MultiColvarFunction :
  public ActionWithValue,
  public vesselbase::ActionWithVessel
  {
private:
  bool usepbc;
  bool readatoms;
  bool verbose_output;
/// Everything for controlling the updating of neighbor lists
  int updateFreq;
  unsigned lastUpdate;
  bool reduceAtNextStep;
/// A DynamicList containing the numbers of 1 - ncolvars
  DynamicList<unsigned> colvar_list;
/// The lists of the atoms involved in each of the individual colvars
/// note these refer to the atoms in all_atoms
  std::vector< DynamicList<unsigned> > colvar_atoms;
/// The multicolvar from which we construct these quantities
  multicolvar::MultiColvar* mycolv;
/// The central atom positions
  multicolvar::StoreCentralAtomsVessel* catoms;
protected:
/// Get the number of functions in the multicolvar we are operating on
  unsigned getNumberOfBaseFunctions() const;
/// Return a pointer to the multicolvar we are using as a base function
  MultiColvar* getPntrToMultiColvar();
/// Add a colvar to the set of colvars we are calculating (in practise just a list of atoms)
  void addColvar( const std::vector<unsigned>& newatoms );
/// Finish off the setup of the VectorFunction
  void completeSetup();
/// Find out if it is time to do neighbor list update
  bool isTimeForNeighborListUpdate() const ;
/// Update the list of atoms after the neighbor list step
  void removeAtomRequest( const unsigned& aa );
/// Get the number of atoms in this particular colvar
  unsigned getNAtoms() const;
/// Get the index of an atom
  unsigned getAtomIndex(const unsigned&) const;
/// Get the position of one of the central atoms
  Vector getPositionOfCentralAtom(const unsigned&) const;
/// Get the separation between a pair of positions
  Vector getSeparation( const Vector& , const Vector& ) const ;
/// Add derivatives of value wrt to an atomic position 
  void addAtomsDerivatives( const unsigned& , const Vector& );
/// Add derivatives of value wrt to the box
  void addBoxDerivatives( const Tensor& );
/// Add derivatives of weight wrt to an atomic position
  void addAtomsDerivativesOfWeight( const unsigned& , const Vector& );
/// Add derivatives of weight wrt to the box
  void addBoxDerivativesOfWeight( const Tensor& );
public:
  MultiColvarFunction(const ActionOptions&);
  static void registerKeywords( Keywords& keys );
/// Calculate the multicolvar
  void calculate();
// Calculate the numerical derivatives for this action
  void calculateNumericalDerivatives( ActionWithValue* a=NULL );
/// Prepare for the calculation
  virtual void prepare();
/// Apply the forces on the values
  void apply();
/// Get the number of derivatives for this action
  unsigned getNumberOfDerivatives();  // N.B. This is replacing the virtual function in ActionWithValue
/// Return the number of derivatives for a given colvar
  unsigned getNumberOfDerivatives( const unsigned& j );
/// Get the number of vectors we are looping over
  unsigned getNumberOfFunctionsInAction();
/// Turn of atom requests when this colvar is deactivated cos its small
  void deactivate_task();
/// Turn on atom requests when the colvar is activated
  void activateValue( const unsigned j );
/// Perform the task
  bool performTask( const unsigned& j );
/// Calculate the weight
  virtual double calculateWeight()=0;
/// Actually compute one of the colvars
  virtual double compute()=0;
};

inline
unsigned MultiColvarFunction::getNumberOfBaseFunctions() const {
  return mycolv->getNumberOfFunctionsInAction();
} 

inline
unsigned MultiColvarFunction::getNumberOfFunctionsInAction(){
  return colvar_atoms.size();
}

inline
MultiColvar* MultiColvarFunction::getPntrToMultiColvar(){
  return mycolv;
}

inline
unsigned MultiColvarFunction::getNumberOfDerivatives(){
  return mycolv->getNumberOfDerivatives();
}

inline
unsigned MultiColvarFunction::getNumberOfDerivatives( const unsigned& j ){
  return mycolv->getNumberOfDerivatives(j);
}

inline
void MultiColvarFunction::deactivate_task(){
  if( !reduceAtNextStep ) return;          // Deactivating tasks only possible during neighbor list update
  colvar_list.deactivate(current);         // Deactivate the colvar from the list
  colvar_atoms[current].deactivateAll();   // Deactivate all atom requests for this colvar
}

inline
void MultiColvarFunction::activateValue( const unsigned j ){
  colvar_atoms[j].activateAll(); 
  colvar_atoms[j].updateActiveMembers();
}

inline
bool MultiColvarFunction::isTimeForNeighborListUpdate() const {
  return reduceAtNextStep;
}

inline
void MultiColvarFunction::removeAtomRequest( const unsigned& i ){
  plumed_massert(reduceAtNextStep,"found removeAtomRequest but not during neighbor list step");
  colvar_atoms[current].deactivate(i);
}

inline
unsigned MultiColvarFunction::getNAtoms() const {
  return colvar_atoms[current].getNumberActive();
}

inline
unsigned MultiColvarFunction::getAtomIndex( const unsigned& iatom ) const {
  plumed_dbg_assert( iatom<colvar_atoms[current].getNumberActive() );
  return colvar_atoms[current][iatom];
}

inline
Vector MultiColvarFunction::getPositionOfCentralAtom( const unsigned& iatom ) const {
  plumed_dbg_assert( iatom<colvar_atoms[current].getNumberActive() );
  return catoms->getPosition( colvar_atoms[current][iatom] );
}

inline
void MultiColvarFunction::addAtomsDerivatives( const unsigned& iatom, const Vector& der ){
  plumed_dbg_assert( iatom<colvar_atoms[current].getNumberActive() );
  catoms->chainRuleForCentralAtom( colvar_atoms[current][iatom], 0, der, this ); 
}

inline
void MultiColvarFunction::addBoxDerivatives( const Tensor& vir ){
  unsigned nstart=mycolv->getNumberOfDerivatives() - 9;
  addElementDerivative( nstart + 0, vir(0,0) );
  addElementDerivative( nstart + 1, vir(0,1) );
  addElementDerivative( nstart + 2, vir(0,2) );
  addElementDerivative( nstart + 3, vir(1,0) );
  addElementDerivative( nstart + 4, vir(1,1) );
  addElementDerivative( nstart + 5, vir(1,2) );
  addElementDerivative( nstart + 6, vir(2,0) );
  addElementDerivative( nstart + 7, vir(2,1) );
  addElementDerivative( nstart + 8, vir(2,2) );
}

inline
void MultiColvarFunction::addAtomsDerivativesOfWeight( const unsigned& iatom, const Vector& der ){
  plumed_dbg_assert( iatom<colvar_atoms[current].getNumberActive() );
  catoms->chainRuleForCentralAtom( colvar_atoms[current][iatom], 1, der, this );
}

inline
void MultiColvarFunction::addBoxDerivativesOfWeight( const Tensor& vir ){
  unsigned nstart=2*mycolv->getNumberOfDerivatives() - 9;
  addElementDerivative( nstart + 0, vir(0,0) );
  addElementDerivative( nstart + 1, vir(0,1) );
  addElementDerivative( nstart + 2, vir(0,2) );
  addElementDerivative( nstart + 3, vir(1,0) );
  addElementDerivative( nstart + 4, vir(1,1) );
  addElementDerivative( nstart + 5, vir(1,2) );
  addElementDerivative( nstart + 6, vir(2,0) );
  addElementDerivative( nstart + 7, vir(2,1) );
  addElementDerivative( nstart + 8, vir(2,2) );
}

}
}
#endif
