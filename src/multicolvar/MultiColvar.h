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
#ifndef __PLUMED_MultiColvar_h
#define __PLUMED_MultiColvar_h

#include "basic/ActionAtomistic.h"
#include "basic/ActionWithValue.h"
#include "ActionWithDistribution.h"
#include <vector>

#define PLUMED_MULTICOLVAR_INIT(ao) Action(ao),MultiColvar(ao)

namespace PLMD {

class MultiColvar;

/**
\ingroup INHERIT
This is the abstract base class to use for creating distributions of colvars and functions
thereof, whtin it there is \ref AddingAMultiColvar "information" as to how to go implementing these types of actions.
*/

class MultiColvar :
  public ActionAtomistic,
  public ActionWithValue,
  public ActionWithDistribution
  {
private:
  bool usepbc;
  bool readatoms;
  bool verbose_output;
  bool needsCentralAtomPosition;
/// The list of all the atoms involved in the colvar
  DynamicList<AtomNumber> all_atoms;
/// The lists of the atoms involved in each of the individual colvars
/// note these refer to the atoms in all_atoms
  std::vector< DynamicList<unsigned> > colvar_atoms;
/// These are used to store the values of CVs etc so they can be retrieved by distribution
/// functions
  std::vector<Vector> pos;
  std::vector<Tensor> central_derivs;
  Value thisval;
  std::vector<Value> catom_pos;
/// Used to make sure we update the correct atoms during neighbor list update
  unsigned current;
/// Read in ATOMS keyword
  void readAtomsKeyword( int& natoms );
/// Read in the various GROUP keywords
  void readGroupsKeyword( int& natoms );
/// Read in the various SPECIES keywords
  void readSpeciesKeyword( int& natoms );
/// Update the atoms request
  void requestAtoms();
protected:
/// Read in all the keywords that can be used to define atoms
  void readAtoms( int& natoms );
/// Read in the atoms that form the backbone of a polymeric chain
  void readBackboneAtoms( const std::vector<std::string>& backnames, std::vector<unsigned>& chain_lengths );
/// Add a colvar to the set of colvars we are calculating (in practise just a list of atoms)
  void addColvar( const std::vector<unsigned>& newatoms );
/// Get the separation between a pair of vectors
  Vector getSeparation( const Vector& vec1, const Vector& vec2 ) const ;
/// Update the list of atoms after the neighbor list step
  void removeAtomRequest( const unsigned& aa );
/// Do we use pbc to calculate this quantity
  bool usesPbc() const ;
/// Add some derivatives for an atom 
  void addAtomsDerivatives(const int&,const Vector&);
/// Add some derivatives to the virial
  void addBoxDerivatives(const Tensor&);
public:
  MultiColvar(const ActionOptions&);
  ~MultiColvar(){};
  static void registerKeywords( Keywords& keys );
/// Calculate the multicolvar
  void calculate();
/// Prepare for the calculation
  void prepare();
/// Apply the forces on the values
  void apply();
/// Get the number of derivatives for this action
  unsigned getNumberOfDerivatives();  // N.B. This is replacing the virtual function in ActionWithValue
/// Return the number of Colvars this is calculating
  unsigned getNumberOfFunctionsInAction();  
/// Return the number of derivatives for a given colvar
  unsigned getNumberOfDerivatives( const unsigned& j );
/// Retrieve the value that was calculated last
  const Value & retreiveLastCalculatedValue();
/// Retrieve the position of the central atom
  void retrieveCentralAtomPos( std::vector<Value>& cpos ) const ;
/// Make sure we calculate the position of the central atom
  void useCentralAtom();
/// Get the weight of the colvar
  virtual void retrieveColvarWeight( const unsigned& i, Value& ww );
/// Merge the derivatives 
  void mergeDerivatives( const unsigned j, const Value& value_in, const double& df, const unsigned& vstart, Vessel* valout );
  void mergeDerivatives( const unsigned j, const Value& value_in, const double& df, Value* valout );
/// Turn of atom requests when this colvar is deactivated cos its small
  void deactivateValue( const unsigned j );
/// Turn on atom requests when the colvar is activated
  void activateValue( const unsigned j );
/// Calcualte the colvar
  bool calculateThisFunction( const unsigned& j );
/// You can use this to screen contributions that are very small so we can avoid expensive (and pointless) calculations
  virtual bool contributionIsSmall( std::vector<Vector>& pos ){ plumed_assert( !isPossibleToSkip() ); return false; }
/// And a virtual function which actually computes the colvar
  virtual double compute( const unsigned& j, const std::vector<Vector>& pos )=0;  
/// A virtual routine to get the position of the central atom - used for things like cv gradient
  virtual void getCentralAtom( const std::vector<Vector>& pos, Vector& cpos, std::vector<Tensor>& deriv ); 
/// Is this a density?
  virtual bool isDensity(){ return false; }
};

inline
unsigned MultiColvar::getNumberOfDerivatives(){
  return 3*getNumberOfAtoms()+9;
} 

inline
unsigned MultiColvar::getNumberOfFunctionsInAction(){
  return colvar_atoms.size();
}

inline
void MultiColvar::deactivateValue( const unsigned j ){
  colvar_atoms[j].deactivateAll();
}

inline
void MultiColvar::activateValue( const unsigned j ){
  colvar_atoms[j].activateAll(); 
  colvar_atoms[j].updateActiveMembers();
}

inline
const Value & MultiColvar::retreiveLastCalculatedValue(){
  return thisval;
  // copy( thisval, myvalue );  
}

inline
unsigned MultiColvar::getNumberOfDerivatives( const unsigned& j ){
  return 3*colvar_atoms[j].getNumberActive() + 9;
}

inline
void MultiColvar::removeAtomRequest( const unsigned& i ){
  plumed_massert(isTimeForNeighborListUpdate(),"found removeAtomRequest but not during neighbor list step");
  colvar_atoms[current].deactivate(i); 
}

inline
bool MultiColvar::usesPbc() const {
  return usepbc;
}

inline
void MultiColvar::addAtomsDerivatives(const int& iatom, const Vector& der){
  plumed_assert( iatom<colvar_atoms[current].getNumberActive() );
  thisval.addDerivative( 3*iatom+0, der[0] );
  thisval.addDerivative( 3*iatom+1, der[1] );
  thisval.addDerivative( 3*iatom+2, der[2] );
}

inline
void MultiColvar::addBoxDerivatives(const Tensor& vir){
  int natoms=colvar_atoms[current].getNumberActive();
  thisval.addDerivative( 3*natoms+0, vir(0,0) );
  thisval.addDerivative( 3*natoms+1, vir(0,1) );
  thisval.addDerivative( 3*natoms+2, vir(0,2) );
  thisval.addDerivative( 3*natoms+3, vir(1,0) );
  thisval.addDerivative( 3*natoms+4, vir(1,1) );
  thisval.addDerivative( 3*natoms+5, vir(1,2) );
  thisval.addDerivative( 3*natoms+6, vir(2,0) );
  thisval.addDerivative( 3*natoms+7, vir(2,1) );
  thisval.addDerivative( 3*natoms+8, vir(2,2) );
}

}

#endif
