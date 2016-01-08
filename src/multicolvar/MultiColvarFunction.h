/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2016 The plumed team
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
#ifndef __PLUMED_multicolvar_MultiColvarFunction_h
#define __PLUMED_multicolvar_MultiColvarFunction_h

#include "tools/Matrix.h"
#include "MultiColvarBase.h"
#include "StoreCentralAtomsVessel.h"

namespace PLMD {
namespace multicolvar {

class MultiColvarFunction : public MultiColvarBase {
private:
/// The multicolvar from which we construct these quantities
  std::vector<multicolvar::MultiColvarBase*> mybasemulticolvars;
/// This is used to keep track of what is calculated where
  std::vector<unsigned> colvar_label;
/// Used for numerical derivatives
  Matrix<double> numder_store;
/// A tempory vector that is used for retrieving vectors
  std::vector<double> tvals;
protected:
/// Get the total number of tasks that this calculation is based on
  unsigned getFullNumberOfBaseTasks() const ;
/// Get the index of the ith colvar we are using
  unsigned getColvarIndex( const unsigned& ) const;
/// Get the position of one of the central atoms
  Vector getPositionOfCentralAtom(const unsigned&) const;
/// Add derivatives of value wrt to an atomic position 
  void addCentralAtomsDerivatives( const unsigned& , const unsigned& , const Vector& );
/// Retrieve the value calculated by the iatom th base task
  void getValueForBaseTask( const unsigned& iatom, std::vector<double>& vals );
/// Retrieve the vector calculated by the iatom th base task
  void getVectorForBaseTask( const unsigned& iatom, std::vector<double>& vecs );
/// Add the value and derivatives of this quantity
  void accumulateWeightedAverageAndDerivatives( const unsigned& iatom, const double& weight );
/// Add derivative wrt to the position of the central atom
  void addDerivativeOfCentralAtomPos( const unsigned& iatom, const Tensor& der );
/// Convert an index in the global array to an index in the individual base colvars
  unsigned convertToLocalIndex( const unsigned& index, const unsigned& mcv_code ) const ;
/// Add derivatives to the orientation
  void addOrientationDerivatives( const unsigned& iatom, const std::vector<double>& der );
/// Build colvars for atoms as if they were symmetry functions
  void buildSymmetryFunctionLists();
/// Build a colvar for each pair of atoms
  void buildAtomListWithPairs( const bool& allow_intra_group );
/// Get the number of base multicolvars 
  unsigned getNumberOfBaseMultiColvars() const ;
/// Get an example of the underlying multicolvar
  MultiColvarBase* getBaseMultiColvar( const unsigned& icolv );
/// Return the base multicolvar index that this colvar is a part of
  unsigned getBaseColvarNumber( const unsigned& iatom ) const ;
public:
  MultiColvarFunction(const ActionOptions&);
  static void registerKeywords( Keywords& keys );
/// Used to make sure we are calculating everything during neighbor list update step
  void unlockContributors();
  void lockContributors();
/// Active element in atoms_with_derivatives
  void atomHasDerivative( const unsigned& iatom );
/// Used to get atom numbers
  unsigned getBaseQuantityIndex( const unsigned& code );
/// Is a task currently being peformed
  bool isCurrentlyActive( const unsigned& code );
/// Finish task list update
  void finishTaskListUpdate();
/// Resize the dynamic arrays 
  void resizeDynamicArrays();
/// Update the atoms that are active
  void updateActiveAtoms();
/// Regular calculate
  void calculate();
/// Calculate the numerical derivatives for this action
  void calculateNumericalDerivatives( ActionWithValue* a=NULL );
/// Calculate the position of the central atom
  Vector calculateCentralAtomPosition();
/// Get the position of the central atom
  virtual Vector getCentralAtom()=0;
/// Add derivatives from storage vessels in MultiColvarBase
  void addStoredDerivative( const unsigned&, const unsigned&, const unsigned&, const double& );
/// This is used in MultiColvarBase only - it is used to setup the link cells
  Vector getPositionOfAtomForLinkCells( const unsigned& iatom );
};

inline
unsigned MultiColvarFunction::getBaseQuantityIndex( const unsigned& code ){
  return code;    
}

inline
bool MultiColvarFunction::isCurrentlyActive( const unsigned& code ){
  return true; // Tasks are all active because of store data vessel
}

inline
unsigned MultiColvarFunction::getFullNumberOfBaseTasks() const {
  return colvar_label.size();
}


inline
unsigned MultiColvarFunction::getNumberOfBaseMultiColvars() const {
  return mybasemulticolvars.size();
}

inline
MultiColvarBase* MultiColvarFunction::getBaseMultiColvar( const unsigned& icolv ){
  plumed_dbg_assert( icolv<mybasemulticolvars.size() );
  return mybasemulticolvars[icolv];
} 

inline
unsigned MultiColvarFunction::convertToLocalIndex( const unsigned& index, const unsigned& mcv_code ) const {
  unsigned t1 = index;
  for(unsigned k=0;k<mcv_code;++k) t1 -= mybasemulticolvars[k]->getFullNumberOfTasks();
  return t1;
}

inline
unsigned MultiColvarFunction::getBaseColvarNumber( const unsigned& iatom ) const {
  return colvar_label[ current_atoms[iatom] ];
}

inline
Vector MultiColvarFunction::getPositionOfAtomForLinkCells( const unsigned& iatom ){
  plumed_dbg_assert( iatom<getFullNumberOfBaseTasks() ); unsigned mmc=colvar_label[ iatom ];
  return mybasemulticolvars[mmc]->getCentralAtomPosition( convertToLocalIndex(iatom,mmc) );
}

inline
Vector MultiColvarFunction::getPositionOfCentralAtom( const unsigned& iatom ) const {
  plumed_dbg_assert( iatom<natomsper ); unsigned mmc = colvar_label[ current_atoms[iatom] ];
  return mybasemulticolvars[mmc]->getCentralAtomPosition( convertToLocalIndex(current_atoms[iatom],mmc) );   
}

inline
void MultiColvarFunction::addCentralAtomsDerivatives( const unsigned& iatom, const unsigned& jout, const Vector& der ){
  if( doNotCalculateDerivatives() ) return ;

  plumed_dbg_assert( iatom<natomsper ); unsigned mmc = colvar_label[ current_atoms[iatom] ]; 
  mybasemulticolvars[mmc]->addCentralAtomDerivativeToFunction( convertToLocalIndex(current_atoms[iatom],mmc), jout, mmc, der, this );
}

inline
void MultiColvarFunction::atomHasDerivative( const unsigned& iatom ){
  plumed_dbg_assert( !doNotCalculateDerivatives() );
  atoms_with_derivatives.activate( iatom );
}

inline
void MultiColvarFunction::addDerivativeOfCentralAtomPos( const unsigned& iatom, const Tensor& der ){
  if( doNotCalculateDerivatives() ) return;

  plumed_dbg_assert( iatom<natomsper ); unsigned mmc = colvar_label[ current_atoms[iatom] ]; Vector tmpder;
  for(unsigned i=0;i<3;++i){
      for(unsigned j=0;j<3;++j) tmpder[j]=der(i,j);
      mybasemulticolvars[mmc]->addCentralAtomDerivativeToFunction( convertToLocalIndex(current_atoms[iatom],mmc), (2+i), mmc, tmpder, this );
  }
}

inline
void MultiColvarFunction::getValueForBaseTask( const unsigned& iatom, std::vector<double>& vals ){
  plumed_dbg_assert( iatom<natomsper ); unsigned mmc = colvar_label[ current_atoms[iatom] ];
  mybasemulticolvars[mmc]->getValueForTask( convertToLocalIndex(current_atoms[iatom],mmc), vals );
}

inline
void MultiColvarFunction::getVectorForBaseTask( const unsigned& iatom, std::vector<double>& vec ){
  plumed_dbg_assert( vec.size()==mybasemulticolvars[0]->getNumberOfQuantities()-5 && tvals.size()>1 );
  getValueForBaseTask( iatom, tvals ); for(unsigned i=0;i<vec.size();++i) vec[i]=tvals[i+1];
}

inline
void MultiColvarFunction::accumulateWeightedAverageAndDerivatives( const unsigned& iatom, const double& weight ){
  if( doNotCalculateDerivatives() ) return;

  plumed_dbg_assert( iatom<natomsper ); unsigned mmc = colvar_label[ current_atoms[iatom] ];
  mybasemulticolvars[mmc]->addWeightedValueDerivatives( convertToLocalIndex(current_atoms[iatom],mmc), mmc, weight, this );
}

inline
void MultiColvarFunction::addOrientationDerivatives( const unsigned& iatom , const std::vector<double>& der ){
  if( doNotCalculateDerivatives() ) return;

  plumed_dbg_assert( iatom<natomsper ); unsigned mmc = colvar_label[ current_atoms[iatom] ];
  unsigned jout=2; if( usespecies && iatom==0 ) jout=1;
  mybasemulticolvars[mmc]->addOrientationDerivativesToBase( convertToLocalIndex(current_atoms[iatom],mmc), jout, mmc, der, this );
}

}
}
#endif
