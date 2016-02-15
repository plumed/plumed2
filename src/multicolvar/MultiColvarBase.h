/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2016 The plumed team
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
#ifndef __PLUMED_multicolvar_MultiColvarBase_h
#define __PLUMED_multicolvar_MultiColvarBase_h

#include "core/ActionAtomistic.h"
#include "core/ActionWithValue.h"
#include "tools/DynamicList.h"
#include "tools/LinkCells.h"
#include "vesselbase/StoreDataVessel.h"
#include "vesselbase/ActionWithVessel.h"
#include <vector>

namespace PLMD {
namespace multicolvar {

class AtomValuePack;
class CatomPack;
class BridgedMultiColvarFunction;
class ActionVolume;

class MultiColvarBase :
  public ActionAtomistic,
  public ActionWithValue,
  public vesselbase::ActionWithVessel
  {
 friend class BridgedMultiColvarFunction;
 friend class VolumeGradientBase;
 friend class MultiColvarFilter;
 friend class AtomValuePack;
private:
/// Use periodic boundary conditions
  bool usepbc;
/// The forces we are going to apply to things
  std::vector<double> forcesToApply;
/// We use this to say that all the atoms in the third block should are in the tasks
  bool allthirdblockintasks;
/// In certain cases we can make three atom link cells faster
  bool uselinkforthree;
/// Stuff for link cells - this is used to make coordination number like variables faster
  LinkCells linkcells;
/// Link cells for third block of atoms
  LinkCells threecells;
/// Bool vector telling us which atoms are required to calculate central atom position
  std::vector<bool> use_for_central_atom;
/// 1/number of atoms involved in central atoms
  double numberForCentralAtom;
/// Ensures that setup is only performed once per loop
  bool setup_completed;
/// Ensures that retrieving of atoms is only done once per calculation loop
  bool atomsWereRetrieved;
protected:
/// This is used to keep track of what is calculated where
  std::vector<unsigned> colvar_label;
/// The multicolvars from which we construct these quantities
  std::vector<MultiColvarBase*> mybasemulticolvars;
/// The vessels in these multicolvars in which the data is stored
  std::vector<vesselbase::StoreDataVessel*> mybasedata;
/// This remembers where the boundaries are for the tasks. It makes link cells work fast
  Matrix<std::pair<unsigned,unsigned> > bookeeping;
/// Function that recursively checks if filters have been used in the input to a multicolvar
/// we need this to ensure that setupLinkCells is run in calculate with some actions
  bool filtersUsedAsInput();
/// Read in a set of multicolvar labels as the input to the action
  bool interpretInputMultiColvars( const std::vector<std::string>& key, const double& wtolerance );
/// Convert an index in the global array to an index in the individual base colvars
  unsigned convertToLocalIndex( const unsigned& index, const unsigned& mcv_code ) const ;
/// This resizes the arrays that are used for link cell update
  void resizeBookeepingArray( const unsigned& num1, const unsigned& num2 );
/// Using the species keyword to read in atoms
  bool usespecies;
/// Number of atoms in each block
  unsigned nblock;
/// This is used when turning cvcodes into atom numbers
  std::vector<unsigned> decoder;
/// Blocks of atom numbers
  std::vector< std::vector<unsigned> > ablocks;
/// Add a task to the list of tasks
  void addTaskToList( const unsigned& taskCode );
/// Finish setting up the multicolvar base
  void setupMultiColvarBase( const std::vector<AtomNumber>& atoms );
/// Add some derivatives to a particular component of a particular atom
  void addAtomDerivatives( const int& , const unsigned& , const Vector& , multicolvar::AtomValuePack& ) const ;
/// Set which atoms are to be used to calculate the central atom position
  void setAtomsForCentralAtom( const std::vector<bool>& catom_ind );
/// Set the value of the cutoff for the link cells
  void setLinkCellCutoff( const double& lcut, double tcut=-1.0 );
/// Setup the link cells and neighbour list stuff
  void setupActiveTaskSet( std::vector<unsigned>& active_tasks, const std::string& input_label );
/// Setup link cells in order to make this calculation faster
  void setupLinkCells();
/// This does setup of link cell stuff that is specific to the non-use of the usespecies keyword
  void setupNonUseSpeciesLinkCells( const unsigned& );
/// Get the separation between a pair of vectors
  Vector getSeparation( const Vector& vec1, const Vector& vec2 ) const ;
/// This sets up the list of atoms that are involved in this colvar
  bool setupCurrentAtomList( const unsigned& taskCode, AtomValuePack& myatoms ) const ;
/// Decode indices if there are 2 or 3 atoms involved
  void decodeIndexToAtoms( const unsigned& taskCode, std::vector<unsigned>& atoms ) const ;
public:
  explicit MultiColvarBase(const ActionOptions&);
  ~MultiColvarBase(){}
  static void registerKeywords( Keywords& keys );
/// Turn on the derivatives 
  virtual void turnOnDerivatives();
/// Do we use pbc to calculate this quantity
  bool usesPbc() const ;
/// Apply PBCs over a set of distance vectors
  void applyPbc(std::vector<Vector>& dlist, unsigned max_index=0) const;
/// Do some setup before the calculation
  void prepare();
/// This is overwritten here in order to make sure that we do not retrieve atoms multiple times
  void retrieveAtoms();
/// Do the calculation
  virtual void calculate();
/// Calculate numerical derivatives
  virtual void calculateNumericalDerivatives( ActionWithValue* a=NULL );
/// Perform one of the tasks
  virtual void performTask( const unsigned& , const unsigned& , MultiValue& ) const ;
/// Update the active atoms
  virtual void updateActiveAtoms( AtomValuePack& myatoms ) const ;
/// This gets the position of an atom for the link cell setup
  virtual Vector getPositionOfAtomForLinkCells( const unsigned& iatom ) const ;
/// Returns the position where we should assume the center is for link cell calculations
  virtual Vector getLinkCellPosition( const std::vector<unsigned>& atoms ) const ;
/// And a virtual function which actually computes the colvar
  virtual double doCalculation( const unsigned& tindex, AtomValuePack& myatoms ) const ;  
/// Get the absolute index of the central atom
  virtual AtomNumber getAbsoluteIndexOfCentralAtom( const unsigned& i ) const ;
/// This is replaced once we have a function to calculate the cv
  virtual double compute( const unsigned& tindex, AtomValuePack& myatoms ) const=0;
/// Apply the forces from this action
  virtual void apply();
/// Get the number of derivatives for this action
  virtual unsigned getNumberOfDerivatives();  // N.B. This is replacing the virtual function in ActionWithValue
/// Checks if an task is being performed at the present time
  virtual bool isCurrentlyActive( const unsigned& bno, const unsigned& code );
///
  virtual CatomPack getCentralAtomPack( const unsigned& basn, const unsigned& curr );
/// Get the index where the central atom is stored
  virtual Vector getCentralAtomPos( const unsigned& curr );
/// You can use this to screen contributions that are very small so we can avoid expensive (and pointless) calculations
  virtual void calculateWeight( const unsigned& taskCode, AtomValuePack& myatoms ) const ;
/// Is this a density?
  virtual bool isDensity() const { return false; }
/// Is the iatom'th stored value currently active
  bool storedValueIsActive( const unsigned& iatom );
/// This is true if multicolvar is calculating a vector or if the multicolvar is the density
  virtual bool hasDifferentiableOrientation() const { return false; }
/// This makes sure we are not calculating the director when we do LocalAverage
  virtual void doNotCalculateDirector(){}
/// Ensure that derivatives are only calculated when needed
  bool doNotCalculateDerivatives() const ;
/// Get the icolv th base multicolvar 
  MultiColvarBase* getBaseMultiColvar( const unsigned& icolv ) const ;
/// Get the number of base multicolvars 
  unsigned getNumberOfBaseMultiColvars() const ;
};

inline
unsigned MultiColvarBase::convertToLocalIndex( const unsigned& index, const unsigned& mcv_code ) const {
  unsigned t1 = index;
  for(unsigned k=0;k<mcv_code;++k) t1 -= mybasemulticolvars[k]->getFullNumberOfTasks();
  return t1;
}

inline
bool MultiColvarBase::isCurrentlyActive( const unsigned& bno, const unsigned& code ){
  if( setup_completed && code<colvar_label.size() ){
     unsigned mmc=colvar_label[code]; 
     return mybasedata[mmc]->storedValueIsActive( convertToLocalIndex(code,mmc) ); 
  }
  return true;
}

inline
AtomNumber MultiColvarBase::getAbsoluteIndexOfCentralAtom( const unsigned& iatom ) const {
  if( iatom<colvar_label.size() ){
      unsigned mmc=colvar_label[ iatom ];
      return mybasemulticolvars[mmc]->getAbsoluteIndexOfCentralAtom( convertToLocalIndex(iatom,mmc) ); 
  }
  return ActionAtomistic::getAbsoluteIndex( getTaskCode(iatom) );
} 

inline
Vector MultiColvarBase::getPositionOfAtomForLinkCells( const unsigned& iatom ) const {
  if( iatom<colvar_label.size()  ){ 
      unsigned mmc=colvar_label[ iatom ];
      return mybasemulticolvars[mmc]->getCentralAtomPos( convertToLocalIndex(iatom,mmc) );
  }
  return ActionAtomistic::getPosition( iatom );
}

inline
Vector MultiColvarBase::getLinkCellPosition( const std::vector<unsigned>& atoms ) const {
  return getPositionOfAtomForLinkCells( atoms[0] );
} 

inline
unsigned MultiColvarBase::getNumberOfDerivatives(){
  return 3*getNumberOfAtoms()+9;
}

inline
bool MultiColvarBase::usesPbc() const {
  return usepbc;
}

inline
bool MultiColvarBase::doNotCalculateDerivatives() const {
  if( !dertime ) return true;
  return ActionWithValue::doNotCalculateDerivatives();
}

inline
unsigned MultiColvarBase::getNumberOfBaseMultiColvars() const {
  return mybasemulticolvars.size(); 
} 

inline 
MultiColvarBase* MultiColvarBase::getBaseMultiColvar( const unsigned& icolv ) const {
  plumed_dbg_assert( icolv<mybasemulticolvars.size() );
  return mybasemulticolvars[icolv]; 
} 


}
}

#endif
