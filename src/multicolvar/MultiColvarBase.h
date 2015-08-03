/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2015 The plumed team
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
#include <vector>

namespace PLMD {
namespace multicolvar {

class AtomValuePack;
class CatomPack;
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
/// The forces we are going to apply to things
  std::vector<double> forcesToApply;
/// Stuff for link cells - this is used to make coordination number like variables faster
  LinkCells linkcells;
/// This remembers where the boundaries are for the tasks. It makes link cells work fast
  Matrix<std::pair<unsigned,unsigned> > bookeeping;
/// Bool vector telling us which atoms are required to calculate central atom position
  std::vector<bool> use_for_central_atom;
/// 1/number of atoms involved in central atoms
  double numberForCentralAtom;
/// A copy of the vessel containg the values of each colvar
//  StoreColvarVessel* myvalues;
/// This resizes the arrays that are used for link cell update
  void resizeBookeepingArray( const unsigned& num1, const unsigned& num2 );
protected:
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
  void setupMultiColvarBase();
/// Set which atoms are to be used to calculate the central atom position
  void setAtomsForCentralAtom( const std::vector<bool>& catom_ind );
/// Set the value of the cutoff for the link cells
  void setLinkCellCutoff( const double& lcut );
/// Setup link cells in order to make this calculation faster
  void setupLinkCells();
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
/// Prepare for the calculation
/// Perform one of the tasks
  virtual void performTask( const unsigned& , const unsigned& , MultiValue& ) const ;
/// This gets the position of an atom for the link cell setup
  virtual Vector getPositionOfAtomForLinkCells( const unsigned& iatom ) const=0;
/// And a virtual function which actually computes the colvar
  virtual double doCalculation( const unsigned& tindex, AtomValuePack& myatoms ) const ;  
/// Update the atoms that have derivatives
  virtual void updateActiveAtoms( AtomValuePack& myatoms ) const=0;
/// This is replaced once we have a function to calculate the cv
  virtual double compute( const unsigned& tindex, AtomValuePack& myatoms ) const=0;
/// Apply the forces from this action
  virtual void apply();
/// Get the number of derivatives for this action
  virtual unsigned getNumberOfDerivatives();  // N.B. This is replacing the virtual function in ActionWithValue
/// Checks if an task is being performed at the present time
  virtual bool isCurrentlyActive( const unsigned& bno, const unsigned& code )=0;
///
  virtual CatomPack getCentralAtomPack( const unsigned& basn, const unsigned& curr );
/// Get the index where the central atom is stored
  virtual Vector getCentralAtomPos( const unsigned& curr );
/// You can use this to screen contributions that are very small so we can avoid expensive (and pointless) calculations
  virtual void calculateWeight( AtomValuePack& myatoms ) const ;
/// Get the list of indices that have derivatives
// virtual void getIndexList( const unsigned& ntotal, const unsigned& jstore, const unsigned& maxder, std::vector<unsigned>& indices );
/// Is this a density?
  virtual bool isDensity() const { return false; }
/// Store central atoms so that this can be used in a function
//  virtual vesselbase::StoreDataVessel* buildDataStashes( const bool& allow_wcutoff, const double& wtol );
/// Calculate and store getElementValue(uder)/getElementValue(vder) and its derivatives in getElementValue(iout)
//  void quotientRule( const unsigned& uder, const unsigned& vder, const unsigned& iout );
/// Activate the atoms that have derivatives from a storeDataVessel
//  void activateIndexes( const unsigned& istart, const unsigned& number, const std::vector<unsigned>& indexes ); 
/// Add central atom derivatives to a multicolvar function
//  void addCentralAtomDerivativeToFunction( const unsigned& iatom, const unsigned& jout, const unsigned& base_cv_no, const Vector& der, MultiColvarFunction* func ); 
/// Get the value for this task
//  virtual void getValueForTask( const unsigned& iatom, std::vector<double>& vals ); 
/// Used to accumulate values
//  virtual void addWeightedValueDerivatives( const unsigned& iatom, const unsigned& base_cv_no, const double& weight, MultiColvarFunction* func );
/// Used for calculating weighted averages
//  virtual void finishWeightedAverageCalculation( MultiColvarFunction* func );
/// Add derivatives to the orientations
//  virtual void addOrientationDerivativesToBase( const unsigned& iatom, const unsigned& jstore, const unsigned& base_cv_no, 
//                                                const std::vector<double>& weight, MultiColvarFunction* func );
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

}
}

#endif
