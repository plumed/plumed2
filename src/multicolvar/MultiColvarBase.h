/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#include "CatomPack.h"
#include <vector>

namespace PLMD {
namespace multicolvar {

class AtomValuePack;
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
/// Number of atoms that are active on this step
  unsigned nactive_atoms;
/// Stuff for link cells - this is used to make coordination number like variables faster
  LinkCells linkcells;
/// Link cells for third block of atoms
  LinkCells threecells;
/// Number of atoms that are being used for central atom position
  unsigned ncentral;
/// Bool vector telling us which atoms are required to calculate central atom position
  std::vector<bool> use_for_central_atom;
/// Vector of tempory holders for central atom values
  std::vector<CatomPack> my_tmp_capacks;
/// 1/number of atoms involved in central atoms
  double numberForCentralAtom;
/// Ensures that setup is only performed once per loop
  bool setup_completed;
/// Ensures that retrieving of atoms is only done once per calculation loop
  bool atomsWereRetrieved;
/// Add derivatives of center of mass position
  void addComDerivatives( const int& ival, const unsigned& iatom, const Vector& der, multicolvar::AtomValuePack& myatoms ) const ;
protected:
/// This is used to keep track of what is calculated where
  std::vector< std::pair<unsigned,unsigned> > atom_lab;
/// The vessels in these multicolvars in which the data is stored
  std::vector<vesselbase::StoreDataVessel*> mybasedata;
/// The multicolvars from which we construct these quantities
  std::vector<MultiColvarBase*> mybasemulticolvars;
/// This remembers where the boundaries are for the tasks. It makes link cells work fast
  Matrix<std::pair<unsigned,unsigned> > bookeeping;
/// Function that recursively checks if filters have been used in the input to a multicolvar
/// we need this to ensure that setupLinkCells is run in calculate with some actions
  bool filtersUsedAsInput();
/// This resizes the arrays that are used for link cell update
  void resizeBookeepingArray( const unsigned& num1, const unsigned& num2 );
/// Are we doing sums of matrix rows
  bool matsums;
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
/// This routine take the vector of input derivatives and adds all the vectors to ivalth output derivatives
/// In other words end-start sets of derivatives come in and one set of derivatives come out
  void mergeInputDerivatives( const unsigned& ival, const unsigned& start, const unsigned& end, const unsigned& jatom,
                              const std::vector<double>& der, MultiValue& myder, AtomValuePack& myatoms ) const ;
/// This routine take the ith set of input derivatives and adds it to each of the (end-start) output derivatives
/// In other words one set of derivatives comes in and end-start sets of derivatives come out
  void splitInputDerivatives( const unsigned& ival, const unsigned& start, const unsigned& end,
                              const unsigned& jatom, const std::vector<double>& der,
                              MultiValue& myder, AtomValuePack& myatoms ) const ;
/// This is used to accumulate functions of the coordination sphere.  Ensures weights are taken into account
  void accumulateSymmetryFunction( const int& ival, const unsigned& iatom, const double& val, const Vector& der, const Tensor& vir, multicolvar::AtomValuePack& myatoms ) const ;
/// Set which atoms are to be used to calculate the central atom position
  void setAtomsForCentralAtom( const std::vector<bool>& catom_ind );
/// Set the value of the cutoff for the link cells
  void setLinkCellCutoff( const double& lcut, double tcut=-1.0 );
/// Setup the link cells and neighbour list stuff
  void setupActiveTaskSet( std::vector<unsigned>& active_tasks, const std::string& input_label );
/// Setup link cells in order to make this calculation faster
  void setupLinkCells();
/// Get the cutoff for the link cells
  double getLinkCellCutoff()  const ;
/// This does setup of link cell stuff that is specific to the non-use of the usespecies keyword
  void setupNonUseSpeciesLinkCells( const unsigned& );
/// This sets up the list of atoms that are involved in this colvar
  bool setupCurrentAtomList( const unsigned& taskCode, AtomValuePack& myatoms ) const ;
/// Decode indices if there are 2 or 3 atoms involved
  void decodeIndexToAtoms( const unsigned& taskCode, std::vector<unsigned>& atoms ) const ;
/// Read in some atoms
  bool parseMultiColvarAtomList(const std::string& key, const int& num, std::vector<AtomNumber>& t);
/// Read in ATOMS keyword
  void readAtomsLikeKeyword( const std::string & key, const int& natoms, std::vector<AtomNumber>& all_atoms );
/// Read in two groups of atoms and setup multicolvars to calculate
  void readTwoGroups( const std::string& key0, const std::string& key1, const std::string& key2, std::vector<AtomNumber>& all_atoms );
/// Read in three groups of atoms
  void readGroupKeywords( const std::string& key0, const std::string& key1, const std::string& key2, const std::string& key3,
                          const bool& no_third_dim_accum, const bool& symmetric, std::vector<AtomNumber>& all_atoms );
/// Read in three groups of atoms and construct CVs involving at least one
  void readThreeGroups( const std::string& key1, const std::string& key2, const std::string& key3,
                        const bool& allow2, const bool& no_third_dim_accum, std::vector<AtomNumber>& all_atoms );
/// Build sets by taking one multicolvar from each base
  void buildSets();
public:
  explicit MultiColvarBase(const ActionOptions&);
  ~MultiColvarBase() {}
  static void registerKeywords( Keywords& keys );
/// Turn on the derivatives
  void turnOnDerivatives() override;
/// Get the separation between a pair of vectors
  Vector getSeparation( const Vector& vec1, const Vector& vec2 ) const ;
/// Do we use pbc to calculate this quantity
  bool usesPbc() const ;
/// Apply PBCs over a set of distance vectors
  void applyPbc(std::vector<Vector>& dlist, unsigned max_index=0) const;
/// Is it safe to use multithreading
  bool threadSafe() const override { return !(mybasemulticolvars.size()>0); }
/// Do some setup before the calculation
  void prepare() override;
/// This is overwritten here in order to make sure that we do not retrieve atoms multiple times
  void retrieveAtoms() override;
/// Do the calculation
  void calculate() override;
/// Calculate numerical derivatives
  void calculateNumericalDerivatives( ActionWithValue* a=NULL ) override;
/// Perform one of the tasks
  void performTask( const unsigned&, const unsigned&, MultiValue& ) const override;
/// Update the active atoms
  virtual void updateActiveAtoms( AtomValuePack& myatoms ) const ;
/// This gets the position of an atom for the link cell setup
  virtual Vector getPositionOfAtomForLinkCells( const unsigned& iatom ) const ;
/// Get the absolute index of the central atom
  virtual AtomNumber getAbsoluteIndexOfCentralAtom( const unsigned& i ) const ;
/// This is replaced once we have a function to calculate the cv
  virtual double compute( const unsigned& tindex, AtomValuePack& myatoms ) const=0;
/// Apply the forces from this action
  void apply() override;
/// Get the number of derivatives for this action
  unsigned getNumberOfDerivatives() override;  // N.B. This is replacing the virtual function in ActionWithValue
/// Checks if an task is being performed at the present time
  virtual bool isCurrentlyActive( const unsigned& code );
/// Add some derivatives to a particular component of a particular atom
  void addAtomDerivatives( const int&, const unsigned&, const Vector&, multicolvar::AtomValuePack& ) const ;
///
  virtual void getCentralAtomPack( const unsigned& basn, const unsigned& curr, CatomPack& mypack);
/// Get the index where the central atom is stored
  virtual Vector getCentralAtomPos( const unsigned& curr );
/// You can use this to screen contributions that are very small so we can avoid expensive (and pointless) calculations
  virtual double calculateWeight( const unsigned& taskCode, const double& weight, AtomValuePack& myatoms ) const ;
/// Is this a density?
  virtual bool isDensity() const { return false; }
/// Is the iatom'th stored value currently active
  bool storedValueIsActive( const unsigned& iatom );
/// This is true if multicolvar is calculating a vector or if the multicolvar is the density
  virtual bool hasDifferentiableOrientation() const { return false; }
/// This makes sure we are not calculating the director when we do LocalAverage
  virtual void doNotCalculateDirector() {}
/// Ensure that derivatives are only calculated when needed
  bool doNotCalculateDerivatives() const override;
/// Get the icolv th base multicolvar
  MultiColvarBase* getBaseMultiColvar( const unsigned& icolv ) const ;
/// Get the number of base multicolvars
  unsigned getNumberOfBaseMultiColvars() const ;
/// Get an input data
  virtual void getInputData( const unsigned& ind, const bool& normed, const multicolvar::AtomValuePack& myatoms, std::vector<double>& orient ) const ;
/// Retrieve the input derivatives
  virtual MultiValue& getInputDerivatives( const unsigned& iatom, const bool& normed, const multicolvar::AtomValuePack& myatoms ) const ;
};

inline
bool MultiColvarBase::isCurrentlyActive( const unsigned& code ) {
  if( setup_completed && atom_lab[code].first>0 ) {
    unsigned mmc=atom_lab[code].first - 1;
    return mybasedata[mmc]->storedValueIsActive( atom_lab[code].second );
  }
  return true;
}

inline
AtomNumber MultiColvarBase::getAbsoluteIndexOfCentralAtom( const unsigned& iatom ) const {
  plumed_dbg_assert( iatom<atom_lab.size() );
  if( atom_lab[iatom].first>0  ) {
    unsigned mmc=atom_lab[iatom].first - 1;
    return mybasemulticolvars[mmc]->getAbsoluteIndexOfCentralAtom( atom_lab[iatom].second );
  }
  plumed_dbg_assert( usespecies );
  return ActionAtomistic::getAbsoluteIndex( atom_lab[getTaskCode(iatom)].second );
}

inline
Vector MultiColvarBase::getPositionOfAtomForLinkCells( const unsigned& iatom ) const {
  plumed_dbg_assert( iatom<atom_lab.size() );
  if( atom_lab[iatom].first>0  ) {
    unsigned mmc=atom_lab[iatom].first - 1;
    return mybasemulticolvars[mmc]->getCentralAtomPos( atom_lab[iatom].second );
  }
  return ActionAtomistic::getPosition( atom_lab[iatom].second );
}

inline
unsigned MultiColvarBase::getNumberOfDerivatives() {
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
