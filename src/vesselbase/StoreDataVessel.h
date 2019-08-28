/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2019 The plumed team
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
#ifndef __PLUMED_vesselbase_StoreDataVessel_h
#define __PLUMED_vesselbase_StoreDataVessel_h

#include <string>
#include <cstring>
#include <vector>
#include "Vessel.h"

namespace PLMD {
namespace vesselbase {

/**
\ingroup TOOLBOX
Objects that inherit from FunctionVessel can be used (in tandem with PLMD::vesselbase::ActionWithVessel)
to store values and derivatives for a set of scalars or vectors that are calculated by a
PLMD::vesselbase::ActionWithVessel.  Functions of these stored quantities can then be calculated in a
second step.
*/

class StoreDataVessel : public Vessel {
  friend class Moments;
private:
/// Do the quantities being stored in here need derivatives
  bool hasderiv;
/// What is the maximum number of vectors we are going to
/// have to store when using lowmem option
  unsigned max_lowmem_stash;
/// The size of the vector we are computing
  unsigned vecsize;
/// The amount of data per vector element
  unsigned nspace;
/// The currently active values
//  std::vector<unsigned> active_val;
/// The active derivative elements
  std::vector<unsigned> active_der;
/// The buffer
  std::vector<double> local_buffer;
/// The actions that are going to use the stored data
  std::vector<ActionWithVessel*> userActions;
/// We create a vector of tempory MultiValues here so as to avoid
/// lots of vector resizing
  unsigned tmp_index;
  std::vector<MultiValue> my_tmp_vals;
protected:
/// Is the weight differentiable
  bool weightHasDerivatives();
/// Are we using low mem option
  bool usingLowMem();
/// Finish the setup of the storage object by setting how much
/// data has to be stored
  void completeSetup( const unsigned&, const unsigned& );
/// Return value of nspace
  unsigned getNumberOfDerivativeSpacesPerComponent() const ;
/// Retrieve the values from the underlying ActionWithVessel
  void storeValues( const unsigned&, MultiValue&, std::vector<double>& ) const ;
/// This stores the data we get from the calculation
  void storeDerivatives( const unsigned&, MultiValue& myvals, std::vector<double>&, std::vector<unsigned>& ) const ;
/// Get the ibuf'th local derivative value
  double getLocalDerivative( const unsigned& ibuf );
/// Set the ibuf'th local derivative value
  void setLocalDerivative( const unsigned& ibuf, const double& val );
public:
  static void registerKeywords( Keywords& keys );
  explicit StoreDataVessel( const VesselOptions& );
/// Get the number of values that have been stored
  virtual unsigned getNumberOfStoredValues() const ;
/// Get the index to store a particular index inside
  unsigned getStoreIndex( const unsigned& ) const ;
/// Get the true index of a quantity from the index it is stored in
  unsigned getTrueIndex( const unsigned& ) const ;
/// Recalculate one of the base quantities
  void recalculateStoredQuantity( const unsigned& myelm, MultiValue& myvals );
/// Set a hard cutoff on the weight of an element
  void setHardCutoffOnWeight( const double& mytol );
/// Add an action that uses this data
  void addActionThatUses( ActionWithVessel* actionThatUses );
/// Return the number of components in the vector
  unsigned getNumberOfComponents() const { return vecsize; }
/// Get the values of all the components in the vector
  void retrieveSequentialValue( const unsigned& myelem, const bool& normed, std::vector<double>& values ) const ;
  void retrieveValueWithIndex( const unsigned& myelem, const bool& normed, std::vector<double>& values ) const ;
  double retrieveWeightWithIndex( const unsigned& myelem ) const ;
/// Get the derivatives for one of the components in the vector
  void retrieveDerivatives( const unsigned& myelem, const bool& normed, MultiValue& myvals );
/// Do all resizing of data
  void resize() override;
///
  std::string description() override { return ""; }
/// Get the number of derivatives for the ith value
  unsigned getNumberOfDerivatives( const unsigned& );
/// Get the size of the derivative list
  unsigned getSizeOfDerivativeList() const ;
/// This stores the data when not using lowmem
  void calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_index ) const override;
/// Final step in gathering data
  void finish( const std::vector<double>& buffer ) override;
/// Is a particular stored value active at the present time
  bool storedValueIsActive( const unsigned& iatom ) const ;
/// Set the active values
  void setActiveValsAndDerivatives( const std::vector<unsigned>& der_index );
/// Activate indexes (this is used at end of chain rule)
  virtual void activateIndices( ActionWithVessel* ) {}
/// Forces on vectors should always be applied elsewhere
  bool applyForce(std::vector<double>&) override { return false; }
///  Get the number of data users
  unsigned getNumberOfDataUsers() const ;
/// Get one of the ith data user
  ActionWithVessel* getDataUser( const unsigned& );
/// Set the number of tempory multivalues we need
  void resizeTemporyMultiValues( const unsigned& nvals );
/// Return a tempory multi value - we do this so as to avoid vector resizing
  MultiValue& getTemporyMultiValue( const unsigned& ind );
};

inline
bool StoreDataVessel::weightHasDerivatives() {
  return getAction()->weightHasDerivatives;
}

inline
bool StoreDataVessel::usingLowMem() {
  return getAction()->lowmem;
}

inline
unsigned StoreDataVessel::getNumberOfDerivativeSpacesPerComponent() const {
  return nspace;
}

inline
bool StoreDataVessel::storedValueIsActive( const unsigned& iatom ) const {
  if( !getAction()->taskIsCurrentlyActive( iatom ) ) return false;
  unsigned jatom = getStoreIndex( iatom );
  plumed_dbg_assert( jatom<getNumberOfStoredValues() );
  return local_buffer[jatom*vecsize*nspace]>epsilon;
}

inline
unsigned StoreDataVessel::getSizeOfDerivativeList() const {
  return active_der.size();
}

inline
unsigned StoreDataVessel::getNumberOfStoredValues() const {
  return getAction()->nactive_tasks;
}

inline
unsigned StoreDataVessel::getStoreIndex( const unsigned& ind ) const {
  if( getAction()->nactive_tasks==getAction()->getFullNumberOfTasks() ) return ind;

  // Binary search for required element - faster scaling than sequential search
  unsigned l=0, r=getAction()->nactive_tasks-1;
  for(unsigned i=0; i<getAction()->nactive_tasks; ++i) {
    plumed_assert( l<=r );
    unsigned m = std::floor( (l + r)/2 );
    if( ind==getAction()->indexOfTaskInFullList[m] ) return m;
    else if( getAction()->indexOfTaskInFullList[m]<ind ) l=m+1;
    else if( getAction()->indexOfTaskInFullList[m]>ind ) r=m-1;
  }
  plumed_merror("requested task is not active");
}

inline
unsigned StoreDataVessel::getTrueIndex( const unsigned& ind ) const {
  return getAction()->indexOfTaskInFullList[ind];
}

inline
void StoreDataVessel::recalculateStoredQuantity( const unsigned& myelem, MultiValue& myvals ) {
  getAction()->performTask( myelem, getAction()->getTaskCode(myelem), myvals );
}

inline
unsigned StoreDataVessel::getNumberOfDataUsers() const {
  return userActions.size();
}

inline
ActionWithVessel* StoreDataVessel::getDataUser( const unsigned& idata ) {
  plumed_dbg_assert( idata<userActions.size() ); return userActions[idata];
}

}
}
#endif


