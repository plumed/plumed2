/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2015 The plumed team
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
protected:
/// Apply a hard cutoff on the weight
  bool hard_cut;
/// The value of the cutoff on the weight
  double wtol;
/// Is the weight differentiable
  bool weightHasDerivatives();
/// Are we using low mem option
  bool usingLowMem();
/// Finish the setup of the storage object by setting how much
/// data has to be stored
  void completeSetup( const unsigned& , const unsigned& );
/// Return value of nspace
  unsigned getNumberOfDerivativeSpacesPerComponent() const ;
/// Retrieve the values from the underlying ActionWithVessel
   void storeValues( const unsigned& , MultiValue& , std::vector<double>& ) const ;
/// This stores the data we get from the calculation
  void storeDerivatives( const unsigned& , MultiValue& myvals, std::vector<double>&, std::vector<unsigned>& ) const ;
/// Get the ibuf'th local derivative value
  double getLocalDerivative( const unsigned& ibuf );
/// Set the ibuf'th local derivative value
  void setLocalDerivative( const unsigned& ibuf, const double& val );
public:
  static void registerKeywords( Keywords& keys );
  explicit StoreDataVessel( const VesselOptions& );
/// Get the number of values that have been stored
  unsigned getNumberOfStoredValues() const ;
/// Set a hard cutoff on the weight of an element
  void setHardCutoffOnWeight( const double& mytol );
/// Is the hard weight cutoff on
  bool weightCutoffIsOn() const ;
/// Return the number of components in the vector
  unsigned getNumberOfComponents() const { return vecsize; }
/// Get the values of all the components in the vector
  void retrieveValue( const unsigned& myelem, const bool& normed, std::vector<double>& values ) const ;
/// Get the derivatives for one of the components in the vector
  void retrieveDerivatives( const unsigned& myelem, const bool& normed, MultiValue& myvals );
/// Do all resizing of data
  virtual void resize();
/// Clear certain data before start of main loop
//  virtual void prepare();
///
  virtual std::string description(){ return ""; }
/// Get the number of derivatives for the ith value
  unsigned getNumberOfDerivatives( const unsigned& );
/// Get the size of the derivative list
  unsigned getSizeOfDerivativeList() const ;
/// This stores the data when not using lowmem
  virtual bool calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_index ) const ;
/// Build index stores
//  void buildIndexStores( const unsigned& current, MultiValue& myvals, std::vector<unsigned>& val_index, std::vector<unsigned>& der_index ) const ;
/// Final step in gathering data
  virtual void finish( const std::vector<double>& buffer );
/// Is a particular stored value active at the present time
  bool storedValueIsActive( const unsigned& iatom ); 
/// Set the active values
  void setActiveValsAndDerivatives( const std::vector<unsigned>& der_index );
/// Activate indexes (this is used at end of chain rule)
  virtual void activateIndices( ActionWithVessel* ){}
/// Forces on vectors should always be applied elsewhere
  virtual bool applyForce(std::vector<double>&){ return false; }
};

inline
bool StoreDataVessel::weightHasDerivatives(){
  return getAction()->weightHasDerivatives;
}

inline
bool StoreDataVessel::usingLowMem(){
  return getAction()->lowmem;
}

inline
unsigned StoreDataVessel::getNumberOfDerivativeSpacesPerComponent() const {
  return nspace;
}

inline
bool StoreDataVessel::storedValueIsActive( const unsigned& iatom ){
  plumed_dbg_assert( iatom<getAction()->getFullNumberOfTasks() );
  if( !hard_cut ) return true; 
  return local_buffer[iatom*vecsize*nspace]>wtol;   // (active_val[iatom]==1);
}

inline
unsigned StoreDataVessel::getSizeOfDerivativeList() const {
  return active_der.size();
}

inline
unsigned StoreDataVessel::getNumberOfStoredValues() const {
  return getAction()->getFullNumberOfTasks();
}

}
}
#endif


