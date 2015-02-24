/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013,2014 The plumed team
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
private:
/// Do the quantities being stored in here need derivatives
  bool hasderiv;
/// What is the maximum number of vectors we are going to 
/// have to store when using lowmem option 
  unsigned max_lowmem_stash;
/// The start point for the data we are storing in this particular object
  unsigned data_start;
/// The size of the vector we are computing 
  unsigned vecsize;
/// The amount of data per vector element
  unsigned nspace;
/// The currently active values 
  std::vector<unsigned> active_val;
/// The active derivative elements
  std::vector<unsigned> active_der;
/// This is a tempory vector that is used to store data
  std::vector<double> fvec;
/// The local derivatives
  std::vector<double> local_derivatives;
/// The final derivatives
  std::vector<double> final_derivatives;
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
  void storeValues( const unsigned& );
/// Set the Task that needs redoing
  void setTaskToRecompute( const unsigned& ivec );
/// Set a component of one of the vectors
  void setComponent( const unsigned& , const unsigned& , const double& );
/// This is the proper chain rule for vectors
  double chainRule( const unsigned&, const unsigned&, const std::vector<double>& );
/// Chain rule the vector and output derivatives to a value
  void chainRule( const unsigned& , const std::vector<double>&, Value* );
/// Get the ibuf'th local derivative value
  double getLocalDerivative( const unsigned& ibuf );
/// Set the ibuf'th local derivative value
  void setLocalDerivative( const unsigned& ibuf, const double& val );
public:
  static void registerKeywords( Keywords& keys );
  StoreDataVessel( const VesselOptions& );
/// Set a hard cutoff on the weight of an element
  void setHardCutoffOnWeight( const double& mytol );
/// Is the hard weight cutoff on
  bool weightCutoffIsOn() const ;
/// Return the number of components in the vector
  unsigned getNumberOfComponents() const { return vecsize; }
/// Do all resizing of data
  virtual void resize();
/// Clear certain data before start of main loop
  virtual void prepare();
/// Get the number of derivatives for the ith value
  unsigned getNumberOfDerivatives( const unsigned& );
/// Get one of the stored indexes
  unsigned getStoredIndex( const unsigned& , const unsigned& );
/// Get a component of the stored vector
  double getComponent( const unsigned& , const unsigned& );
/// Recalculate a vector - used in lowmem mode
  virtual void recompute( const unsigned& , const unsigned& );
/// This reperforms the task in the underlying action
  virtual void performTask( const unsigned& );
/// This reperforms the task
  virtual void finishTask( const unsigned& ){};
/// Chain rule and store output in local array called final_derivatives
/// with vectors this does chain rule for dot products
  void chainRule( const unsigned& , const std::vector<double>& );
/// Get the ider'th final derivative value
  double getFinalDerivative( const unsigned& ider ) const ;
/// This stores the data when not using lowmem
  bool calculate();
/// This stores the data we get from the calculation
  void storeDerivativesLowMem( const unsigned& );
/// This stores the data we get from the calculation
  void storeDerivativesHighMem( const unsigned& );
/// Final step in gathering data
  virtual void finish();
/// Is a particular stored value active at the present time
  bool storedValueIsActive( const unsigned& iatom ); 
/// Activate indexes (this is used at end of chain rule)
  virtual void activateIndices( ActionWithVessel* ){}
/// Get the list of indices that we are storing data for
  virtual void getIndexList( const unsigned& , const unsigned& , const unsigned& , std::vector<unsigned>& );
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
void StoreDataVessel::performTask( const unsigned& ivec ){
  if( usingLowMem() ) getAction()->performTask();
}

inline
double StoreDataVessel::getComponent( const unsigned& ival, const unsigned& jcomp ){
  plumed_dbg_assert( ival<getAction()->getFullNumberOfTasks() && jcomp<vecsize );
  return getBufferElement( ival*(vecsize*nspace) + jcomp*nspace ); 
}

inline
void StoreDataVessel::setComponent( const unsigned& ival, const unsigned& jcomp, const double& val ){
  plumed_dbg_assert( ival<getAction()->getFullNumberOfTasks() && jcomp<vecsize );
  setBufferElement( ival*(vecsize*nspace) + jcomp*nspace, val );
}

inline
unsigned StoreDataVessel::getNumberOfDerivativeSpacesPerComponent() const {
  return nspace;
}

inline
unsigned StoreDataVessel::getNumberOfDerivatives( const unsigned& ival ){
  plumed_dbg_assert( ival<getAction()->getFullNumberOfTasks() );
  return active_der[ival];
}

inline
unsigned StoreDataVessel::getStoredIndex( const unsigned& ival, const unsigned& jindex ){
  plumed_dbg_assert( ival<getAction()->getFullNumberOfTasks() && jindex<active_der[ival] );

  unsigned kder;
  if( getAction()->lowmem ) kder = max_lowmem_stash + ival*getAction()->getNumberOfDerivatives();
  else kder = getAction()->getFullNumberOfTasks() + ival*(nspace-1); 

  return active_der[kder + jindex];
}

inline
double StoreDataVessel::getLocalDerivative( const unsigned& ibuf ){
  plumed_dbg_assert( getAction()->lowmem && ibuf<local_derivatives.size() );
  return local_derivatives[ibuf];
}

inline
double StoreDataVessel::getFinalDerivative( const unsigned& ider ) const {
  plumed_dbg_assert( ider<final_derivatives.size() );
  return final_derivatives[ider];
}

inline
void StoreDataVessel::setLocalDerivative( const unsigned& ibuf, const double& val ){
  plumed_dbg_assert( getAction()->lowmem && ibuf<local_derivatives.size() );
  local_derivatives[ibuf]=val;
}

inline
bool StoreDataVessel::storedValueIsActive( const unsigned& iatom ){
  plumed_dbg_assert( iatom<active_val.size() );
  return (active_val[iatom]==1);
}

}
}
#endif


