/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
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
/// What is the maximum number of vectors we are going to 
/// have to store when using lowmem option 
  unsigned max_lowmem_stash;
/// The start point for the data we are storing in this particular object
  unsigned data_start;
/// The size of the vector we are computing 
  unsigned vecsize;
/// The amount of data per vector element
  unsigned nspace;
/// The active derivative elements
  std::vector<unsigned> active_der;
/// This is a tempory vector that is used to store data
  std::vector<double> fvec;
/// The local derivatives
  std::vector<double> local_derivatives;
/// The final derivatives
  std::vector<double> final_derivatives;
/// Chain rule and store output in local array called final_derivatives
/// with vectors this does chain rule for dot products
  void chainRule( const unsigned& , const std::vector<double>& );
/// This is the proper chain rule for vectors
  double chainRule( const unsigned&, const unsigned&, const std::vector<double>& );
protected:
/// Is the weight differentiable
  bool weightHasDerivatives();
/// Are we using low mem option
  bool usingLowMem();
/// Finish the setup of the storage object by setting how much
/// data has to be stored
  void completeSetup( const unsigned& , const unsigned& );
/// Get the number of derivatives for the ith value
  unsigned getNumberOfDerivatives( const unsigned& );
/// Retrieve the values from the underlying ActionWithVessel
  void storeValues( const unsigned& );
/// Get a component of the stored vector
  double getComponent( const unsigned& , const unsigned& );
/// Chain rule the vector and output derivatives to a value
  void chainRule( const unsigned& , const std::vector<double>&, Value* );
/// Calculate derivatives using chain rule and save in std::vector
  void chainRule( const unsigned& ival, const std::vector<double>& df, std::vector<double>& derout );
/// Multiply the vector by a scalar that is a function of the components of the vector
/// (used to normalize vectors)
  void transformComponents( const unsigned& jstore, const double& weight, double& wdf, const std::vector<double>& dfvec );
public:
  static void registerKeywords( Keywords& keys );
  StoreDataVessel( const VesselOptions& );
/// Return the number of components in the vector
  unsigned getNumberOfComponents() const { return vecsize; }
/// Do all resizing of data
  virtual void resize();
/// Clear certain data before start of main loop
  void prepare();
/// Recalculate a vector - used in lowmem mode
  virtual void recompute( const unsigned& , const unsigned& );
/// This reperforms the task in the underlying action
  virtual void performTask( const unsigned& );
/// This reperforms the task
  virtual void finishTask( const unsigned& ){};
/// This stores the data when not using lowmem
  bool calculate();
/// This stores the data we get from the calculation
  void storeDerivativesLowMem( const unsigned& );
/// This stores the data we get from the calculation
  void storeDerivativesHighMem( const unsigned& );
/// Chain rule the vector and output derivatives to an ActionWithVessel element
  void chainRule( const unsigned& , const unsigned&, const std::vector<double>&, ActionWithVessel* );
/// Chain rule the vector and output derivatives to an ActionWithVessel element
  void chainRule( const unsigned& , const unsigned&, const unsigned&, const std::vector<double>&, ActionWithVessel* );
/// This adds the derivatives to the elements
  void chainRuleForComponent( const unsigned& , const unsigned& , const unsigned& , const double& , ActionWithVessel* );
/// Final step in gathering data
  virtual void finish();
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
  plumed_dbg_assert( ival<getAction()->getNumberOfTasks() && jcomp<vecsize );
  return getBufferElement( ival*(vecsize*nspace) + jcomp*nspace ); 
}

inline
unsigned StoreDataVessel::getNumberOfDerivatives( const unsigned& ival ){
  plumed_dbg_assert( ival<getAction()->getNumberOfTasks() );
  return active_der[ival];
}

}
}
#endif


