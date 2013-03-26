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
#ifndef __PLUMED_vesselbase_ActionWithVessel_h
#define __PLUMED_vesselbase_ActionWithVessel_h

#include "core/ActionWithValue.h"
#include "core/ActionAtomistic.h"
#include "tools/Exception.h"
#include "tools/DynamicList.h"
#include <vector>

namespace PLMD{
class Value;

namespace vesselbase{

class Vessel;
class BridgeVessel;

/**
\ingroup MULTIINHERIT
This is used to create PLMD::Action objects that are computed by calculating the same function multiple
times.  This is used in PLMD::MultiColvar.
*/

class ActionWithVessel : public virtual Action {
friend class Vessel;
friend class ShortcutVessel;
friend class StoreValuesVessel;
friend class FunctionVessel;
friend class BridgeVessel;
private:
/// This is used to ensure that we have properly read the action
  bool read;
/// Do all calculations in serial
  bool serial;
/// The tolerance on the accumulators for neighbour list
  double tolerance;
/// The value of the current element in the sum
  std::vector<double> thisval;
/// A boolean that makes sure we don't accumulate very wrong derivatives
  std::vector<bool> thisval_wasset;
/// Vector of derivatives for the object
  std::vector<double> derivatives;
/// The buffers we use for mpi summing DistributionFunction objects
  unsigned current_buffer_start;
  unsigned current_buffer_stride;
  std::vector<double> buffer;
/// Pointers to the functions we are using on each value
  std::vector<Vessel*> functions;
/// Tempory storage for forces
  std::vector<double> tmpforces;
/// Avoid hiding of base class virtual function
  using Action::deactivate;
protected:
/// Does the weight have derivatives
  bool weightHasDerivatives;
/// The numerical index of the task we are curently performing
  unsigned current;
/// The list of tasks we have to perform
  DynamicList<unsigned> taskList;
/// Add a vessel to the list of vessels
  void addVessel( const std::string& name, const std::string& input, const int numlab=0, const std::string thislab="" );
  void addVessel( Vessel* vv );
/// Add a bridging vessel to the list of vessels
  void addBridgingVessel( ActionWithVessel* tome, BridgeVessel* bv );
/// Complete the setup of this object (this routine must be called after construction of ActionWithValue)
  void readVesselKeywords();
/// Return the value of the tolerance
  double getTolerance() const ;
/// Get the number of vessels
  unsigned getNumberOfVessels() const;
/// Get a pointer to the ith vessel
   Vessel* getPntrToVessel( const unsigned& i );
/// Calculate the values of all the vessels
  void runAllTasks();
/// Deactivate the task we are currently working on
  void deactivateCurrentTask();
/// Finish running all the calculations
  void finishComputations();
/// Resize all the functions when the number of derivatives change
  void resizeFunctions();
/// Set the derivative of the jth element wrt to a numbered element
  void setElementDerivative( const unsigned&, const double& );
/// This loops over all the vessels calculating them and also 
/// sets all the element derivatives equal to zero
  bool calculateAllVessels();
/// Retrieve the forces from all the vessels (used in apply)
  void getForcesFromVessels( std::vector<double>& forcesToApply );
/// This is used to accumulate the derivatives when we merge using chainRuleForElementDerivatives
  void accumulateDerivative( const unsigned& ider, const double& df );
/// Clear tempory data that is calculated for each task
  void clearAfterTask();
public:
  static void registerKeywords(Keywords& keys);
  ActionWithVessel(const ActionOptions&ao);
  ~ActionWithVessel();
/// Activate the jth colvar
/// Deactivate the current task in future loops
  virtual void deactivate_task()=0;
/// Merge the derivatives
  void chainRuleForElementDerivatives( const unsigned&, const unsigned&, const double& , Vessel* );
  void chainRuleForElementDerivatives( const unsigned&, const unsigned& , const unsigned& , const unsigned& , const double& , Vessel* );
  virtual void mergeDerivatives( const unsigned& ider, const double& df );
  virtual unsigned getFirstDerivativeToMerge();
  virtual unsigned getNextDerivativeToMerge( const unsigned& );
  virtual void buildDerivativeIndexArrays( std::vector< DynamicList<unsigned> >& active_der );
  virtual void clearDerivativesAfterTask( const unsigned& );
/// Are the base quantities periodic
  virtual bool isPeriodic()=0;
/// What are the domains of the base quantities
  virtual void retrieveDomain( std::string& min, std::string& max);
/// Get the number of functions from which we are calculating the distribtuion
  virtual unsigned getNumberOfFunctionsInAction()=0;
/// Get the number of derivatives for final calculated quantity 
  virtual unsigned getNumberOfDerivatives()=0;
/// Do any jobs that are required before the task list is undertaken
  virtual void doJobsRequiredBeforeTaskList();
/// Calculate one of the functions in the distribution
  virtual bool performTask( const unsigned& j )=0;
/// Return a pointer to the field 
  Vessel* getVessel( const std::string& name );
///  Add some derivative of the quantity in the sum wrt to a numbered element
  void addElementDerivative( const unsigned&, const double& );
/// Set the value of the element
  void setElementValue( const unsigned& , const double& );
/// Get the value of this element
  double getElementValue( const unsigned& ival ) const ;
/// Retrieve the derivative of the quantity in the sum wrt to a numbered element
  double getElementDerivative( const unsigned& ) const ;
};

inline
double ActionWithVessel::getTolerance() const {
  return tolerance;
}

inline
unsigned ActionWithVessel::getNumberOfVessels() const {
  return functions.size();
}

inline
Vessel* ActionWithVessel::getPntrToVessel( const unsigned& i ){
  plumed_dbg_assert( i<functions.size() );
  return functions[i];
}

inline
double ActionWithVessel::getElementValue(const unsigned& ival) const {
  return thisval[ival];
}

inline
void ActionWithVessel::setElementValue( const unsigned& ival, const double& val ){
  // Element 0 is reserved for the value we are accumulating
  // Element 1 is reserved for the normalization constant for calculating AVERAGES, normalized HISTOGRAMS
  plumed_dbg_massert( !thisval_wasset[ival], "In action named " + getName() + " with label " + getLabel() );
  thisval[ival]=val;
  thisval_wasset[ival]=true;
}

inline
unsigned ActionWithVessel::getFirstDerivativeToMerge(){
  return 0;
}

inline
unsigned ActionWithVessel::getNextDerivativeToMerge( const unsigned& ider){
  return ider+1;
}

inline
double ActionWithVessel::getElementDerivative( const unsigned& ider ) const {
  plumed_dbg_assert( ider<derivatives.size() );
  return derivatives[ider];
}

inline
void ActionWithVessel::addElementDerivative( const unsigned& ider, const double& der ){
#ifndef NDEBUG
  unsigned ndertmp=getNumberOfDerivatives();
  if( ider>=ndertmp && ider<2*ndertmp ) plumed_dbg_massert( weightHasDerivatives, "In " + getLabel() );
#endif
  plumed_dbg_assert( ider<derivatives.size() );
  derivatives[ider] += der;
}

inline
void ActionWithVessel::setElementDerivative( const unsigned& ider, const double& der ){
  plumed_dbg_assert( ider<derivatives.size() );
  derivatives[ider] = der;
}

inline
void ActionWithVessel::accumulateDerivative( const unsigned& ider, const double& der ){
  plumed_dbg_assert( ider<getNumberOfDerivatives() );
  buffer[current_buffer_start + current_buffer_stride*ider] += der;
}


} 
}
#endif
