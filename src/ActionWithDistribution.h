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
#ifndef __PLUMED_ActionWithDistribution_h
#define __PLUMED_ActionWithDistribution_h

#include "ActionWithValue.h"
#include "ActionAtomistic.h"
#include "Value.h"
#include "PlumedException.h"
#include "DynamicList.h"
#include <vector>

namespace PLMD{

class Vessel;

/**
\ingroup MULTIINHERIT
This is used to create PLMD::Action objects that are computed by calculating the same function multiple
times.  This is used in PLMD::MultiColvar.
*/

class ActionWithDistribution : public virtual Action {
friend class Vessel;
friend class FieldVessel;
private:
/// This is used to ensure that we have properly read the action
  bool read;
/// This tells us we are calculating all values (not doing anything to the distribution)
  bool all_values;
/// Do all calculations in serial
  bool serial;
/// Everything for controlling the updating of neighbor lists
  int updateFreq;
  unsigned lastUpdate;
  bool reduceAtNextStep;
/// The tolerance on the accumulators for neighbour list
  double tolerance;
/// The buffers we use for mpi summing DistributionFunction objects
  std::vector<double> buffer;
/// Pointers to the functions we are using on each value
  std::vector<Vessel*> functions;
/// The list of quantities that should be calculated
  DynamicList<unsigned> members;
/// Deactivave the jth object in the list
  void deactivate( const unsigned& );
/// Activate all the values in the list
  void activateAll();
protected:
/// Complete the setup of this object (this routine must be called after construction of ActionWithValue)
  void requestDistribution();
/// Return the value of the tolerance
  double getTolerance() const ;
/// Find out if it is time to do neighbor list update
  bool isTimeForNeighborListUpdate() const ;
/// Get the number of members that are currently active
  unsigned getNumberOfActiveMembers() const ;
/// Get the jth active member
  unsigned getActiveMember(const unsigned& ) const ;
/// Update the list of active members
  void updateActiveMembers();
/// Get the number of vessels
  unsigned getNumberOfVessels() const;
/// Get a pointer to the ith vessel
   Vessel* getPntrToVessel( const unsigned& i );
/// Calculate the values of all the vessels
  void calculateAllVessels( const int& stepn );
/// Resize all the functions when the number of derivatives change
  void resizeFunctions();
public:
  static void registerKeywords(Keywords& keys);
/// By calling this function during register keywords you tell plumed to use a
/// a default method to parallelize the calculation.
  static void autoParallelize(Keywords& keys);
  ActionWithDistribution(const ActionOptions&ao);
  ~ActionWithDistribution();
/// Activate the jth colvar
  virtual void activateValue( const unsigned j )=0;
/// Ensure that nothing gets done for your deactivated colvars
  virtual void deactivateValue( const unsigned j )=0;
/// Merge the derivatives
  virtual void mergeDerivatives( const unsigned j, const Value& value_in, const double& df, Value& value_out )=0;
/// Can we skip the calculations of quantities
  virtual bool isPossibleToSkip(); 
/// Are the base quantities periodic
  virtual bool isPeriodic()=0;
/// What are the domains of the base quantities
  virtual void retrieveDomain( double& min, double& max);
/// Retrieve the previously calculated value and derivatives
  virtual void retreiveLastCalculatedValue( Value& myvalue )=0;
/// Get the number of functions from which we are calculating the distribtuion
  virtual unsigned getNumberOfFunctionsInAction()=0;
/// Get the number of derivatives for final calculated quantity 
  virtual unsigned getNumberOfDerivatives()=0;
/// Get number of derivatives for ith function
  virtual unsigned getNumberOfDerivatives( const unsigned& i )=0;
/// Calculate one of the functions in the distribution
  virtual bool calculateThisFunction( const unsigned& j )=0;
/// Return a pointer to the field 
  Vessel* getVessel( const std::string& name );
};

inline
double ActionWithDistribution::getTolerance() const {
  return tolerance;
}

inline
bool ActionWithDistribution::isTimeForNeighborListUpdate() const {
  return reduceAtNextStep;
}

inline
bool ActionWithDistribution::isPossibleToSkip(){
  return false;
}

inline
unsigned ActionWithDistribution::getNumberOfActiveMembers() const {
  return members.getNumberActive();
}

inline
unsigned ActionWithDistribution::getActiveMember(const unsigned& m ) const {
  plumed_assert( m<members.getNumberActive() );
  return members[m];
}

inline
void ActionWithDistribution::deactivate( const unsigned& m ){
  members.deactivate(m); 
  deactivateValue(m);
}

inline
void ActionWithDistribution::updateActiveMembers(){
  members.mpi_gatherActiveMembers( comm );
}

inline
unsigned ActionWithDistribution::getNumberOfVessels() const {
  return functions.size();
}

inline
Vessel* ActionWithDistribution::getPntrToVessel( const unsigned& i ){
  plumed_assert( i<functions.size() );
  return functions[i];
}

} 
#endif
