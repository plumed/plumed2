#ifndef __PLUMED_ActionWithDistribution_h
#define __PLUMED_ActionWithDistribution_h

#include "ActionWithValue.h"
#include "ActionAtomistic.h"
#include "Value.h"
#include "PlumedException.h"
#include "DistributionFunctions.h"
#include "DynamicList.h"
#include <vector>

namespace PLMD{

class ActionWithDistribution : public virtual Action {
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
/// Accumulators for the values
  std::vector<double> totals;
/// Pointers to the values for this object
  std::vector<Value*> final_values;
/// Pointers to the functions we are using on each value
  std::vector<DistributionFunction*> functions;
/// The list of quantities that should be calculated
  DynamicList members;
protected:
/// Add a distribution function to the list
  void addDistributionFunction( std::string name, DistributionFunction* fun );
/// Read the keywords for the distribution (this routine must be called after construction of ActionWithValue)
  void readDistributionKeywords();
/// Find out if we are running the calculation without mpi
  bool getSerial() const ;
/// This resets members so that we calculate all functions - this is used for neighbour list update
//  void resetMembers();
/// Find out if it is time to do neighbor list update
  bool isTimeForNeighborListUpdate() const ;
/// Get the frequency with which to update neighbor lists
  int getUpdateFreq() const ;
public:
  static void registerKeywords(Keywords& keys);
/// By calling this function during register keywords you tell plumed to use a
/// a default method to parallelize the calculation.
  static void autoParallelize(Keywords& keys);
  ActionWithDistribution(const ActionOptions&ao);
  ~ActionWithDistribution();
/// Prepare everything for neighbour list update
  virtual void prepare();
/// Calculate the values of the object
  void calculate();
/// Are we using distributions 
  bool usingDistributionFunctions() const;
/// Overwrite this in your inherited actions if neighbour list update is more complicated
/// than just calculating everything and seeing whats big.
  virtual void prepareForNeighborListUpdate(){};
/// Overwrite this in your inherited actions if neighbour list update is more complicated
/// than just calculating everything and seeing whats big.
  virtual void completeNeighborListUpdate(){};
/// Merge the derivatives
  virtual void mergeDerivatives( const unsigned j, Value* value_in, Value* value_out )=0;
/// Get the number of functions from which we are calculating the distribtuion
  virtual unsigned getNumberOfFunctionsInDistribution()=0;
/// Calculate one of the functions in the distribution
  virtual void calculateThisFunction( const unsigned& j, Value* value_in, std::vector<Value>& aux )=0;
/// Get the number of derivatives for a given function
  virtual unsigned getThisFunctionsNumberOfDerivatives( const unsigned& j )=0;
};

inline
bool ActionWithDistribution::getSerial() const {
  return serial;
}

inline
bool ActionWithDistribution::usingDistributionFunctions() const {
  return !all_values;
}

inline
bool ActionWithDistribution::isTimeForNeighborListUpdate() const {
  return reduceAtNextStep;
}

inline
int ActionWithDistribution::getUpdateFreq() const {
  return updateFreq;
}

}

#endif

