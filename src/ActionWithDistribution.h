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
/// This resets members so that we calculate all functions - this is used for neighbour list update
  void resetMembers();
public:
  static void registerKeywords(Keywords& keys);
  ActionWithDistribution(const ActionOptions&ao);
  ~ActionWithDistribution();
/// Calculate the values of the object
  void calculate();
/// Are we using distributions 
  bool usingDistributionFunctions() const;
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
bool ActionWithDistribution::usingDistributionFunctions() const {
  return !all_values;
}

}

#endif

