#ifndef __PLUMED_ActionWithDistribution_h
#define __PLUMED_ActionWithDistribution_h

#include "ActionWithValue.h"
#include "ActionAtomistic.h"
#include "Value.h"
#include "Field.h"
#include "PlumedException.h"
#include "DistributionFunctions.h"
#include "DynamicList.h"
#include <vector>

namespace PLMD{

/**
\ingroup MULTIINHERIT
This is used to create PLMD::Action objects that are computed by calculating the same function multiple
times.  This is used in PLMD::MultiColvar.
*/

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
/// The tolerance on the accumulators for neighbour list
  double tolerance;
/// The buffers we use for mpi summing DistributionFunction objects
  std::vector<double> buffer;
/// Pointers to the values for this object
  std::vector<Value*> final_values;
/// Pointers to the functions we are using on each value
  std::vector<DistributionFunction*> functions;
/// The list of quantities that should be calculated
  DynamicList<unsigned> members;
/// This holds everything for the field
   Field* myfield;
/// Calculate if this is a function of the distribution
  void calculateFunctions();
/// Setup stuff for the derivatives of the field
  void setNumberOfFieldDerivatives( const unsigned D );
protected:
/// Add a distribution function to the list (this routine must be called after construction of ActionWithValue)
  void addDistributionFunction( std::string name, DistributionFunction* fun );
/// Setup a field cv
  void addField( std::string key, Field* ff );
/// Complete the setup of this object (this routine must be called after construction of ActionWithValue)
  void requestDistribution();
/// Find out if we are running the calculation without mpi
  bool getSerial() const ;
/// Find out if it is time to do neighbor list update
  bool isTimeForNeighborListUpdate() const ;
/// Get the frequency with which to update neighbor lists
  int getUpdateFreq() const ;
/// Get the jth active member
  unsigned getActiveMember( const unsigned& j ) const ;
/// Get the tolerance 
  double getTolerance() const ;
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
/// Overwrite this in your inherited actions if neighbor list update is more complicated
/// than just calculating everything and seeing whats big.
  virtual void prepareForNeighborListUpdate(){};
/// Overwrite this in your inherited actions if neighbor list update is more complicated
/// than just calculating everything and seeing whats big.
  virtual void completeNeighborListUpdate(){};
/// Ensure that nothing gets done for your deactivated colvars
  virtual void deactivateValue( const unsigned j )=0;
/// Merge the derivatives
  virtual void mergeDerivatives( const unsigned j, Value* value_in, Value* value_out )=0;
/// Are the base quantities periodic
  virtual bool isPeriodic(const unsigned nn)=0;
/// What are the domains of the base quantities
  virtual void retrieveDomain( const unsigned nn, double& min, double& max);
/// Get the number of functions from which we are calculating the distribtuion
  virtual unsigned getNumberOfFunctionsInDistribution()=0;
/// Calculate one of the functions in the distribution
  virtual void calculateThisFunction( const unsigned& j, Value* value_in, std::vector<Value>& aux )=0;
/// Get the number of derivatives for a given function
  virtual unsigned getThisFunctionsNumberOfDerivatives( const unsigned& j )=0;
// --- Abstract functions for field cvs -- //

/// Setup the values for this field cv and set the value of SIGMA
  virtual void derivedFieldSetup( const double sig )=0;
/// Get the number of derivatives of the field
  virtual unsigned getNumberOfFieldDerivatives()=0;
/// Calculate the contribution of a particular value to the field
  virtual void calculateFieldContribution( const unsigned& j, const std::vector<double>& hisp, Value* tmpvalue, Value& tmpstress, std::vector<Value>& tmpder )=0;
//  virtual void mergeFieldDerivatives( const std::vector<double>& der, Value* value_out )=0;
/// Return a pointer to the field 
  Field* getField();
};

inline
Field* ActionWithDistribution::getField() {
  return myfield;
} 

inline
bool ActionWithDistribution::getSerial() const {
  return serial;
}

inline
bool ActionWithDistribution::isTimeForNeighborListUpdate() const {
  return reduceAtNextStep;
}

inline
int ActionWithDistribution::getUpdateFreq() const {
  return updateFreq;
}

inline
unsigned ActionWithDistribution::getActiveMember( const unsigned& j ) const {
  plumed_assert( j<members.getNumberActive() );
  return members[j];
}

inline
double ActionWithDistribution::getTolerance() const {
  return tolerance;
}

}
#endif

