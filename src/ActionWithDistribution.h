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
/// This tells us we are using the field cvs
  bool use_field;
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
// ---- STUFF FOR FIELDS ------ //

   // A class for storing the field
   class FieldClass {
   private:
   /// Total number of points
     unsigned npoints;
   /// Number of high dimensionality vectors
     unsigned ndX;
   /// Number of low dimensionality vectors
     unsigned ndx;
   /// Number of doubles for each point in the grid
     unsigned nper;
   /// The sizes of all the base quantities
     std::vector<unsigned> baseq_nder;
   /// The start in baseq_buffer for each function
     std::vector<unsigned> baseq_starts;
   /// Storage space for the input data
     std::vector<double> baseq_buffer;
   /// Storage space for the grid
     std::vector<double> grid_buffer;
   public:
   /// Setup the field
     void setup( const unsigned nfunc, const unsigned np, const unsigned D, const unsigned d );
   /// Set everything for the field equal to zero
     void clear();
   /// Set the sizes of all the base quantity buffers
     void resizeBaseQuantityBuffers( const std::vector<unsigned>& cv_sizes );
   /// Set the derivatives for one of the base quantities
     void setBaseQuantity( const unsigned nn, Value* val );
   /// Gather the values and derivatives for all the base quantities
     void gatherBaseQuantities( PlumedCommunicator& comm );
   /// Extract one of the base quantities
     void extractBaseQuantity( const unsigned nn, Value* val );
   /// Get number of high dimensional derivatives 
     unsigned get_NdX() const ;
   };
   // This holds everything for the field
   FieldClass myfield;

/// Calculate if this is a function of the distribution
  void calculateFunctions();
/// Calculate the field if this is a field
  void calculateField();
protected:
/// Add a distribution function to the list (this routine must be called after construction of ActionWithValue)
  void addDistributionFunction( std::string name, DistributionFunction* fun );
/// Setup a field cv
  void setupField( unsigned ldim );
/// Complete the setup of this object (this routine must be called after construction of ActionWithValue)
  void requestDistribution();
/// Find out if we are running the calculation without mpi
  bool getSerial() const ;
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
/// Overwrite this in your inherited actions if neighbor list update is more complicated
/// than just calculating everything and seeing whats big.
  virtual void prepareForNeighborListUpdate(){};
/// Overwrite this in your inherited actions if neighbor list update is more complicated
/// than just calculating everything and seeing whats big.
  virtual void completeNeighborListUpdate(){};
/// Ensure that nothing gets done for your deactivated colvars
  virtual void deactivate( const unsigned j )=0;
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
/// Set the kkth output of the field 
  virtual void setFieldOutputValue( const unsigned& kk, Value* tmpvalue )=0;
};

inline
bool ActionWithDistribution::getSerial() const {
  return serial;
}

inline
bool ActionWithDistribution::usingDistributionFunctions() const {
  return ( !all_values && !use_field );
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
unsigned ActionWithDistribution::FieldClass::get_NdX() const {
  return ndX;
}

}

#endif

