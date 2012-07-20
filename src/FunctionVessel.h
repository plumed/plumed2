#ifndef __PLUMED_FunctionVessel_h
#define __PLUMED_FunctionVessel_h

#include <string>
#include <cstring>
#include <vector>
#include "VesselValueAccess.h"

namespace PLMD {

class NormedSumVessel : public VesselAccumulator {
private:
  bool donorm;
  Value myvalue, myvalue2;
  Value myweight, myweight2;
protected:
/// We are normalizing the values
  void useNorm();
public:
  NormedSumVessel( const VesselOptions& );
/// This retrieves data from action and calculates the average
  bool calculate( const unsigned& , const double& );
/// This does the final step of the calculation
  void finish( const double& tolerance );
/// This gets the weight
  virtual void getWeight( const unsigned& , Value& )=0;  
/// This gets each value
  virtual void compute( const unsigned& , const unsigned& , Value& )=0;
};

class SumVessel : public VesselAccumulator {
private:
  Value myvalue, myvalue2;
public:
  SumVessel( const VesselOptions& );
/// This retrieves data from action and calculates
  bool calculate( const unsigned& , const double& );
/// Compute the ith component and the derivatives
  virtual double compute( const unsigned& , const double& , double& )=0;
/// This does the final step of the calculation
  void finish( const double& tolerance );
/// Do any final compuations
  virtual double final_computations( const unsigned& , const double& , double& );
};

}
#endif
