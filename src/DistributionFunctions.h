#ifndef __PLUMED_DistributionFunction_h
#define __PLUMED_DistributionFunction_h

#include <string>
#include <sstream> 
#include <cstring>
#include <cassert>
#include <vector>
#include "Value.h"
#include "Tools.h"
#include "SwitchingFunction.h"

namespace PLMD{

//+DEVELDOC TOOLBOX DistributionFunction
/**
DistributionFunction is an abstract base class.  The classes that inherit
from it can be used to calculate functions of a distribution of values such
as the number of values less than a target, the minimum, the average and so 
on.  This class is used in PLMD::ActionWithDistribution.  
*/
//+ENDDEVELDOC

class DistributionFunction {
private:
  bool fine;
  std::string errorMsg;
public:
  DistributionFunction( const std::vector<std::string>& parameters );
  double calculate( Value* value_in, Value* value_out );
  void finish( const double& total, Value* value_out );
  void error(const std::string& msg);
  bool check() const;
  std::string errorMessage() const ;
  virtual double compute( const double val, double& df )=0;
  virtual double last_step( const double total, double& df )=0;
  virtual std::string message()=0;
}; 

inline 
void DistributionFunction::error( const std::string& msg ){
  fine=false; errorMsg=msg;
}

inline
std::string DistributionFunction::errorMessage() const {
  return errorMsg;
}

inline
bool DistributionFunction::check() const {
  return fine;
}

class less_than : public DistributionFunction {
private:
  double r_0;
  unsigned nn,mm;
  SwitchingFunction sf;
public:
  less_than( const std::vector<std::string>& parameters );
  double compute( const double val, double& df );
  double last_step( const double total, double& df );
  std::string message();
};

}

#endif
