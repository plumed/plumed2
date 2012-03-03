#ifndef __PLUMED_DistributionFunction_h
#define __PLUMED_DistributionFunction_h

#include <string>
#include <sstream> 
#include <cstring>
#include <cassert>
#include <vector>
#include <iostream>
#include "Value.h"
#include "Tools.h"
#include "SwitchingFunction.h"
#include "HistogramBead.h"

namespace PLMD{

//+DEVELDOC TOOLBOX DistributionFunction
/**
DistributionFunction is an abstract base class.  The classes that inherit
from it can be used to calculate functions of a distribution of values such
as the number of values less than a target, the minimum, the average and so 
on.  This class is used in PLMD::ActionWithDistribution.  
*/
//+ENDDEVELDOC

class DistributionFunctionDocs {
/// The documentation for this particular thing
  std::map<std::string,std::string> alldocs;
public:
  DistributionFunctionDocs();
/// Print the documentation for the object named 
  void printDocs(const std::string& name);
};

class DistributionFunction {
private:
  bool fine;
  double last_add;
  std::string errorMsg;
protected:
  void copyDerivatives( Value* value_in, Value* value_out );
public:
  DistributionFunction( const std::vector<std::string>& parameters );
  virtual void calculate( Value* value_in, std::vector<Value>& aux, Value* value_out )=0;
  virtual void finish( const double& total, Value* value_out )=0;
  virtual bool sizableContribution(Value* val, const double& tol);
  virtual void gather( PlumedCommunicator& comm ){};
  void error(const std::string& msg);
  bool check() const;
  std::string errorMessage() const ;
  virtual std::string message()=0;
}; 

inline
bool DistributionFunction::sizableContribution( Value* val, const double& tol ){
  return (val->get()>=tol);
}

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

class sum : public DistributionFunction {
private:
  double nval;
  double prev_nval;
public:
  static void writeDocs( std::string& docs );
  sum( const std::vector<std::string>& parameters );
  void calculate( Value* value_in, std::vector<Value>& aux, Value* value_out );
  inline bool sizableContribution(Value* val, const double& tol){ return true; }
  void gather( PlumedCommunicator& comm );
  void finish( const double& total, Value* value_out );
  std::string message();
};

class mean : public DistributionFunction {
private:
  double nval;
  double prev_nval;
public:
  static void writeDocs( std::string& docs );
  mean( const std::vector<std::string>& parameters );
  void calculate( Value* value_in, std::vector<Value>& aux, Value* value_out );
  void finish( const double& total, Value* value_out );
  inline bool sizableContribution(Value* val, const double& tol){ return true; }
  void gather( PlumedCommunicator& comm );
  std::string message();
};

class less_than : public DistributionFunction {
private:
  double r_0;
  unsigned nn,mm;
  double total, last_add;
  SwitchingFunction sf;
public:
  static void writeDocs( std::string& docs );
  less_than( const std::vector<std::string>& parameters );
  void calculate( Value* value_in, std::vector<Value>& aux, Value* value_out );
  void finish( const double& total, Value* value_out );
  std::string message();
};

class more_than : public DistributionFunction {
private:
  double r_0;
  unsigned nn,mm;
  SwitchingFunction sf;
public:
  static void writeDocs( std::string& docs ); 
  more_than( const std::vector<std::string>& parameters );
  void calculate( Value* value_in, std::vector<Value>& aux, Value* value_out );
  void finish( const double& total, Value* value_out );
  std::string message();
};

class within : public DistributionFunction {
private:
  double a,b,sigma;
  HistogramBead hist;
public:
  static void writeDocs( std::string& docs );
  within( const std::vector<std::string>& parameters );
  void calculate( Value* value_in, std::vector<Value>& aux, Value* value_out );
  void finish( const double& total, Value* value_out );
  std::string message();
};

class min : public DistributionFunction {
private:
  double beta;
public:
  static void writeDocs( std::string& docs );
  min( const std::vector<std::string>& parameters );
  void calculate( Value* value_in, std::vector<Value>& aux, Value* value_out );
  void finish( const double& total, Value* value_out );
  std::string message();
};

}

#endif
