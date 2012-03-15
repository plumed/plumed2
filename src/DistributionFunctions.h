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

class ActionWithDistribution;

class DistributionFunction {
private:
  bool fine;
  double last_add;
  std::string errorMsg;
  std::vector<bool> hasDeriv;
  std::vector<Value*> thesevalues;
  std::vector<Value*> accumulators;
protected:
  void addAccumulator( const bool wderiv );
  void copyValue( const unsigned nn, Value* value_in );
  void extractDerivatives( const unsigned nn, Value* value_in );
  void setValue( const unsigned nn, const double f );
  void chainRule( const unsigned nn, const double f );
  void multiplyValue( const unsigned nn, Value* val );
  Value* getPntrToAccumulator( const unsigned nn );
public:
  DistributionFunction( const std::vector<std::string>& parameters );
  ~DistributionFunction();
  void setNumberOfDerivatives( const unsigned nder );
  void reset();
  void mergeDerivatives( const unsigned kk, ActionWithDistribution& action );
  unsigned requiredBufferSpace() const ;
  void copyDataToBuffers( unsigned& bufsize, std::vector<double>& buffers ) const ;
  void retrieveDataFromBuffers( unsigned& bufsize, const std::vector<double>& buffers );
  virtual void calculate( Value* value_in, std::vector<Value>& aux )=0;
  virtual void finish( Value* value_out )=0;
  virtual bool sizableContribution(const double& tol);
  void error(const std::string& msg);
  bool check() const;
  std::string errorMessage() const ;
  virtual std::string message()=0;
}; 

inline
bool DistributionFunction::sizableContribution( const double& tol ){
  return (thesevalues[0]->get()>=tol);
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

inline
Value* DistributionFunction::getPntrToAccumulator( const unsigned nn ){
  plumed_assert( nn<accumulators.size() );
  return accumulators[nn];
}

class sum : public DistributionFunction {
private:
  double nval;
  double prev_nval;
public:
  static void writeDocs( std::string& docs );
  sum( const std::vector<std::string>& parameters );
  void calculate( Value* value_in, std::vector<Value>& aux );
  inline bool sizableContribution( const double& tol ){ return true; }
  void finish( Value* value_out );
  std::string message();
};

class mean : public DistributionFunction {
private:
  double nval;
public:
  static void writeDocs( std::string& docs );
  mean( const std::vector<std::string>& parameters );
  void calculate( Value* value_in, std::vector<Value>& aux );
  void finish( Value* value_out );
  inline bool sizableContribution( const double& tol ){ return true; }
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
  void calculate( Value* value_in, std::vector<Value>& aux );
  void finish( Value* value_out );
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
  void calculate( Value* value_in, std::vector<Value>& aux );
  void finish( Value* value_out );
  std::string message();
};

class within : public DistributionFunction {
private:
  double a,b,sigma;
  HistogramBead hist;
public:
  static void writeDocs( std::string& docs );
  within( const std::vector<std::string>& parameters );
  void calculate( Value* value_in, std::vector<Value>& aux );
  void finish( Value* value_out );
  std::string message();
};

class min : public DistributionFunction {
private:
  double beta;
public:
  static void writeDocs( std::string& docs );
  min( const std::vector<std::string>& parameters );
  void calculate( Value* value_in, std::vector<Value>& aux );
  void finish( Value* value_out );
  std::string message();
};

class cvdens : public DistributionFunction {
private:
  bool isDensity;
  double ax, bx, xsig, ay, by, ysig, az, bz, zsig;
  std::vector<unsigned> dir;
  std::vector<HistogramBead> beads;
public:
  static void writeDocs( std::string& docs );
  cvdens( const std::vector<std::string>& parameters );
  void calculate( Value* value_in, std::vector<Value>& aux );
  void finish( Value* value_out ); 
  std::string message();
};

}

#endif
