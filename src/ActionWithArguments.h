#ifndef __PLUMED_ActionWithArguments_h
#define __PLUMED_ActionWithArguments_h

#include "Action.h"
#include "Value.h"
#include <vector>

namespace PLMD{

//+DEVELDOC MULTI-INHERIT ActionWithArguments
/**
This is used to create PLMD::Action objects that take the output from some other Action as input.  
This is used in PLMD::Function and PLMD::Bias
*/
//+ENDDEVELDOC

/// Action which takes other Action's as arguments.
/// Arguments are objects of type PLMD::Value, and
/// are addressed using the ARG= keyword on the directive line
class ActionWithArguments:
  public virtual Action
{
  std::vector<Value*> arguments;
  bool lockRequestArguments;

protected:
///
  double getProjection(unsigned i,unsigned j)const;
public:
/// Registers the list of keywords
  static void registerKeywords( Keywords& keys );
/// Returns the value of an argument
  double getArgument( const unsigned n ) const;
/// Return a pointer to specific argument
  Value* getPntrToArgument( const unsigned n );
/// Returns the number of arguments
  unsigned getNumberOfArguments() const ;
/// Takes the difference taking into account pbc for arg i
  double difference(int, double, double) const;
/// Parse a list of arguments
  void parseArgumentList(const std::string&key,std::vector<Value*>&args);
/// Setup the dependencies
  void requestArguments(const std::vector<Value*> &arg);
public:
  ActionWithArguments(const ActionOptions&);
  virtual ~ActionWithArguments(){};
/// Registers the list of keywords
  static void registerKeywords( Keywords& keys );
/// Calculate the numerical derivatives
  void calculateNumericalDerivatives();
  void lockRequests();
  void unlockRequests();
};


inline
Value* ActionWithArguments::getPntrToArgument( const unsigned n ){
  return arguments[n];
}

inline
double ActionWithArguments::getArgument(const unsigned n) const {
  return arguments[n]->get();
}

inline
unsigned ActionWithArguments::getNumberOfArguments()const{
  return arguments.size();
}

inline
double ActionWithArguments::difference(int i,double d1,double d2)const{
  return arguments[i]->difference(d1,d2);
}

inline
void ActionWithArguments::lockRequests(){
  lockRequestArguments=true;
}

inline
void ActionWithArguments::unlockRequests(){
  lockRequestArguments=false;
}

}

#endif
