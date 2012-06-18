#ifndef __PLUMED_ActionWithArguments_h
#define __PLUMED_ActionWithArguments_h

#include "Action.h"
#include "Value.h"
#include <vector>

namespace PLMD{

/**
\ingroup MULTIINHERIT
This is used to create PLMD::Action objects that take the output from some other Action as input.  
This is used in PLMD::Function and PLMD::Bias
 PLMD::Action objects that inherit from PLMD::ActionWithArguments take 
 values and components calculated in other PLMD::Action objects and
 use this information to calculate some new function.  If you have 
 only one list of arguments you should use the reserved keyword <b> ARG </b> 
 when you use parseArgumentList.
*/

class ActionWithArguments:
  public virtual Action
{
  std::vector<Value*> arguments;
  bool lockRequestArguments;
public:
/// Get the scalar product between the gradients of two variables
  double getProjection(unsigned i,unsigned j)const;
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
/// Takes one value and brings it back into the pbc of argument i 
  double bringBackInPbc(int i,double d1)const;
/// Parse a list of arguments
  void parseArgumentList(const std::string&key,std::vector<Value*>&args);
/// Setup the dependencies
  void requestArguments(const std::vector<Value*> &arg);
public:
  ActionWithArguments(const ActionOptions&);
  virtual ~ActionWithArguments(){};
/// Calculate the numerical derivatives
/// N.B. only pass an ActionWithValue to this routine if you know exactly what you 
/// are doing.  The default will be correct for the vast majority of cases
  void calculateNumericalDerivatives( ActionWithValue* a=NULL );
  void lockRequests();
  void unlockRequests();
/// Returns an array of pointers to the arguments
  std::vector<Value*>    & getArguments();
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
double ActionWithArguments::bringBackInPbc(int i,double d1)const{
  return arguments[i]->bringBackInPbc(d1);
}

inline
void ActionWithArguments::lockRequests(){
  lockRequestArguments=true;
}

inline
void ActionWithArguments::unlockRequests(){
  lockRequestArguments=false;
}

inline
std::vector<Value*> & ActionWithArguments::getArguments(){
  return arguments;
}

}

#endif
