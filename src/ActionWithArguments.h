#ifndef __PLUMED_ActionWithArguments_h
#define __PLUMED_ActionWithArguments_h

#include "Action.h"
#include "Value.h"
#include <vector>

namespace PLMD{

class Value;

/// Action which takes other Action's as arguments.
/// Arguments are objects of type PLMD::Value, and
/// are addressed using the ARG= keyword on the directive line
class ActionWithArguments:
  public virtual Action
{
  std::vector<Value*> arguments;
  bool lockRequestArguments;

protected:
                           ActionWithArguments(const ActionOptions&);
  virtual                 ~ActionWithArguments(){};
public:
/// Returns an array of pointers to the arguments
  std::vector<Value*>    & getArguments();
/// Returns the value of an argument
  double                   getArgument(int)const;
/// Returns the number of arguments
  unsigned                 getNumberOfArguments()const;
///
  void                     calculateNumericalDerivatives();
/// Takes the difference taking into account pbc for arg i
  double                   difference(int,double,double)const;
/// Parse a list of arguments
  void                     parseArgumentList(const std::string&key,std::vector<Value*>&args);
  void                     requestArguments(const std::vector<Value*> &arg);
  void lockRequests();
  void unlockRequests();
};


inline
std::vector<Value*> & ActionWithArguments::getArguments(){
  return arguments;
}

inline
double ActionWithArguments::getArgument(int i)const{
  return arguments[i]->get();
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
