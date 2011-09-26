#ifndef __PLUMED_ActionWithArguments_h
#define __PLUMED_ActionWithArguments_h

#include "ActionWithValue.h"
#include "Value.h"
#include <vector>

namespace PLMD{

class Value;

/// Action which takes other Action's as arguments.
/// Arguments are objects of type PLMD::Value, and
/// are addressed using the ARG= keyword on the directive line
class ActionWithArguments : public ActionWithValue {
private:
  std::vector<Value*> arguments;
  bool lockRequestArguments;
protected:
/// Reads in the input to action with arguments
  void readActionWithArguments( const std::vector<double>& domain );
/// Applies forces to the arguments 
  void applyForces( const std::vector<double>& forces );
/// Returns the number of derivatives in a particular argument
  unsigned getNumberOfDerivatives( const unsigned& n ) const ;
/// Prints the arguments names out on a file
  void printArgumentNames( FILE* fp );
/// Returns the value of an argument
  double getArgument(int) const;
/// Returns the jth derivative of the ith argument
  double getArgumentDerivative( const unsigned& i, const unsigned& j ) const ;
/// Returns the number of arguments
  unsigned getNumberOfArguments()const;
/// Takes the difference taking into account pbc for arg i
  double difference(int,double,double)const;
public:
  ActionWithArguments(const ActionOptions&);
  virtual ~ActionWithArguments(){};
/// Calculate numerical derivatives
  void calculateNumericalDerivatives();
/// Parse a list of arguments
  void lockRequests();
  void unlockRequests();
};

inline
unsigned ActionWithArguments::getNumberOfDerivatives( const unsigned& n ) const {
  assert( n<arguments.size() );
  return arguments[n]->derivatives.size();
}

inline
double ActionWithArguments::getArgument(int i)const{
  assert( i<arguments.size() );
  return arguments[i]->value;
}

inline
double ActionWithArguments::getArgumentDerivative( const unsigned& i, const unsigned& j ) const {
  assert( i<arguments.size() ); 
  return arguments[i]->getDerivative(j);
}

inline
unsigned ActionWithArguments::getNumberOfArguments()const{
  return arguments.size();
}

inline
void ActionWithArguments::applyForces( const std::vector<double>& forces ){
  assert( forces.size()==arguments.size() );
  for(unsigned i=0;i<arguments.size();++i){ arguments[i]->inputForce+=forces[i]; }
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
