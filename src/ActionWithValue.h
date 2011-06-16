#ifndef __PLUMED_ActionWithValue_h
#define __PLUMED_ActionWithValue_h

#include "Action.h"
#include <map>
#include <vector>
#include <cassert>
#include "Value.h"

namespace PLMD{

/// Action which can take one or more values.
/// This object contains an array of PLMD::Value, one for each component.
/// It also stores all the derivatives of these values wrt the parameters
/// Parameters are other values (from other Action s) or atomic positions.
class ActionWithValue:
  public virtual Action
{
  int numberOfParameters;
  std::vector<Value*> values;
  void assertUnique(const std::string&name);
  int getValueIndex(const std::string&name)const;
public:
  ActionWithValue(const ActionOptions&ao);
  ~ActionWithValue();
/// Add a new value without derivatives.
/// This should be used for values which are only evaluated (e.g. for printing)
/// but for which we do not make derivatives available so that forces cannot
/// be applied
  void addValue(const std::string&name);
/// Add a new value with derivatives.
/// This should be used for values for which we make derivatives available
/// so that forces can be applied
  void addValueWithDerivatives(const std::string&name);
/// Check if a value with a given name is already used
  bool hasNamedValue(const std::string&name)const;
/// Return a pointer to the value by name
  Value* getValue(const std::string&name)const;
/// Return a pointer to the value by index
/// This should be an indexes growing for new inserted values.
/// E.g., the default value (no name) is number 0, ...
  Value* getValue(int i)const;
/// Returns an array of strings with the names of the values
  std::vector<std::string> getValueNames()const;
/// Returns the number of values defined
  int getNumberOfValues();
/// Set the number of parameters on which this Action depends.
/// Example: for a Bias, this is the number of arguments, for a Colvar
/// is 3*Natoms+cell variables
  void setNumberOfParameters(int n);
/// Returns the number of parameters on which this Action depends.
  int getNumberOfParameters()const;
/// Returns the total force applied on i-th value 
  double getForce(int n);
/// Add a force to the i-th value 
  void addInputForces(int i,double f);
/// Clear the forces on the values
  void clearInputForces();
/// Clear the derivatives of values wrt parameters
  void clearDerivatives();
/// Set the value
  void setValue(Value*,double);
/// Set the default value (the one without name)
  void setValue(double);
};

inline
void ActionWithValue::setValue(Value*v,double d){
  v->setValue(d);
}

inline
void ActionWithValue::setValue(double d){
  values[0]->setValue(d);
}

inline
void ActionWithValue::addInputForces(int i,double f){values[i]->addForce(f);}

inline
double ActionWithValue::getForce(int n){return values[n]->getForce();}

inline
void ActionWithValue::assertUnique(const std::string&name){
  assert(!hasNamedValue(name));
}

inline
int ActionWithValue::getNumberOfValues(){return values.size();}

inline
int ActionWithValue::getNumberOfParameters()const{
  return numberOfParameters;
}



}

#endif
