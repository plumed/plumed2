#ifndef __PLUMED_ActionWithValue_h
#define __PLUMED_ActionWithValue_h

#include "Action.h"
#include "Value.h"
#include <vector>
#include <cassert>

namespace PLMD{

/// Action which can take one or more values.
/// This object contains an array of PLMD::Value, one for each component.
/// It also stores all the derivatives of these values wrt the parameters
/// Parameters are other values (from other Action s) or atomic positions.
class ActionWithValue : public Action {
  int numberOfParameters;
  std::vector<Value*> values;
  void assertUnique(const std::string&name);
  void check(const std::string&name);
  int getValueIndex(const std::string&name)const;
  bool numericalDerivatives;
  bool hasMultipleValues;
  bool hasUnnamedValue;
protected:
/// Enforce the use of numerical derivatives.
/// This may be useful during the implementation of new collective
/// variables. Before implementing the derivatives, the used can
/// just tell plumed to use finite difference irrespectively of
/// the NUMERICAL_DERIVATIVES keyword in the input file
  void enforceNumericalDerivatives();
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
/// Clear the forces on the values
  void clearInputForces();
/// Clear the derivatives of values wrt parameters
  void clearDerivatives();
/// Set the value
  void setValue(Value*,double);
/// Set the default value (the one without name)
  void setValue(double);
/// Check if numerical derivatives should be used
  bool checkNumericalDerivatives()const;
};

inline
void ActionWithValue::setValue(Value*v,double d){
  v->set(d);
}

inline
void ActionWithValue::setValue(double d){
  values[0]->set(d);
}

inline
double ActionWithValue::getForce(int n){
  return values[n]->getForce();
}

inline
void ActionWithValue::assertUnique(const std::string&name){
  assert(!hasNamedValue(name));
}

inline
int ActionWithValue::getNumberOfValues(){
  return values.size();
}

inline
int ActionWithValue::getNumberOfParameters()const{
  return numberOfParameters;
}

inline
bool ActionWithValue::checkNumericalDerivatives()const{
  return numericalDerivatives;
}


}

#endif
