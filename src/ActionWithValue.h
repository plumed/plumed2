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
  void addValue(const std::string&name);
  void addValueWithDerivatives(const std::string&name);
  bool hasNamedValue(const std::string&name)const;
  Value* getValue(const std::string&name)const;
  Value* getValue(int i)const;
  std::vector<std::string> getValueNames()const;
  int getNumberOfValues();
  void setNumberOfParameters(int n);
  int getNumberOfParameters()const;
  double getForce(int n);
  void addInputForces(int i,double f);
  void clearInputForces();
  void clearDerivatives();
  void setValue(Value*,double);
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
