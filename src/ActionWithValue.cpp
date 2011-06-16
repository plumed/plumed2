#include "ActionWithValue.h"
#include "PlumedMain.h"

using namespace std;
using namespace PLMD;



ActionWithValue::ActionWithValue(const ActionOptions&ao):
  Action(ao),
  numberOfParameters(0)
{}

ActionWithValue::~ActionWithValue(){
  for(unsigned i=0;i<values.size();++i)delete values[i];
}

void ActionWithValue::addValue(const std::string&name){
  assertUnique(name);
  values.push_back(new Value(*this,name));
}

void ActionWithValue::addValueWithDerivatives(const std::string&name){
  assertUnique(name);
  Value* v=new Value(*this,name);
  v->enableDerivatives();
  values.push_back(v);
}

bool ActionWithValue::hasNamedValue(const std::string&name)const{
  for(unsigned i=0;i<values.size();++i){
    if(name==values[i]->getName()) return true;
   }
  return false;
}

int ActionWithValue::getValueIndex(const std::string&name)const{
  for(unsigned i=0;i<values.size();++i) if(name==values[i]->getName()) return i;
  assert(0);
  return -1; // otherwise the compiler complains
}

Value* ActionWithValue::getValue(const std::string&name)const{
  return values[getValueIndex(name)];
}

Value* ActionWithValue::getValue(int i)const{
  return values[i];
}

std::vector<std::string> ActionWithValue::getValueNames()const{
  std::vector<std::string> ret;
  for(unsigned i=0;i<values.size();i++) ret.push_back(values[i]->getName());
  return ret;
}

void ActionWithValue::setNumberOfParameters(int n){
  numberOfParameters=n;
  for(unsigned i=0;i<values.size();i++) values[i]->setNumberOfParameters(n);
}

void ActionWithValue::clearInputForces(){
  for(unsigned i=0;i<values.size();i++) values[i]->clearInputForce();
}
void ActionWithValue::clearDerivatives(){
  for(unsigned i=0;i<values.size();i++) values[i]->clearDerivatives();
}


