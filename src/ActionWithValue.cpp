#include "ActionWithValue.h"
#include "PlumedMain.h"
#include "PlumedException.h"

using namespace std;
using namespace PLMD;

void ActionWithValue::registerKeywords(Keywords& keys){
  keys.addFlag("NUMERICAL_DERIVATIVES", false, "calculate the derivatives for these quantities numerically");
}

void ActionWithValue::noAnalyticalDerivatives(Keywords& keys){
   keys.remove("NUMERICAL_DERIVATIVES");
   keys.addFlag("NUMERICAL_DERIVATIVES",true,"analytical derivatives are not implemented for this keyword so numerical derivatives are always used");
}

ActionWithValue::ActionWithValue(const ActionOptions&ao):
  Action(ao),
  numberOfParameters(0),
  numericalDerivatives(false),
  hasMultipleValues(false),
  hasUnnamedValue(false)
{
  if( keywords.exists("NUMERICAL_DERIVATIVES") ) parseFlag("NUMERICAL_DERIVATIVES",numericalDerivatives);
  if(numericalDerivatives) log.printf("  using numerical derivatives\n");
}

ActionWithValue::~ActionWithValue(){
  for(unsigned i=0;i<values.size();++i)delete values[i];
}

void ActionWithValue::check(const std::string&name){
  assertUnique(name);
  if(name==""){
    hasUnnamedValue=true;
    plumed_massert(!hasMultipleValues,"cannot use unnamed components for a multicomponent Action");
  }else{
    hasMultipleValues=true;
    plumed_massert(!hasUnnamedValue,"cannot use unnamed components for a multicomponent Action");
  }
}

void ActionWithValue::addValue(const std::string&name){
  check(name);
  values.push_back(new Value(*this,name));
}

void ActionWithValue::addValueWithDerivatives(const std::string&name){
  check(name);
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
  plumed_merror("value not found" + name);
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


