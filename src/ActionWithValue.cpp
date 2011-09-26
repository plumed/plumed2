#include "ActionWithValue.h"
#include "PlumedMain.h"

using namespace std;
using namespace PLMD;

ActionWithValue::ActionWithValue(const ActionOptions&ao):
  Action(ao),
  numericalDerivatives(false)
{
  registerKeyword(0, "NUMERICAL_DERIVATIVES", "calculate the derivatives for these quantities numerically"); 
}

ActionWithValue::~ActionWithValue(){
  for(unsigned i=0;i<values.size();++i) delete values[i];
}

void ActionWithValue::readActionWithValue( const unsigned& nd, const std::vector<double>& d ){
  parseFlag("NUMERICAL_DERIVATIVES",numericalDerivatives);
  if (numericalDerivatives) log.printf("  using numerical derivatives\n");
  domain.resize(2); domain[0]=d[0]; domain[1]=d[1]; nderivatives=nd;
}

void ActionWithValue::noAnalyticalDerivatives(){
  numericalDerivatives=true;
  warning("Numerical derivatives will be used as analytical derivatives are not available");
}

void ActionWithValue::addValue( const std::string& name, const bool& ignorePeriod, const bool& hasDerivatives ){
  std::string thename; thename=getLabel() + "." + name;
  for(unsigned i=0;i<values.size();++i){
     if( values[i]->myname==thename ) assert(false);
  }
  
  unsigned nder;
  if( hasDerivatives ){ nder=nderivatives; } else{ nder=0; }
  if( ignorePeriod ){
      std::vector<double> fdomain(2,0.0);
      values.push_back(new Value(*this, thename, nder, fdomain) );
  } else {
      values.push_back(new Value(*this, thename, nder, domain) );
  } 
}

void ActionWithValue::clearInputForces(){
  for(unsigned i=0;i<values.size();i++){
     values[i]->hasForce=false; values[i]->inputForce=0.0;
  }
}

void ActionWithValue::clearDerivatives(){
  for(unsigned i=0;i<values.size();i++) values[i]->clearDerivatives();
}

Value* ActionWithValue::getValuePointer( const unsigned& i ){ 
  assert( i<values.size() );
  return values[i]; 
}

Value* ActionWithValue::getValuePointer( const std::string& name ){
  std::string thename=getLabel() + "." +  name;
  for(unsigned i=0;i<values.size();++i){
     if( values[i]->myname==thename ) return values[i];
  }
  error("there is no component with label " + name );
  return values[0];
}

/*void ActionWithValue::check(const std::string&name){
  assertUnique(name);
  if(name==""){
    hasUnnamedValue=true;
    assert(!hasMultipleValues);
  }else{
    hasMultipleValues=true;
    assert(!hasUnnamedValue);
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
} */


