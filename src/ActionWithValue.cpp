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
  numericalDerivatives(false)
{
  if( keywords.exists("NUMERICAL_DERIVATIVES") ) parseFlag("NUMERICAL_DERIVATIVES",numericalDerivatives);
  if(numericalDerivatives) log.printf("  using numerical derivatives\n");
}

ActionWithValue::~ActionWithValue(){
  for(unsigned i=0;i<values.size();++i)delete values[i];
}

void ActionWithValue::clearInputForces(){
  for(unsigned i=0;i<values.size();i++) values[i]->clearInputForce();
}
void ActionWithValue::clearDerivatives(){
  for(unsigned i=0;i<values.size();i++) values[i]->clearDerivatives();
} 

// -- These are the routine for copying the value pointers to other classes -- //

bool ActionWithValue::exists( const std::string& name ) const {
  for(unsigned i=0;i<values.size();++i){
     if (values[i]->name==name) return true;
  }
  return false;
}

Value* ActionWithValue::copyOutput( const std::string& name ) const {
  for(unsigned i=0;i<values.size();++i){
     if (values[i]->name==name) return values[i];
  }
  plumed_massert(0,"there is no pointer with name " + name);
}

Value* ActionWithValue::copyOutput( const unsigned& n ) const {
  plumed_massert(n<values.size(),"you have requested a pointer that is out of bounds");
  return values[n];
}

// -- HERE WE HAVE THE STUFF FOR THE DEFAULT VALUE -- //

void ActionWithValue::addValue(){
  plumed_massert(values.size()==0,"You have already added the default value for this action");
  values.push_back(new Value(getLabel(), false ) );
}

void ActionWithValue::addValueWithDerivatives(){
  plumed_massert(values.size()==0,"You have already added the default value for this action");
  values.push_back(new Value(getLabel(), true ) );
}

void ActionWithValue::setNotPeriodic(){
  plumed_massert(values.size()==1,"The number of components is not equal to one");
  plumed_massert(values[0]->name==getLabel(), "The value you are trying to set is not the default");
  values[0]->min=0; values[0]->max=0;
  values[0]->setupPeriodicity();
}

void ActionWithValue::setPeriodic( const double min, const double max ){
  plumed_massert(values.size()==1,"The number of components is not equal to one");
  plumed_massert(values[0]->name==getLabel(), "The value you are trying to set is not the default");
  values[0]->min=min; values[0]->max=max;
  values[0]->setupPeriodicity();
}

Value* ActionWithValue::getPntrToValue(){
  plumed_massert(values.size()==1,"The number of components is not equal to one");
  plumed_massert(values[0]->name==getLabel(), "The value you are trying to retrieve is not the default");
  return values[0];
}

// -- HERE WE HAVE THE STUFF FOR NAMED VALUES / COMPONENTS -- //

void ActionWithValue::addComponent( const std::string& name ){
  std::string thename; thename=getLabel() + "." + name;
  for(unsigned i=0;i<values.size();++i){
     plumed_massert(values[i]->name!=getLabel(),"Cannot mix single values with components");
     plumed_massert(values[i]->name!=thename,"there is already a value with this name");
  }
  values.push_back(new Value(thename, false ) );
}

void ActionWithValue::addComponentWithDerivatives( const std::string& name ){
  std::string thename; thename=getLabel() + "." + name;
  for(unsigned i=0;i<values.size();++i){
     plumed_massert(values[i]->name!=getLabel(),"Cannot mix single values with components");
     plumed_massert(values[i]->name!=thename,"there is already a value with this name");
  }
  values.push_back(new Value(thename, true ) );
}

int ActionWithValue::getComponent( const std::string& name ) const {
  plumed_massert( !exists( getLabel() ), "You should not be calling this routine if you are using a value");
  plumed_massert(name!=getLabel(),"You should never be calling this routine to retrieve the value");
  std::string thename; thename=getLabel() + "." + name;
  for(unsigned i=0;i<values.size();++i){
     if (values[i]->name==thename) return i;
  }
  plumed_massert(0,"there is no component with name" + name);
  return -1;
}

void ActionWithValue::componentIsNotPeriodic( const std::string& name ){
  int kk=getComponent(name);
  values[kk]->min=0; values[kk]->max=0;
  values[kk]->setupPeriodicity();
}

void ActionWithValue::setGradientsIfNeeded(){
  if(isOptionOn("GRADIENTS")) {
 	 for(unsigned i=0;i<values.size();i++) values[i]->setGradients();
  }
}

void ActionWithValue::componentIsPeriodic( const std::string& name, const double min, const double max ){
  int kk=getComponent(name);
  values[kk]->min=min; values[kk]->max=max;
  values[kk]->setupPeriodicity();
}

Value* ActionWithValue::getPntrToComponent( const std::string& name ){
  int kk=getComponent(name);
  return values[kk]; 
}

Value* ActionWithValue::getPntrToComponent( int n ){
  plumed_massert(n<values.size(),"you have requested a pointer that is out of bounds");
  return values[n];
}
