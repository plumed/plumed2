#include "Value.h"
#include "ActionWithValue.h"
#include "PlumedException.h"

using namespace PLMD;

Value::Value(ActionWithValue&action,const std::string& name):
  action(action),
  value(0.0),
  name(name),
  deriv(false),
  periodicity(unset),
  min(0.0),
  max(0.0),
  max_minus_min(0.0),
  inv_max_minus_min(0.0)
{}

bool Value::isPeriodic()const{
  plumed_massert(periodicity!=unset,"periodicity should be set");
  return periodicity==periodic;
}

void Value::getDomain(double&min,double&max)const{
  min=this->min;
  max=this->max;
}

void Value::setPeriodicity(bool p){
  if(p) periodicity=periodic;
  else periodicity=notperiodic;
}

void Value::setDomain(double min,double max){
  this->min=min;
  this->max=max;
  max_minus_min=max-min;
  if(max_minus_min!=0.0) inv_max_minus_min=1.0/max_minus_min;
}


const std::string Value::getFullName()const{
  if(name.length()==0) return action.getLabel();
  else return action.getLabel()+"."+name;
}

void Value::enableDerivatives()
{
  deriv=true;derivatives.resize(action.getNumberOfParameters());
}



