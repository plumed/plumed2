#include "Value.h"
#include "ActionWithValue.h"

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
  assert(periodicity!=unset);
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
  return action.getLabel()+"."+name;
}

void Value::enableDerivatives()
{
  deriv=true;derivatives.resize(action.getNumberOfParameters());
}

double Value::difference(double d1,double d2)const{
  assert(periodicity!=unset);
  if(periodicity==periodic){
    double s=(d2-d1)*inv_max_minus_min;
    s=Tools::pbc(s);
    return s*max_minus_min;
  }else{
    return d2-d1;
  }
}





