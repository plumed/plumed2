#include "Value.h"
#include "ActionWithValue.h"

using namespace PLMD;

const std::string Value::getFullName()const{
  return action.getLabel()+"."+name;
}

void Value::enableDerivatives()
{
  deriv=true;derivatives.resize(action.getNumberOfParameters());
}


