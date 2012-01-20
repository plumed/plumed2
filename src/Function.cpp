#include "Function.h"
#include "Colvar.h"

using namespace PLMD;
using namespace std;

Function::Function(const ActionOptions&ao):
Action(ao),
ActionWithValue(ao),
ActionWithArguments(ao)
{
    setNumberOfParameters(getNumberOfArguments());
}

void Function::apply(){

  vector<double>   f(getNumberOfArguments(),0.0);
  bool at_least_one_forced=false;

  for(int i=0;i<getNumberOfValues();++i){
    if(!getValue(i)->checkForced()) continue;
    at_least_one_forced=true;
    const vector<double> & derivatives(getValue(i)->getDerivatives());
    for(unsigned j=0;j<derivatives.size();j++){
      f[j]+=getForce(i)*derivatives[j];
    }
  }
  if(at_least_one_forced) for(unsigned i=0;i<getNumberOfArguments();++i){
    getArguments()[i]->addForce(f[i]);
  }
}

void Function::registerKeywords(Keywords& keys){
  Action::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
}


