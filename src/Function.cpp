#include "Function.h"
#include "Colvar.h"
#include <cassert>

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

  for(int i=0;i<getNumberOfValues();++i){
    const vector<double> & derivatives(getValue(i)->getDerivatives());
    for(unsigned j=0;j<derivatives.size();j++){
      f[j]+=getForce(i)*derivatives[j];
    }
  }
  for(unsigned i=0;i<getNumberOfArguments();++i){
    getArguments()[i]->addForce(f[i]);
  }
}



