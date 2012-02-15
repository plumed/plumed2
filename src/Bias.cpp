#include "Bias.h"
#include "Colvar.h"

using namespace PLMD;
using namespace std;

Bias::Bias(const ActionOptions&ao):
Action(ao),
ActionPilot(ao),
ActionWithValue(ao),
ActionWithArguments(ao),
outputForces(getNumberOfArguments(),0.0)
{
}

void Bias::registerKeywords( Keywords& keys ){
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.add("compulsory","STRIDE","1","the frequency with which the forces due to the bias should be calculated.  This can be used to correctly set up multistep algorithms");
}


void Bias::apply(){
  if(onStep()) for(unsigned i=0;i<getNumberOfArguments();++i){
    getPntrToArgument(i)->addForce(getStride()*outputForces[i]);
  }
}


