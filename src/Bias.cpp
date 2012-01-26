#include "Bias.h"
#include "Colvar.h"

using namespace PLMD;
using namespace std;

Bias::Bias(const ActionOptions&ao):
Action(ao),
ActionPilot(ao),
ActionWithValue(ao),
ActionWithArguments(ao),
outputForces(getArguments().size(),0.0)
{
}

void Bias::registerKeywords( Keywords& keys ){
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
}


void Bias::apply(){
  if(onStep()) for(unsigned i=0;i<getNumberOfArguments();++i){
    getArguments()[i]->addForce(getStride()*outputForces[i]);
  }
}


