#include "ActionPilot.h"
#include "PlumedMain.h"

using namespace PLMD;
using namespace std;

void ActionPilot::registerKeywords(Keywords& keys){
  keys.add("compulsory","STRIDE","the frequency with which this action is to be performed"); 
}

ActionPilot::ActionPilot(const ActionOptions&ao):
Action(ao),
stride(1)
{
  parse("STRIDE",stride);
  log.printf("  with stride %d\n",stride);
}

bool ActionPilot::onStep()const{
  return getStep()%stride==0;
}

int ActionPilot::getStride()const{
  return stride;
}


