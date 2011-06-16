#include "ActionPilot.h"
#include "PlumedMain.h"
#include <cassert>

using namespace PLMD;
using namespace std;

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


