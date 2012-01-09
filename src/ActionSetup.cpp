#include "ActionSetup.h"
#include "PlumedMain.h"
#include "ActionSet.h"
#include "PlumedException.h"

using namespace PLMD;

ActionSetup::ActionSetup(const ActionOptions&ao):
  Action(ao)
{
  const ActionSet& actionset(plumed.getActionSet());
  for(ActionSet::const_iterator p=actionset.begin();p!=actionset.end();++p){
// check that all the preceeding actions are ActionSetup
    plumed_massert(dynamic_cast<ActionSetup*>(*p),
      "Action " + getLabel() + " is a setup action, and should be only preceeded by other setup actions");
  }
}

