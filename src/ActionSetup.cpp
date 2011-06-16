#include "ActionSetup.h"
#include "PlumedMain.h"

using namespace PLMD;

ActionSetup::ActionSetup(const ActionOptions&ao):
  Action(ao)
{
  const ActionSet& actionset(plumed.getActionSet());
  for(ActionSet::const_iterator p=actionset.begin();p!=actionset.end();++p){
// check that all the preceeding actions are ActionSetup
    assert(dynamic_cast<ActionSetup*>(*p));
  }
}

