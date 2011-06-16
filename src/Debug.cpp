#include "ActionRegister.h"
#include "ActionPilot.h"
#include "PlumedMain.h"

using namespace std;

namespace PLMD{

//+PLUMEDOC GENERIC DEBUG
/**
  Action to set some debug flags

  This action can be used to enable debug options.
  It accepts a STRIDE option as the \ref Bias actions,
  which is used for debug options which are repeated
  along the simulation.
  It accepts the following flags:
- logActivity: writes (every STRIDE steps) the list
  of which actions are active and which ar inactive in the
  plumed log, in term of a sequence of + (active) and - (inactive)
- NOVIRIAL: switches off (for the entire simulation) the contribution
  of virial computed in plumed
*/
//+ENDPLUMEDOC
class Debug:
  public ActionPilot
{
  bool logActivity;
  bool novirial;
public:
  Debug(const ActionOptions&ao);
  void calculate(){};
  void apply();
};

PLUMED_REGISTER_ACTION(Debug,"DEBUG")

Debug::Debug(const ActionOptions&ao):
Action(ao),
ActionPilot(ao),
logActivity(false),
novirial(false){
  parseFlag("logActivity",logActivity);
  if(logActivity) log.printf("  logging activity\n");
  parseFlag("NOVIRIAL",novirial);
  if(novirial) log.printf("  Switching off virial contribution\n");
  if(novirial) plumed.novirial=true;
  checkRead();
}

void Debug::apply(){
  if(logActivity){
    const ActionSet&actionSet(plumed.getActionSet());
    int a=0;
    for(ActionSet::const_iterator p=actionSet.begin();p!=actionSet.end();++p){
      if(dynamic_cast<Debug*>(*p))continue;
      if((*p)->isActive()) a++;
    };
    if(a>0){
      log.printf("activity at step %i: ",getStep());
      for(ActionSet::const_iterator p=actionSet.begin();p!=actionSet.end();++p){
        if(dynamic_cast<Debug*>(*p))continue;
        if((*p)->isActive()) log.printf("+");
        else                 log.printf("-");
      };
      log.printf("\n");
    };
  };

}

}

