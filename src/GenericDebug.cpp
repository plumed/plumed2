#include "ActionRegister.h"
#include "ActionPilot.h"
#include "ActionSet.h"
#include "PlumedMain.h"

using namespace std;

namespace PLMD{

//+PLUMEDOC GENERIC DEBUG
/**
  Set some debug options

\par syntax
\verbatim
DEBUG [STRIDE=s] [NOVIRIAL] [logActivity]
\endverbatim
The NOVIRIAL flag switches off (for the entire simulation) the contribution
of virial computed in plumed. The logActivity keyword writes in the log
the list of which objects are active and which are inactive
as a sequence of + (active) and - (inactive). Logging is done with stride s.
*/
//+ENDPLUMEDOC
class GenericDebug:
  public ActionPilot
{
  bool logActivity;
  bool logRequestedAtoms;
  bool novirial;
public:
  GenericDebug(const ActionOptions&ao);
/// Register all the relevant keywords for the action  
  static void registerKeywords( Keywords& keys );
  void calculate(){};
  void apply();
};

PLUMED_REGISTER_ACTION(GenericDebug,"DEBUG")

void GenericDebug::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords(keys);
  keys.add("optional","logActivity","write in the log which actions are inactive and which are inactive");
  keys.add("optional","logRequestedAtoms","write in the log which atoms have been requested at a given time");
  keys.add("optional","NOVIRIAL","switch off the virial contribution for the entirity of the simulation");
}

GenericDebug::GenericDebug(const ActionOptions&ao):
Action(ao),
ActionPilot(ao),
logActivity(false),
logRequestedAtoms(false),
novirial(false){
  parseFlag("logActivity",logActivity);
  if(logActivity) log.printf("  logging activity\n");
  parseFlag("logRequestedAtoms",logRequestedAtoms);
  if(logRequestedAtoms) log.printf("  logging requested atoms\n");
  parseFlag("NOVIRIAL",novirial);
  if(novirial) log.printf("  Switching off virial contribution\n");
  if(novirial) plumed.novirial=true;
  checkRead();
}

void GenericDebug::apply(){
  if(logActivity){
    const ActionSet&actionSet(plumed.getActionSet());
    int a=0;
    for(ActionSet::const_iterator p=actionSet.begin();p!=actionSet.end();++p){
      if(dynamic_cast<GenericDebug*>(*p))continue;
      if((*p)->isActive()) a++;
    };
    if(a>0){
      log.printf("activity at step %i: ",getStep());
      for(ActionSet::const_iterator p=actionSet.begin();p!=actionSet.end();++p){
        if(dynamic_cast<GenericDebug*>(*p))continue;
        if((*p)->isActive()) log.printf("+");
        else                 log.printf("-");
      };
      log.printf("\n");
    };
  };
  if(logRequestedAtoms){
    log.printf("requested atoms at step %i: ",getStep());
    int* l;
    int n;
    plumed.cmd("createFullList",&n);
    plumed.cmd("getFullList",&l);
    for(int i=0;i<n;i++) log.printf(" %d",l[i]);
    log.printf("\n");
    plumed.cmd("clearFullList");
  }

}

}

