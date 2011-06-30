#include "ActionRegister.h"
#include "ActionPilot.h"
#include "PlumedMain.h"

namespace PLMD{

using namespace std;

//+PLUMEDOC GENERIC FLUSH
/**
Periodically flush open files.

This Action is used to flush the open files periodically.
Similarly the other Actions, it understands the keyword STRIDE,
which is the number of timesteps between flushings
\verbatim
# This is flushing all output files every 100 steps
FLUSH STRIDE=100
\endverbatim
*/
//+ENDPLUMEDOC

class Flush:
  public ActionPilot
{
public:
  Flush(const ActionOptions&ao):
    Action(ao),
    ActionPilot(ao)
  {
    checkRead();
  }
  void calculate(){};
  void apply(){
    const ActionSet & actionSet(plumed.getActionSet());
    for(ActionSet::const_iterator p=actionSet.begin();p!=actionSet.end();++p)
    (*p)->fflush();
  }
};

PLUMED_REGISTER_ACTION(Flush,"FLUSH")

}


