#include "ActionRegister.h"
#include "ActionPilot.h"
#include "PlumedMain.h"
#include "ActionSet.h"

namespace PLMD{

using namespace std;

//+PLUMEDOC GENERIC FLUSH
/*
This command instructs plumed to flush all the open files with a user specified frequency.  This
is useful for preventing data loss that would otherwise arrise as a consequence of the code
storing data for printing in the buffers

\par Examples
A command like this in the input will instruct plumed to flush all the output files every 100 steps
\verbatim
FLUSH STRIDE=100
\endverbatim
*/
//+ENDPLUMEDOC

class GenericFlush:
  public ActionPilot
{
public:
  GenericFlush(const ActionOptions&ao):
    Action(ao),
    ActionPilot(ao)
  {
    checkRead();
  }
  static void registerKeywords( Keywords& keys );
  void calculate(){};
  void apply(){
    log.flush();
    const ActionSet & actionSet(plumed.getActionSet());
    for(ActionSet::const_iterator p=actionSet.begin();p!=actionSet.end();++p)
    (*p)->fflush();
  }
};

PLUMED_REGISTER_ACTION(GenericFlush,"FLUSH")

void GenericFlush::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  keys.add("compulsory","STRIDE","the frequency with which all the open files should be flushed");
  keys.remove("LABEL");
}

}


