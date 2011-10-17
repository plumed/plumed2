#include "ActionRegister.h"
#include "PlumedMain.h"

namespace PLMD{

using namespace std;

//+PLUMEDOC GENERIC FLUSH
/**
Periodically flush open files.

\par Example
This input instructs plumed to flush all the output files every 100 steps
\verbatim
FLUSH STRIDE=100
\endverbatim
*/
//+ENDPLUMEDOC

class GenericFlush:
  public Action
{
public:
  GenericFlush(const ActionOptions&ao) : Action(ao) {
    readAction(); checkRead();
  }
/// Do nothing
  void calculate(){};
/// Flush all open files
  void apply(){
    const ActionSet & actionSet(plumed.getActionSet());
    for(ActionSet::const_iterator p=actionSet.begin();p!=actionSet.end();++p) (*p)->fflush();
  }
/// You can't do numerical derivatives for an action flush
  void calculateNumericalDerivatives(){ assert(false); }
};

PLUMED_REGISTER_ACTION(GenericFlush,"FLUSH")

}


