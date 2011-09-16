#include "ActionRegister.h"
#include "PlumedMain.h"

namespace PLMD{

using namespace std;

//+PLUMEDOC GENERIC FLUSH
/**
Periodically flush open files.

\par Syntax
\verbatim
FLUSH [STRIDE=s]
\endverbatim
This directive is used to flush all the open files periodically.
It understands the keyword STRIDE,
which is the number of timesteps between flushings

\par Example
This input is flushing all the output files every 100 steps
\verbatim
FLUSH STRIDE=100
\endverbatim
*/
//+ENDPLUMEDOC

class GenericFlush:
  public Action
{
public:
  GenericFlush(const ActionOptions&ao):
    Action(ao)
  {
    strideKeywordIsCompulsory();
    readAction();
    checkRead();
  }
  void calculate(){};
  void apply(){
    const ActionSet & actionSet(plumed.getActionSet());
    for(ActionSet::const_iterator p=actionSet.begin();p!=actionSet.end();++p)
    (*p)->fflush();
  }
};

PLUMED_REGISTER_ACTION(GenericFlush,"FLUSH")

}


