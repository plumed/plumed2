#include "../src/core/ActionRegister.h"
#include "../src/core/ActionSetup.h"

using namespace std;

namespace PLMD{

/// Action which does nothing
///
/// It is here to test linking of external libraries
class DoNothing:
  public ActionSetup
{
public:
  DoNothing(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(DoNothing,"DONOTHING")

DoNothing::DoNothing(const ActionOptions&ao):
Action(ao),
ActionSetup(ao)
{
}

}

