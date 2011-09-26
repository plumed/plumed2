#include "ActionWithExternalArguments.h"

using namespace std;
using namespace PLMD;

ActionWithExternalArguments::ActionWithExternalArguments(const ActionOptions& ao ) :
ActionWithValue(ao)
{
}

void ActionWithExternalArguments::readActionWithExternalArguments(const unsigned& nd, const std::vector<double>& d){
  readActionWithValue(nd,d);
}
