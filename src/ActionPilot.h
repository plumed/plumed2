#ifndef __PLUMED_ActionPilot_h
#define __PLUMED_ActionPilot_h

#include "Action.h"

namespace PLMD{

/**
\ingroup MULTIINHERIT
This is used to create PLMD::Action objects that are run with some set frequency.
Any PLMD::Action 
that does not inherit from PLMD::Action is only run when some other Action requires the output from 
it in order to run.  This class is used in PLMD::Bias
 Action which drives the execution of other Action's.
 Action's of this kind are executed with a fixed stride
 which is specified on the directive line with a STRIDE= keyword
*/
class ActionPilot:
  public virtual Action
{
  int stride; // multiple time step
protected:
  int getStride()const;
public:
  ActionPilot(const ActionOptions&);
/// Create the keywords for actionPilot
  static void registerKeywords(Keywords& keys);
/// Check if the action is active on this step
  bool onStep()const;
};

}

#endif

