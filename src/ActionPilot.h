#ifndef __PLUMED_ActionPilot_h
#define __PLUMED_ActionPilot_h

#include "Action.h"

namespace PLMD{

/// Action which drives the execution of other Action's.
/// Action's of this kind are executed with a fixed stride
/// which is specified on the directive line with a STRIDE= keyword
class ActionPilot:
  public virtual Action
{
  int stride; // multiple time step
protected:
  int getStride()const;
public:
  ActionPilot(const ActionOptions&);
  bool onStep()const;
};

}

#endif

