#ifndef __PLUMED_ActionSetup_h
#define __PLUMED_ActionSetup_h

#include "Action.h"

namespace PLMD{

/// Action executed only at startup
class ActionSetup : public Action {
protected:
/// Read everything into the action setup
  void readActionSetup();
public:
/// Constructor
  ActionSetup(const ActionOptions&ao);
/// Do nothing.
  void calculate(){};
/// Do nothing.
  void apply(){};
/// You can't do this for an action setup
  void calculateNumericalDerivatives(){ assert(false); }
};

}

#endif
