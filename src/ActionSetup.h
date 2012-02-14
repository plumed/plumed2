#ifndef __PLUMED_ActionSetup_h
#define __PLUMED_ActionSetup_h

#include "Action.h"

namespace PLMD{

//+DEVELDOC MULTI-INHERIT ActionSetup
/**
Used to create a PLMD::Action that do something during setup only e.g. PLMD::SetupUnits
*/
//+ENDDEVELDOC

/// Action executed only at startup
class ActionSetup :
  public virtual Action {
public:
/// Constructor
  ActionSetup(const ActionOptions&ao);
/// Creator of keywords
  static void registerKeywords( Keywords& keys ); 
/// Do nothing.
  void calculate(){};
/// Do nothing.
  void apply(){};
};

}

#endif
