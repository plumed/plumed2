#ifndef __PLUMED_ActionWithExternalArguments_h
#define __PLUMED_ActionWithExternalArguments_h

#include "ActionWithValue.h"
#include "PlumedMain.h"
#include <vector>
#include <set>

namespace PLMD {

/// Action which can takes in data from the MD code
class ActionWithExternalArguments : public ActionWithValue {
public:
  ActionWithExternalArguments(const ActionOptions&ao);
  virtual void clearOutputForces()=0;
  virtual void retrieveData()=0;
};

}

#endif

