#ifndef __PLUMED_ActionWithExternalArguments_h
#define __PLUMED_ActionWithExternalArguments_h

#include "ActionWithValue.h"
#include "PlumedMain.h"
#include <vector>
#include <set>

namespace PLMD {

/// Action which can takes in data from the MD code
class ActionWithExternalArguments : public ActionWithValue {
protected:
/// Read in keywords related to external arguments
  void readActionWithExternalArguments( const unsigned& nd, const std::vector<double>& d );
public:
  ActionWithExternalArguments(const ActionOptions&ao);
/// Clear the data in the output force arrays
//  virtual void clearOutputForces()=0;
/// Retrieve the data from the MD code that is used to calculate this action
  virtual void retrieveData()=0;
};

}

#endif

