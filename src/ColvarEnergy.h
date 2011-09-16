#ifndef __PLUMED_ColvarEnergy_h
#define __PLUMED_ColvarEnergy_h

#include "ActionWithExternalArguments.h"

#include <string>
#include <cmath>
#include <cassert>

namespace PLMD{

class ColvarEnergy : public ActionWithExternalArguments {
private:
  double energy;
  double forceOnEnergy;
public:
  ColvarEnergy(const ActionOptions&);
// active methods:
  virtual void clearOutputForces();
  virtual void retrieveData();
  virtual void calculate();
  virtual void apply();
};

}

#endif
