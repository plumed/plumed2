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
public:
  ColvarEnergy(const ActionOptions&);
/// Clear the output forces
  void clearOutputForces(){};
/// Get the energy from atoms
  void retrieveData();
/// Transfer the energy to the output variables
  void calculate();
/// Apply force to the energy
  void apply();
/// You can't calculate numerical derivatives
  void calculateNumericalDerivatives(){ assert(false); }
};

}

#endif
