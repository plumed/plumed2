/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2021 The VES code team
   (see the PEOPLE-VES file at the root of this folder for a list of names)

   See http://www.ves-code.org for more information.

   This file is part of VES code module.

   The VES code module is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The VES code module is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the VES code module.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include "Optimizer.h"
#include "CoeffsVector.h"

#include "core/ActionRegister.h"


namespace PLMD {
namespace ves {

//+PLUMEDOC VES_OPTIMIZER OPT_ROBBINS_MONRO_SGD
/*
Robbins-Monro stochastic gradient decent.

\attention
__This optimizer is only included for reference. We recommend to use the averaged stochastic gradient decent optimizer (\ref OPT_AVERAGED_SGD)__.

\par Examples

*/
//+ENDPLUMEDOC

class Opt_RobbinsMonroSGD : public Optimizer {
private:
  double decay_constant_;
public:
  static void registerKeywords(Keywords&);
  explicit Opt_RobbinsMonroSGD(const ActionOptions&);
  void coeffsUpdate(const unsigned int c_id = 0);
};


PLUMED_REGISTER_ACTION(Opt_RobbinsMonroSGD,"OPT_ROBBINS_MONRO_SGD")


void Opt_RobbinsMonroSGD::registerKeywords(Keywords& keys) {
  Optimizer::registerKeywords(keys);
  Optimizer::useDynamicStepSizeKeywords(keys);
  Optimizer::useMultipleWalkersKeywords(keys);
  Optimizer::useMaskKeywords(keys);
  Optimizer::useRestartKeywords(keys);
  // Optimizer::useMonitorAveragesKeywords(keys);
  Optimizer::useDynamicTargetDistributionKeywords(keys);
  keys.add("optional","DECAY_CONSTANT","the decay constant used for the step size.");
}


Opt_RobbinsMonroSGD::Opt_RobbinsMonroSGD(const ActionOptions&ao):
  PLUMED_VES_OPTIMIZER_INIT(ao),
  decay_constant_(1.0)
{
  parse("DECAY_CONSTANT",decay_constant_);
  if(decay_constant_<1.0) {
    plumed_merror("the value given in DECAY_CONSTANT doesn't make sense, it should be larger than 1.0");
  }
  if(decay_constant_>1.0) {
    log.printf("  using a decay constant of %f\n",decay_constant_);
  }
  checkRead();
}


void Opt_RobbinsMonroSGD::coeffsUpdate(const unsigned int c_id) {
  // getIterationCounterDbl() gives n-1 as it is updated afterwards.
  double current_stepsize =  StepSize(c_id) /(1.0 + getIterationCounterDbl()/decay_constant_);
  setCurrentStepSize(current_stepsize,c_id);
  Coeffs(c_id) += - current_stepsize * CoeffsMask(c_id) * Gradient(c_id);
  //
  double aver_decay = 1.0 / ( getIterationCounterDbl() + 1.0 );
  AuxCoeffs(c_id) += aver_decay * ( Coeffs(c_id)-AuxCoeffs(c_id) );
}


}
}
