/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2026 Pengchao Zhang

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "function/Function.h"
#include "core/ActionRegister.h"
#include <algorithm>
#include <cmath>

namespace PLMD {
namespace pathbias {

//+PLUMEDOC FUNCTION PROBABILITY_MIX
/*
Mix two frozen path potentials with a fixed global normalization.

ARG contains two nonperiodic energy scalars V_c,V_A. The output is
-kBT log[(1-eta) exp(-V_c/kBT)+eta exp(-V_A/kBT)/C], where COUPLING is
eta in [0,1) and LOG_NORMALIZER is log C. With C=Z_A/Z_c this samples the
normalized mixture (1-eta)q_c+eta q_A. An approximate frozen C still defines
a conservative potential but changes the actual mixture fraction. C is
global, not the coordinate-dependent normalizer of CONDITIONAL_PATH.
Derivatives are (1-chi,chi). All replicas must use identical parameters and
the same shared full-path inputs. This function neither learns nor biases
directly; apply the result using BIASVALUE.

```plumed
c: DISTANCE ATOMS=1,2 NOPBC
vc: CUSTOM ARG=c FUNC=x*x PERIODIC=NO
va: CUSTOM ARG=c FUNC=2*x*x PERIODIC=NO
total: PROBABILITY_MIX ARG=vc,va KBT=2.5 COUPLING=0.5 LOG_NORMALIZER=0
bias: BIASVALUE ARG=total
```

The example illustrates composition only; its C is not a fitted partition
function ratio. KBT and both arguments must use the same energy units.
At eta=0 this function returns V_c with derivatives (1,0), but PLUMED still
evaluates upstream actions: generate a centroid-only graph to avoid evaluating
an invalid or expensive inactive field. Pure V_A is a separate graph.
*/
//+ENDPLUMEDOC

class ProbabilityMix : public function::Function {
  double kbt=0;
  double coupling=0;
  double logNormalizer=0;
public:
  explicit ProbabilityMix(const ActionOptions&);
  static void registerKeywords(Keywords&);
  void calculate() override;
};

PLUMED_REGISTER_ACTION(ProbabilityMix,"PROBABILITY_MIX")

void ProbabilityMix::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys);
  keys.add("compulsory","KBT","positive thermal energy in the argument energy units");
  keys.add("compulsory","COUPLING","mixture fraction in [0,1)");
  keys.add("compulsory","LOG_NORMALIZER","fixed dimensionless global log partition-function ratio");
  keys.setValueDescription("scalar","the complete-path mixed bias energy");
}

ProbabilityMix::ProbabilityMix(const ActionOptions& ao):Action(ao),Function(ao) {
  parse("KBT",kbt);
  parse("COUPLING",coupling);
  parse("LOG_NORMALIZER",logNormalizer);
  if(getNumberOfArguments()!=2) {
    error("ARG requires the centroid and arithmetic path energies");
  }
  for(unsigned i=0; i<2; ++i) {
    if(getPntrToArgument(i)->getRank()!=0 || getPntrToArgument(i)->isPeriodic()) {
      error("arguments must be nonperiodic scalars");
    }
  }
  if(!std::isfinite(kbt) || kbt<=0 || !std::isfinite(coupling) || coupling<0 || coupling>=1 || !std::isfinite(logNormalizer)) {
    error("KBT must be positive, COUPLING in [0,1), and all parameters finite");
  }
  checkRead();
  addValueWithDerivatives();
  setNotPeriodic();
}

void ProbabilityMix::calculate() {
  const double vc=getArgument(0);
  if(!std::isfinite(vc)) {
    error("centroid energy must be finite");
  }
  double value=vc,chi=0;
  if(coupling>0) {
    const double va=getArgument(1);
    const double lhs=std::log1p(-coupling)-vc/kbt;
    const double rhs=std::log(coupling)-va/kbt-logNormalizer;
    if(!std::isfinite(va) || !std::isfinite(lhs) || !std::isfinite(rhs)) {
      error("mixture component log densities must be finite");
    }
    const double maximum=std::max(lhs,rhs);
    const double sum=std::exp(lhs-maximum)+std::exp(rhs-maximum);
    value=-kbt*(maximum+std::log(sum));
    chi=std::exp(rhs-maximum)/sum;
    if(!std::isfinite(value)) {
      error("mixed energy exceeds floating-point range");
    }
  }
  setValue(value);
  setDerivative(0,1-chi);
  setDerivative(1,chi);
}

}
}
