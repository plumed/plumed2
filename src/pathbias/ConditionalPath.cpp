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

//+PLUMEDOC FUNCTION CONDITIONAL_PATH
/*
Evaluate a fixed conditional path log mixture, with its full chain rule.

For positive path score a and positive conditional normalizer m(c), this
function returns log R = log[(1-lambda)+lambda*a/m(c)]. ARG must contain
exactly three nonperiodic scalars: c, log(a), log(m(c)), in that order.
COUPLING is lambda in [0,1). LOWER and UPPER bound the admitted domain of c;
evaluation outside this closed interval fails instead of extrapolating.

The derivatives with respect to the three arguments are (0,chi,-chi),
where chi=lambda*a/(m*R). The explicit c argument only checks the domain.
The normalizer argument MUST retain its dependence on c: freezing a fitted
function does not mean dropping its derivative. All functions must be fixed
during production. This action does not learn m or preserve external file
identity across restarts. Verify frozen input/model hashes before restarting.

A bias V=B_c-kBT*log R is built using CUSTOM and BIASVALUE. kBT must be in
the same energy units as B_c. With an exact m=E_0[a|c], the selected c marginal
of the fixed B_c ensemble is preserved. Approximate m still defines a fixed
potential with full-path weight exp(beta*V), but does not imply that marginal
identity or any efficiency improvement. Every bead in one path shares the
same path weight; beads are not independent samples.

The following two-replica example uses a mean Cartesian distance component,
then squares it, rather than averaging the squared components. It assumes
consistent nonperiodic coordinate images across the replicas. This scalar
example does not implement arbitrary molecular Cartesian-centroid CVs.

```plumed
#SETTINGS NREPLICAS=2
d: DISTANCE ATOMS=1,2 COMPONENTS NOPBC
mean: ENSEMBLE ARG=d.x
c: CUSTOM ARG=mean.d.x FUNC=x*x PERIODIC=NO
h: CUSTOM ARG=d.x FUNC=exp(-x*x/2) PERIODIC=NO
fraction: ENSEMBLE ARG=h
loga: CUSTOM ARG=fraction.h FUNC=log(0.1+x) PERIODIC=NO
logm: CUSTOM ARG=c FUNC=-1+0.2*x PERIODIC=NO
ratio: CONDITIONAL_PATH ARG=c,loga,logm COUPLING=0.5 LOWER=0 UPPER=4
correction: CUSTOM ARG=ratio FUNC=-2.5*x PERIODIC=NO
bias: BIASVALUE ARG=correction
PRINT ARG=c,loga,logm,ratio,bias.* FILE=COLVAR
```

The illustrative logm above is not a fitted conditional reference. Supply a
qualified frozen normalizer and a separate fixed centroid bias for production.
The molecular dynamics adapter must provide consistent dynamical force
scaling and single-path energy accounting. This action does not apply
additional factors of the bead count.
*/
//+ENDPLUMEDOC

class ConditionalPath : public function::Function {
  double coupling;
  double lower;
  double upper;
public:
  explicit ConditionalPath(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  void calculate() override;
};

PLUMED_REGISTER_ACTION(ConditionalPath,"CONDITIONAL_PATH")

void ConditionalPath::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys);
  keys.add("compulsory","COUPLING","mixture fraction in the half-open interval [0,1)");
  keys.add("compulsory","LOWER","inclusive lower limit of the conditioning argument");
  keys.add("compulsory","UPPER","inclusive upper limit of the conditioning argument");
  keys.setValueDescription("scalar","the dimensionless conditional path log mixture");
}

ConditionalPath::ConditionalPath(const ActionOptions& ao):
  Action(ao),
  Function(ao),
  coupling(0),
  lower(0),
  upper(0) {
  parse("COUPLING",coupling);
  parse("LOWER",lower);
  parse("UPPER",upper);
  if(getNumberOfArguments()!=3) {
    error("ARG must contain conditioning coordinate, log score, and log normalizer");
  }
  for(unsigned i=0; i<3; ++i) {
    if(getPntrToArgument(i)->getRank()!=0 || getPntrToArgument(i)->isPeriodic()) {
      error("arguments must be nonperiodic scalars");
    }
  }
  if(!std::isfinite(coupling) || coupling<0 || coupling>=1) {
    error("COUPLING must be finite and in [0,1)");
  }
  if(!std::isfinite(lower) || !std::isfinite(upper) || lower>=upper) {
    error("LOWER and UPPER must be finite and strictly increasing");
  }
  checkRead();
  addValueWithDerivatives();
  setNotPeriodic();
  log.printf("  fixed conditional path mixture with coupling %g in [%g,%g]\n",
             coupling,lower,upper);
}

void ConditionalPath::calculate() {
  const double c=getArgument(0), loga=getArgument(1), logm=getArgument(2);
  if(!std::isfinite(c) || !std::isfinite(loga) || !std::isfinite(logm)) {
    error("conditional path arguments must be finite");
  }
  if(c<lower || c>upper) {
    error("conditioning coordinate is outside the frozen model domain");
  }
  double logratio=0, chi=0;
  if(coupling>0) {
    const double lhs=std::log1p(-coupling);
    const double rhs=std::log(coupling)+loga-logm;
    if(!std::isfinite(rhs)) {
      error("conditional log ratio exceeds floating-point range");
    }
    logratio=std::max(lhs,rhs)+std::log1p(std::exp(-std::abs(lhs-rhs)));
    chi=std::exp(rhs-logratio);
  }
  setValue(logratio);
  setDerivative(0,0);
  setDerivative(1,chi);
  setDerivative(2,-chi);
}

}
}
