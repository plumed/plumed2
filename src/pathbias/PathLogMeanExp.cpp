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
#include "tools/Communicator.h"
#include <algorithm>
#include <cmath>
#include <vector>

namespace PLMD {
namespace pathbias {

//+PLUMEDOC FUNCTION PATH_LOGMEANEXP
/*
Compute a stable log arithmetic mean across synchronized replicas.

ARG is one finite, nonperiodic dimensionless scalar ell_b on each replica.
The result is log(sum_b exp(ell_b)/P). Its local derivative is
exp(ell_b)/sum_c exp(ell_c), with no additional factor of 1/P. Only replica
roots communicate across replicas; spatial MPI ranks are not extra beads.
EXPECTED_REPLICAS is required to catch a missing or incorrect replica group.
All replicas must execute the same graph at the same step and time.
Numerical derivatives at the action level are unsupported: finite-difference
tests must perturb one replica's input at a time.

For a common frozen energy field v(s), use ell=-v/kBT and apply -kBT times
this result as a complete-path bias. This averages the effective target to
unbiased probability ratio, not the density estimate or its logarithm.

```plumed
#SETTINGS NREPLICAS=2
d: DISTANCE ATOMS=1,2 NOPBC
v: CUSTOM ARG=d FUNC=x*x PERIODIC=NO
ell: CUSTOM ARG=v FUNC=-x/2.5 PERIODIC=NO
lm: PATH_LOGMEANEXP ARG=ell EXPECTED_REPLICAS=2
energy: CUSTOM ARG=lm FUNC=-2.5*x PERIODIC=NO
bias: BIASVALUE ARG=energy
```

Verify identical field parameters/state hashes before a run or restart.
This action has no bias history and does not make adaptive OPES learning
compatible with this Hamiltonian. If v is an active frozen bias action,
apply the correction V_path-v instead of adding V_path on top of it.
The molecular dynamics adapter owns dynamical force scaling and single-path
energy accounting. Uniform bead observables still use one path weight
exp(V_path/kBT); the local derivative weights are not observable weights.
*/
//+ENDPLUMEDOC

class PathLogMeanExp : public function::Function {
  int replicas=1;
  int replica=0;
  void collectiveCheck(int invalid);
public:
  explicit PathLogMeanExp(const ActionOptions&);
  static void registerKeywords(Keywords&);
  void calculate() override;
};

PLUMED_REGISTER_ACTION(PathLogMeanExp,"PATH_LOGMEANEXP")

void PathLogMeanExp::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys);
  keys.add("compulsory","EXPECTED_REPLICAS","number of synchronized replicas in this path");
  keys.setValueDescription("scalar","the dimensionless log arithmetic mean of replica ratios");
}

void PathLogMeanExp::collectiveCheck(int invalid) {
  comm.Sum(invalid);
  if(comm.Get_rank()==0) {
    multi_sim_comm.Sum(invalid);
  }
  comm.Bcast(invalid,0);
  if(invalid) {
    error("invalid scalar, replica count, numerical derivative request, or nonfinite input on at least one replica");
  }
}

PathLogMeanExp::PathLogMeanExp(const ActionOptions& ao):Action(ao),Function(ao) {
  if(comm.Get_rank()==0) {
    replicas=multi_sim_comm.Get_size();
    replica=multi_sim_comm.Get_rank();
  }
  comm.Bcast(replicas,0);
  comm.Bcast(replica,0);
  int expected=0;
  parse("EXPECTED_REPLICAS",expected);
  int invalid=(expected<1 || expected!=replicas || getNumberOfArguments()!=1);
  if(getNumberOfArguments()==1) {
    invalid=invalid || getPntrToArgument(0)->getRank()!=0 || getPntrToArgument(0)->isPeriodic();
  }
  invalid=invalid || checkNumericalDerivatives();
  collectiveCheck(invalid);
  checkRead();
  addValueWithDerivatives();
  setNotPeriodic();
  log.printf("  stable log mean over %d synchronized replicas\n",replicas);
}

void PathLogMeanExp::calculate() {
  const double value=getArgument(0);
  collectiveCheck(!std::isfinite(value));
  std::vector<double> values(3*replicas,0.0);
  if(comm.Get_rank()==0) {
    values[replica]=value;
    values[replicas+replica]=static_cast<double>(getStep());
    values[2*replicas+replica]=getTime();
    multi_sim_comm.Sum(values);
  }
  comm.Bcast(values,0);
  for(int b=0; b<replicas; ++b) {
    if(values[replicas+b]!=values[replicas] || values[2*replicas+b]!=values[2*replicas]) {
      error("replicas must evaluate at identical steps and times");
    }
  }
  const double maximum=*std::max_element(values.begin(),values.begin()+replicas);
  double sum=0;
  for(int b=0; b<replicas; ++b) {
    sum+=std::exp(values[b]-maximum);
  }
  setValue(maximum+std::log(sum/replicas));
  setDerivative(0,std::exp(value-maximum)/sum);
}

}
}
