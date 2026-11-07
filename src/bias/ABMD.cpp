/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

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
#include "Bias.h"
#include "tools/Random.h"
#include "tools/Checkpoint.h"
#include "core/ActionRegister.h"
#include <cmath>
#include <ctime>

namespace PLMD {
namespace bias {

//+PLUMEDOC BIAS ABMD
/*
Adds a ratchet-and-pawl like restraint on one or more variables.

This action can be used to evolve a system towards a target value in
CV space using an harmonic potential that moves with the thermal fluctuations of the CV.
This method was introduced in the papers cited below.  In it the biasing potential in this
method is as follows:

$$
V(\rho(t)) = \left \{ \begin{array}{ll} \frac{K}{2}\left(\rho(t)-\rho_m(t)\right)^2, &\rho(t)>\rho_m(t)\\
              0, & \rho(t)\le\rho_m(t), \end{array} \right .
$$

where

$$
\rho(t)=\left(CV(t)-TO\right)^2
$$

and

$$
\rho_m(t)=\min_{0\le\tau\le t}\rho(\tau)+\eta(t)
$$.

The method is based on the introduction of a biasing potential which is zero when
the system is moving towards the desired arrival point and which damps the
fluctuations when the system attempts to move in the opposite direction. As in the
case of the ratchet and pawl system, propelled by thermal motion of the solvent
molecules, the biasing potential does not exert work on the system. $\eta(t)$ is
an additional white noise acting on the minimum position of the bias.

## The LINEAR flag

In the potential above the distance from the target enters through the squared distance $\rho$:
$\rho-\rho_m=(d-d_m)(d+d_m)$, with $d=|CV-TO|$ and $d_m$ the smallest distance reached, so that the effective stiffness
of the ratchet depends on how far the system is from the target (the bias is strong far from the target and
becomes quartic, hence very soft, close to it) and $K$ has the units of an energy divided by the fourth power of the CV.
With the flag `LINEAR` the ratchet is applied to the distance itself:

$$
V(d(t)) = \left \{ \begin{array}{ll} \frac{K}{2}\left(d(t)-d_m(t)\right)^2, &d(t)>d_m(t)\\
              0, & d(t)\le d_m(t), \end{array} \right .
\qquad
d(t)=|CV(t)-TO|,\quad d_m(t)=\min_{0\le\tau\le t}d(\tau)+\eta(t)
$$

In this way the target only fixes the direction in which the pawl lets the system move: the bias is a harmonic
restraint on the fluctuations of the CV away from the closest point reached so far, with the same stiffness at any distance
from the target, and $K$ has the units of an energy divided by the square of the CV (the values used with and without `LINEAR` are not interchangeable).
Note that with `LINEAR` the quantities `MIN`, `NOISE` and the components `_min` are distances (in the units of the CV) and not squared distances.

## Examples

The following input sets up two biases, one on the distance between atoms 3 and 5
and another on the distance between atoms 2 and 4. The two target values are defined
using TO and the two strength using KAPPA. The total energy of the bias is printed.

```plumed
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
ABMD ARG=d1,d2 TO=1.0,1.5 KAPPA=5.0,5.0 LABEL=abmd
PRINT ARG=abmd.bias,abmd.d1_min,abmd.d2_min
```

The same ratchet, applied to the distances and not to the squared distances, is obtained with the flag `LINEAR`:

```plumed
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
ABMD ARG=d1,d2 TO=1.0,1.5 KAPPA=5.0,5.0 LINEAR LABEL=abmd
PRINT ARG=abmd.bias,abmd.d1_min,abmd.d2_min
```

## Restart

With [CHECKPOINT](CHECKPOINT.md) the running minimum, the noise intensity and the random number generators are stored in the checkpoint, so that the restart is exact.
The checkpoint also records whether `LINEAR` was used, since the running minimum is a squared distance without it and a distance with it:
restarting from a checkpoint written with a different setting of `LINEAR` is an error.
With [RESTART](RESTART.md) and no checkpoint the running minimum is not restored: it starts again from `MIN`, or from the value of the argument at the first step.
To continue a run with the same minimum, give the last value of the `_min` components with `MIN`, or use [CHECKPOINT](CHECKPOINT.md).

*/
//+ENDPLUMEDOC

class ABMD : public Bias {
  std::vector<double> to;
  std::vector<double> min;
  std::vector<double> kappa;
  std::vector<double> temp;
  std::vector<int> seed;
  std::vector<Random> random;
  bool linear=false;
public:
  explicit ABMD(const ActionOptions&);
  void calculate() override;
  static void registerKeywords(Keywords& keys);
  bool hasCheckpointState() const override {
    return true;
  }
  void writeCheckpoint(CheckpointWriter& w) override;
  void readCheckpoint(CheckpointReader& r) override;
};

PLUMED_REGISTER_ACTION(ABMD,"ABMD")

void ABMD::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.add("compulsory","TO","The array of target values");
  keys.add("compulsory","KAPPA","The array of force constants.");
  keys.add("optional","MIN","Array of starting values for the bias (set rho_m(t), otherwise it is set using the current value of ARG)");
  keys.add("optional","NOISE","Array of white noise intensities (add a temperature to the ABMD)");
  keys.addFlag("LINEAR",false,"Apply the ratchet to the distance from the target, V=K/2(d-d_m)^2, instead of to the squared distance, V=K/2(d^2-d_m^2)^2: the stiffness does not depend on the distance from the target. MIN, NOISE and the _min components are then distances and not squared distances");
  keys.add("optional","SEED","Array of seeds for the white noise (add a temperature to the ABMD)");
  keys.addOutputComponent("force2","default","scalar","the instantaneous value of the squared force due to this bias potential");
  keys.addOutputComponent("_min","default","scalar","one or multiple instances of this quantity can be referenced elsewhere in the input file. "
                          " These quantities will be named with the arguments of the bias followed by "
                          "the character string _min. These quantities tell the user the minimum value assumed by rho_m(t).");
  keys.addDOI("10.1063/1.478259");
  keys.addDOI("10.1016/j.bpj.2010.01.047");
  keys.addDOI("10.1063/1.3523345");
}

ABMD::ABMD(const ActionOptions&ao):
  PLUMED_BIAS_INIT(ao),
  to(getNumberOfArguments(),0.0),
  min(getNumberOfArguments(),-1.0),
  kappa(getNumberOfArguments(),0.0),
  temp(getNumberOfArguments(),0.0),
  seed(getNumberOfArguments(),std::time(0)),
  random(getNumberOfArguments()) {
  // Note : parseVector will check that number of arguments is correct
  parseVector("KAPPA",kappa);
  parseVector("MIN",min);
  if(min.size()==0) {
    min.assign(getNumberOfArguments(),-1.0);
  }
  if(min.size()!=getNumberOfArguments()) {
    error("MIN array should have the same size as ARG array");
  }
  // optional keywords that are not given leave an empty vector: restore the defaults
  // (the code below read temp[i] and seed[i] out of bounds without NOISE or SEED)
  parseVector("NOISE",temp);
  if(temp.size()==0) {
    temp.assign(getNumberOfArguments(),0.0);
  }
  parseVector("SEED",seed);
  if(seed.size()==0) {
    seed.assign(getNumberOfArguments(),std::time(0));
  }
  parseVector("TO",to);
  parseFlag("LINEAR",linear);
  checkRead();

  if(linear) {
    log.printf("  the ratchet is applied to the distance from the target (LINEAR): V=K/2(d-d_m)^2\n");
  }
  log.printf("  min");
  for(unsigned i=0; i<min.size(); i++) {
    log.printf(" %f",min[i]);
  }
  log.printf("\n");
  log.printf("  to");
  for(unsigned i=0; i<to.size(); i++) {
    log.printf(" %f",to[i]);
  }
  log.printf("\n");
  log.printf("  with force constant");
  for(unsigned i=0; i<kappa.size(); i++) {
    log.printf(" %f",kappa[i]);
  }
  log.printf("\n");

  for(unsigned i=0; i<getNumberOfArguments(); i++) {
    std::string str_min=getPntrToArgument(i)->getName()+"_min";
    addComponent(str_min);
    componentIsNotPeriodic(str_min);
    if(min[i]!=-1.0) {
      getPntrToComponent(str_min)->set(min[i]);
    }
  }
  for(unsigned i=0; i<getNumberOfArguments(); i++) {
    random[i].setSeed(-seed[i]);
  }
  addComponent("force2");
  componentIsNotPeriodic("force2");
}


// The running minimum, the noise intensity (it is set to zero at run time when
// the CV is close enough to the target) and the generators are the state of the
// bias; the components are written by Bias
void ABMD::writeCheckpoint(CheckpointWriter& w) {
  w.write("linear",static_cast<long long>(linear));
  w.write("min",min);
  w.write("temp",temp);
  for(unsigned i=0; i<random.size(); ++i) {
    w.write("random"+std::to_string(i),random[i].stateToString());
  }
  writeComponentValues(w);
}

void ABMD::readCheckpoint(CheckpointReader& r) {
  long long lin=0;
  r.read("linear",lin);
  if((lin!=0)!=linear) {
    error("the checkpoint was written "+std::string(lin!=0?"with":"without")+" LINEAR, but this input "+std::string(linear?"has":"does not have")+
          " it: the running minimum has a different meaning (distance or squared distance)");
  }
  std::vector<double> m,t;
  r.read("min",m);
  r.read("temp",t);
  if(m.size()!=getNumberOfArguments() || t.size()!=getNumberOfArguments()) {
    error("the checkpoint was written with a different number of arguments");
  }
  min=m;
  temp=t;
  for(unsigned i=0; i<random.size(); ++i) {
    std::string state;
    r.read("random"+std::to_string(i),state);
    random[i].stateFromString(state);
  }
  readComponentValues(r);
}

void ABMD::calculate() {
  double ene=0.0;
  double totf2=0.0;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    const double cv=difference(i,to[i],getArgument(i));
    // the quantity on which the ratchet acts: the squared distance from the target or, with LINEAR, the distance
    const double rho=linear ? std::abs(cv) : cv*cv;
    const double k=kappa[i];
    double noise=0.;
    double diff=temp[i];
    if(diff>0) {
      noise = 2.*random[i].Gaussian()*diff;
      if(rho<=diff) {
        diff=0.;
        temp[i]=0.;
      }
    }

    // min < 0 means that the variable has not been used in the input file, so the current position of the CV is used
    // rho < min means that the collective variable is nearer to the target value than at any other previous time so
    // min is set to the CV value
    if(min[i]<0.||rho<min[i]) {
      min[i] = rho;
      setOutputForce(i,0.0);
    } else {
      // otherwise a noise is added to the minimum value
      min[i] += noise;
      // derivative of rho=cv^2 or of d=|cv| with respect to the CV
      const double f = linear ? -k*(rho-min[i])*(cv<0.0?-1.0:1.0) : -2.*k*(rho-min[i])*cv;
      setOutputForce(i,f);
      ene += 0.5*k*(rho-min[i])*(rho-min[i]);
      totf2+=f*f;
    }
    getPntrToComponent(i+1)->set(min[i]);
  }
  setBias(ene);
  getPntrToComponent("force2")->set(totf2);
}

}
}


