/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include "core/ActionRegister.h"

namespace PLMD {
namespace bias {

//+PLUMEDOC BIAS BIASVALUE
/*
Takes the value of one variable and use it as a bias

This is the simplest possible bias: when you use this method the bias potential is set equal to the scalar-valued argument that was input.
This action is useful for creating custom biasing potential, e.g. if you apply a function (see [function module](module_function.md))
to some collective variable then use the value of this function directly as a bias.

This action is also useful for testing that forces are being calculated correctly.  If I wanted to test whether the derivatives calculated
by the [DISTANCE](DISTANCE.md) action are computed correctly I would write a PLUMED input like this one:

```plumed
d1: DISTANCE ATOMS=1,2
BIASVALUE ARG=d1
```

I would then run a calculation using [driver](driver.md) with the following input:

```plumed
plumed driver --ixyz traj.xyz --debug-forces forces.num
```

This outputs a file with two columns that both contain the forces that are acting upon the atoms.  The first set of forces output are the
analytic forces that are implemented within PLUMED.  The second set of forces are calcluated numerically using finite difference.  If you
method is implemented correctly these two sets of forces should be the same.

## Examples

The following input tells plumed to use the value of the distance between atoms 3 and 5
and the value of the distance between atoms 2 and 4 as biases.
It then tells plumed to print the energy of the restraint

```plumed
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=3,6 LABEL=d2
BIASVALUE ARG=d1,d2 LABEL=b
PRINT ARG=d1,d2,b.d1_bias,b.d2_bias
```

Another thing one can do is ask the system to follow a circle in sin/cos in accordance with a time dependence

```plumed
t: TIME
# this just print cos and sin of time
cos: MATHEVAL ARG=t VAR=t FUNC=cos(t) PERIODIC=NO
sin: MATHEVAL ARG=t VAR=t FUNC=sin(t) PERIODIC=NO
c1: COM ATOMS=1,2
c2: COM ATOMS=3,4
d: DISTANCE COMPONENTS ATOMS=c1,c2
PRINT ARG=t,cos,sin,d.x,d.y,d.z STRIDE=1 FILE=colvar FMT=%8.4f
# this calculates sine and cosine of a projected component of distance
mycos:  MATHEVAL ARG=d.x,d.y  VAR=x,y   FUNC=x/sqrt(x*x+y*y) PERIODIC=NO
mysin:  MATHEVAL ARG=d.x,d.y  VAR=x,y   FUNC=y/sqrt(x*x+y*y) PERIODIC=NO
# this creates a moving spring so that the system follows a circle-like dynamics
# but it is not a bias, it is a simple value now
vv1:  MATHEVAL ARG=mycos,mysin,cos,sin VAR=mc,ms,c,s  FUNC=100*((mc-c)^2+(ms-s)^2) PERIODIC=NO
# this takes the value calculated with matheval and uses as a bias
cc: BIASVALUE ARG=vv1
# some printout
PRINT ARG=t,cos,sin,d.x,d.y,d.z,mycos,mysin,cc.vv1_bias STRIDE=1 FILE=colvar FMT=%8.4f
```

*/
//+ENDPLUMEDOC

class BiasValue : public Bias {
public:
  explicit BiasValue(const ActionOptions&);
  void calculate() override;
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(BiasValue,"BIASVALUE")

void BiasValue::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.addInputKeyword("compulsory","ARG","scalar/vector","the labels of the scalar/vector arguments whose values will be used as a bias on the system");
  // Should be _bias below
  keys.addOutputComponent("_bias","default","scalar","one or multiple instances of this quantity can be referenced elsewhere in the input file. "
                          "these quantities will named with  the arguments of the bias followed by "
                          "the character string _bias. These quantities tell the user how much the bias is "
                          "due to each of the colvars.");
  keys.reset_style("NUMERICAL_DERIVATIVES","hidden");
}

BiasValue::BiasValue(const ActionOptions&ao):
  PLUMED_BIAS_INIT(ao) {
  checkRead();
  // add one bias for each argument
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    std::string ss=getPntrToArgument(i)->getName()+"_bias";
    addComponent(ss);
    componentIsNotPeriodic(ss);
  }
}

void BiasValue::calculate() {
  double bias=0.0;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    double val;
    val=getArgument(i);
    getPntrToComponent(i+1)->set(val);
    setOutputForce(i,-1.);
    bias+=val;
  }
  setBias(bias);
}

}
}
