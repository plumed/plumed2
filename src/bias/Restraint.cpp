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
#include "core/ActionRegister.h"

namespace PLMD {
namespace bias {

//+PLUMEDOC BIAS RESTRAINT
/*
Adds harmonic and/or linear restraints on one or more variables.

Either or both of SLOPE and KAPPA must be present to specify the linear and harmonic force constants
respectively. If the input arguments are scalars the resulting potential is given by:

$$
  \sum_i \frac{k_i}{2} (x_i-a_i)^2 + m_i*(x_i-a_i)
$$.

The following example input illustrates how you can apply a potential of this type.  This input
tells plumed to restrain the distance between atoms 3 and 5 and the distance between atoms 2 and 4, at different equilibrium
values, and to print the energy of the restraint

```plumed
d1: DISTANCE ATOMS=3,5
d2: DISTANCE ATOMS=2,4
restraint: RESTRAINT ARG=d1,d2 AT=1.0,1.5 KAPPA=150.0,150.0
PRINT ARG=restraint.bias
```

If by contrast the inputs are vectors as in the example below:

```plumed
d1: DISTANCE ATOMS1=3,5 ATOMS2=5,6 ATOMS3=6,7
d2: DISTANCE ATOMS1=2,4 ATOMS2=4,6 ATOMS3=6,8
restraint: RESTRAINT ARG=d1,d2 AT=1.0,1.5 KAPPA=150.0,150.0
PRINT ARG=restraint.bias
```

Then the resulting potential is given by:

$$
\sum_i \sum_j \frac{k_i}{2} (x_{ij}-a_i)^2 + m_i*(x_{ij}-a_i)
$$

The sum over $i$ here is runs over the two arguments, while the sum over $j$ runs over the three components of the input vectors.
Notice, that regardless of whether the input is a scalar, vector or matrix the number of $k_i$, $a_i$ and $m_i$ values  must be
equal to the number of arguments to the action.

*/
//+ENDPLUMEDOC

class Restraint : public Bias {
  std::vector<double> at;
  std::vector<double> kappa;
  std::vector<double> slope;
  Value* valueForce2;
public:
  explicit Restraint(const ActionOptions&);
  void calculate() override;
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(Restraint,"RESTRAINT_SCALAR")

void Restraint::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.setDisplayName("RESTRAINT");
  keys.add("compulsory","SLOPE","0.0","specifies that the restraint is linear and what the values of the force constants on each of the variables are");
  keys.add("compulsory","KAPPA","0.0","specifies that the restraint is harmonic and what the values of the force constants on each of the variables are");
  keys.add("compulsory","AT","the position of the restraint");
  keys.addOutputComponent("force2","default","scalar","the instantaneous value of the squared force due to this bias potential");
}

Restraint::Restraint(const ActionOptions&ao):
  PLUMED_BIAS_INIT(ao),
  at(getNumberOfArguments()),
  kappa(getNumberOfArguments(),0.0),
  slope(getNumberOfArguments(),0.0) {
  parseVector("SLOPE",slope);
  parseVector("KAPPA",kappa);
  parseVector("AT",at);
  checkRead();

  log.printf("  at");
  for(unsigned i=0; i<at.size(); i++) {
    log.printf(" %f",at[i]);
  }
  log.printf("\n");
  log.printf("  with harmonic force constant");
  for(unsigned i=0; i<kappa.size(); i++) {
    log.printf(" %f",kappa[i]);
  }
  log.printf("\n");
  log.printf("  and linear force constant");
  for(unsigned i=0; i<slope.size(); i++) {
    log.printf(" %f",slope[i]);
  }
  log.printf("\n");

  addComponent("force2");
  componentIsNotPeriodic("force2");
  valueForce2=getPntrToComponent("force2");
}


void Restraint::calculate() {
  double ene=0.0;
  double totf2=0.0;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    const double cv=difference(i,at[i],getArgument(i));
    const double k=kappa[i];
    const double m=slope[i];
    const double f=-(k*cv+m);
    ene+=0.5*k*cv*cv+m*cv;
    setOutputForce(i,f);
    totf2+=f*f;
  }
  setBias(ene);
  valueForce2->set(totf2);
}

}


}
