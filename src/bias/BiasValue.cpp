/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
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
#include "ActionRegister.h"


using namespace std;


namespace PLMD {
namespace bias {

//+PLUMEDOC BIAS BIASVALUE
/*
Takes the value of one variable and use it as a bias

This is the simplest possible bias: the bias potential is equal to a collective variable.
It is useful to create custom biasing potential, e.g. applying a function (see \ref Function)
to some collective variable then using the value of this function directly as a bias.

\par Examples

The following input tells plumed to use the value of the distance between atoms 3 and 5
and the value of the distance between atoms 2 and 4 as biases.
It then tells plumed to print the energy of the restraint
\plumedfile
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=3,6 LABEL=d2
BIASVALUE ARG=d1,d2 LABEL=b
PRINT ARG=d1,d2,b.d1_bias,b.d2_bias
\endplumedfile

Another thing one can do is asking one system to follow
a circle in sin/cos according a  time dependence

\plumedfile
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
\endplumedfile

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
  keys.use("ARG");
  // Should be _bias below
  keys.addOutputComponent("_bias","default","one or multiple instances of this quantity can be referenced elsewhere in the input file. "
                          "these quantities will named with  the arguments of the bias followed by "
                          "the character string _bias. These quantities tell the user how much the bias is "
                          "due to each of the colvars.");
}

BiasValue::BiasValue(const ActionOptions&ao):
  PLUMED_BIAS_INIT(ao)
{
  checkRead();
  // add one bias for each argument
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    string ss=getPntrToArgument(i)->getName()+"_bias";
    addComponent(ss); componentIsNotPeriodic(ss);
  }
}

void BiasValue::calculate() {
  double bias=0.0;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    double val; val=getArgument(i);
    getPntrToComponent(i+1)->set(val);
    setOutputForce(i,-1.);
    bias+=val;
  }
  setBias(bias);
}

}
}
