/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The PLUMED team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of PLUMED, version 2.0.

   PLUMED is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   PLUMED is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with PLUMED.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Bias.h"
#include "ActionRegister.h"

#include <cassert>

using namespace std;


namespace PLMD{

//+PLUMEDOC BIAS UPPER_WALLS
/*
Add an upper wall on one or more variables.

The expression for the bias due to the wall is given by:

\f[
  \sum_i {k_i}((x_i-a_i+o_i)/s_i)^e_i
\f]

\par Examples
The following input tells plumed to restrain the distance between atoms 3 and 5
and the distance between atoms 2 and 4, at different equilibrium
values. It then tells plumed to print the energy of the restraint
\verbatim
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
UPPER_WALLS ARG=d1,d2 AT=1.0,1.5 KAPPA=150.0,150.0 EXP=2,2 EPS=1,1 OFFSET 0,0 LABEL=uwall
LOWER_WALLS ARG=d1,d2 AT=0.0,1.0 KAPPA=150.0,150.0 EXP=2,2 EPS=1,1 OFFSET 0,0 LABEL=lwall
PRINT ARG=uwall.bias,lwall.bias
\endverbatim
(See also \ref DISTANCE and \ref PRINT).

*/
//+ENDPLUMEDOC

class BiasUWalls : public Bias{
  std::vector<double> at;
  std::vector<double> kappa;
  std::vector<double> exp;
  std::vector<double> eps;
  std::vector<double> offset;
public:
  BiasUWalls(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(BiasUWalls,"UPPER_WALLS")

void BiasUWalls::registerKeywords(Keywords& keys){
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","AT","the positions of the wall. The a_i in the expression for a wall.");
  keys.add("compulsory","KAPPA","the force constant for the wall.  The k_i in the expression for a wall.");
  keys.add("compulsory","OFFSET","0.0","the offset for the start of the wall.  The o_i in the expression for a wall.");
  keys.add("compulsory","EXP","2.0","the powers for the walls.  The e_i in the expression for a wall.");
  keys.add("compulsory","EPS","1.0","the values for s_i in the expression for a wall");
}

BiasUWalls::BiasUWalls(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao),
at(getNumberOfArguments(),0),
kappa(getNumberOfArguments(),0.0),
exp(getNumberOfArguments(),2.0),
eps(getNumberOfArguments(),1.0),
offset(getNumberOfArguments(),0.0)
{
  parseVector("OFFSET",kappa);
  plumed_assert(offset.size()==getNumberOfArguments());
  parseVector("EPS",kappa);
  plumed_assert(eps.size()==getNumberOfArguments());
  parseVector("EXP",kappa);
  plumed_assert(exp.size()==getNumberOfArguments());
  parseVector("KAPPA",kappa);
  plumed_assert(kappa.size()==getNumberOfArguments());
  parseVector("AT",at);
  plumed_assert(at.size()==getNumberOfArguments());
  checkRead();

  log.printf("  at");
  for(unsigned i=0;i<at.size();i++) log.printf(" %f",at[i]);
  log.printf("\n");
  log.printf("  with an offset");
  for(unsigned i=0;i<offset.size();i++) log.printf(" %f",offset[i]);
  log.printf("\n");
  log.printf("  with force constant");
  for(unsigned i=0;i<kappa.size();i++) log.printf(" %f",kappa[i]);
  log.printf("\n");
  log.printf("  and exponent");
  for(unsigned i=0;i<exp.size();i++) log.printf(" %f",exp[i]);
  log.printf("\n");
  log.printf("  rescaled");
  for(unsigned i=0;i<eps.size();i++) log.printf(" %f",eps[i]);
  log.printf("\n");

  addComponent("bias");
  addComponent("force2");
}

void BiasUWalls::calculate(){
  double ene=0.0;
  double totf2=0.0;
  for(unsigned i=0;i<getNumberOfArguments();++i){
    const double cv=difference(i,at[i],getArgument(i));
    const double k=kappa[i];
    const double exponent=exp[i];
    const double epsilon=eps[i];
    const double off=offset[i];
    const double uscale = (cv+off)/epsilon;
    if(uscale>0.) {
      const double f=-(k/epsilon)*exponent*pow(uscale, exponent-1);
      ene+=k*pow(uscale, exponent);
      setOutputForce(i,f);
      totf2+=f*f;
    }
  };
  getPntrToComponent("bias")->set(ene);
  getPntrToComponent("force2")->set(totf2);
}

}
