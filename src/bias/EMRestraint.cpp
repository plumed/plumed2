/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016 The plumed team
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
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include <cmath>

using namespace std;

namespace PLMD{
namespace bias{

//+PLUMEDOC BIAS EMRESTRAINT
/*
Put the doc here

*/
//+ENDPLUMEDOC


class EMrestraint : public Bias
{
  const double sqrt2pi;
  // temperature in kbt
  double kbt_;
  // exp data points
  vector<double> ovdd_;
  // output
  Value* valueBias;
  
public:
  EMrestraint(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(EMrestraint,"EMRESTRAINT")

void EMrestraint::registerKeywords(Keywords& keys){
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","EXPVALUES","experimental data points");
  componentsAreNotOptional(keys); 
  keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");
}

EMrestraint::EMrestraint(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao), 
sqrt2pi(2.506628274631001)
{

  parseVector("EXPVALUES", ovdd_);
  checkRead();

  // check if experimental data points are as many as arguments
  if(ovdd_.size()!=getNumberOfArguments()) error("Wrong number of experimental data points\n");

  // get temperature
  kbt_ = plumed.getAtoms().getKbT();

  log.printf("  temperature of the system %f\n",kbt_);
  log.printf("  number of experimental data points %u\n",ovdd_.size());

  addComponent("bias");   componentIsNotPeriodic("bias");

  valueBias=getPntrToComponent("bias");

}


void EMrestraint::calculate(){

  
  // cycle on arguments 
  double ene = 0.0;
  for(unsigned i=0;i<getNumberOfArguments();++i){
    // increment energy
    ene += getArgument(i); 
    // set derivatives
    setOutputForce(i, -0.5*kbt_);
  };

  // add normalizations and priors
  ene = 0.5*ene;
  // set value of the bias
  valueBias->set(kbt_*ene);
}


}
}


