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
  // temperature in kbt
  double kbt_;
  // exp data points
  vector<double> ovdd_;
  // output
  Value* valueBias;
  // parallel stuff
  bool serial_;
  unsigned size_;
  unsigned rank_;
  
public:
  EMrestraint(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(EMrestraint,"EMRESTRAINT")

void EMrestraint::registerKeywords(Keywords& keys){
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.addFlag("SERIAL",false,"perform the calculation in serial - for debug purpose");
  componentsAreNotOptional(keys); 
  keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");
}

EMrestraint::EMrestraint(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao),
serial_(false)
{

  // serial or parallel
  parseFlag("SERIAL",serial_);
  if(serial_){
    size_=1; rank_=0;
  } else {
    size_=comm.Get_size(); rank_=comm.Get_rank();
  }

  checkRead();

  // last half of the arguments inside ovdd_;
  for(unsigned i=getNumberOfArguments()/2; i<getNumberOfArguments();++i){
    ovdd_.push_back(getArgument(i));
  }

  // get temperature
  kbt_ = plumed.getAtoms().getKbT();

  log.printf("  temperature of the system %f\n",kbt_);
  log.printf("  number of experimental data points %u\n",static_cast<unsigned>(ovdd_.size()));
  if(serial_) log.printf("  serial calculation\n");

  addComponent("bias");   componentIsNotPeriodic("bias");

  valueBias=getPntrToComponent("bias");
  
}


void EMrestraint::calculate(){
   
  double ene = 0.0;
  vector<double> ene_der(getNumberOfArguments()/2);
  
  // cycle on arguments 
  for(unsigned i=rank_;i<getNumberOfArguments()/2;i=i+size_){
    // check for zero overlaps
    double ovmd = getArgument(i);
    if(ovmd > 0.0){
     // individual term
     ene_der[i] = std::log(ovmd/ovdd_[i]);
     // increment energy
     ene += ene_der[i] * ene_der[i];
    }
  };
  
  if(!serial_){
   comm.Sum(&ene, 1);
   comm.Sum(&ene_der[0], ene_der.size());
  }

  // constant factor
  double fact = kbt_ * 0.5 * static_cast<double>(ovdd_.size());

  // get derivatives
  for(unsigned i=0;i<getNumberOfArguments()/2;++i){
    // check for zero overlaps
    double ovmd = getArgument(i);
    if(ovmd > 0.0 && ene > 0.0){
     // calculate derivative
     double der = 2.0 * fact / ene * ene_der[i] / ovmd;
     // set forces
     setOutputForce(i, -der);
    }
  }

  // set value of the bias
  valueBias->set(fact * std::log(ene));
}


}
}


