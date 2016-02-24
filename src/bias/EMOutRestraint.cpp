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

//+PLUMEDOC BIAS EMOUTRESTRAINT
/*
Put the doc here

*/
//+ENDPLUMEDOC


class EMOutrestraint : public Bias
{

  const double sqrt2_div_pi;
  // sigma is data uncertainty
  double sigma_;
  double sigma_min_;
  double sigma_max_;
  double Dsigma_;
  // Monte Carlo stuff
  unsigned MCsteps_;
  unsigned MCstride_;
  unsigned MCaccept_;
  long int MCfirst_;
  // temperature in kbt
  double kbt_;
  // exp data points
  vector<double> ovdd_;
  // output
  Value* valueBias;
  Value* valueSigma;
  Value* valueAccept;
  
  void doMonteCarlo();
  double getEnergy(double sigma);
  
public:
  EMOutrestraint(const ActionOptions&);
  void calculate();
  void update();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(EMOutrestraint,"EMOUTRESTRAINT")

void EMOutrestraint::registerKeywords(Keywords& keys){
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","SIGMA0","initial value of the uncertainty parameter");
  keys.add("compulsory","SIGMA_MIN","minimum value of the uncertainty parameter");
  keys.add("compulsory","SIGMA_MAX","maximum value of the uncertainty parameter");
  keys.add("compulsory","DSIGMA","maximum MC move of the uncertainty parameter");
  keys.add("optional","MC_STEPS","number of MC steps");
  keys.add("optional","MC_STRIDE","MC stride");
  componentsAreNotOptional(keys); 
  keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");
  keys.addOutputComponent("sigma", "default","uncertainty parameter");
  keys.addOutputComponent("accept","default","MC acceptance");
}

EMOutrestraint::EMOutrestraint(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao),
sqrt2_div_pi(0.45015815807855),
MCsteps_(1), 
MCstride_(1), 
MCaccept_(0), 
MCfirst_(-1)
{
  parse("SIGMA0",sigma_);
  parse("SIGMA_MIN",sigma_min_);
  parse("SIGMA_MAX",sigma_max_);
  parse("DSIGMA",Dsigma_);
  parse("MC_STEPS",MCsteps_);
  parse("MC_STRIDE",MCstride_);
  checkRead();

  // last half of the arguments inside ovdd_;
  for(unsigned i=getNumberOfArguments()/2; i<getNumberOfArguments();++i){
    ovdd_.push_back(getArgument(i));
  }

  // get temperature
  kbt_ = plumed.getAtoms().getKbT();
  
  // adjust for multiple-time steps
  MCstride_ *= getStride();

  log.printf("  temperature of the system %f\n",kbt_);
  log.printf("  number of experimental data points %u\n",static_cast<unsigned>(ovdd_.size()));
  log.printf("  initial data uncertainty %f\n",sigma_);
  log.printf("  minimum data uncertainty %f\n",sigma_min_);
  log.printf("  maximum data uncertainty %f\n",sigma_max_);
  log.printf("  maximum MC move of data uncertainty %f\n",Dsigma_);
  log.printf("  MC steps %u\n",MCsteps_);
  log.printf("  MC stride %u\n",MCstride_);

  addComponent("bias");   componentIsNotPeriodic("bias");
  addComponent("sigma");  componentIsNotPeriodic("sigma");
  addComponent("accept"); componentIsNotPeriodic("accept");
  valueBias=getPntrToComponent("bias");
  valueSigma=getPntrToComponent("sigma");
  valueAccept=getPntrToComponent("accept");
  
  // initialize random seed
  unsigned int iseed;
  if(comm.Get_rank()==0) iseed = time(NULL);
  else iseed = 0;     
  comm.Sum(&iseed, 1);
  srand(iseed);
  
}

void EMOutrestraint::doMonteCarlo(){
 // store old energy
 double old_energy = getEnergy(sigma_);
 // cycle on MC steps 
 for(unsigned i=0;i<MCsteps_;++i){
  // propose move
  double r = static_cast<double>(rand()) / RAND_MAX;
  double ds = -Dsigma_ + r * 2.0 * Dsigma_;
  double new_sigma = sigma_ + ds;
  // check boundaries
  if(new_sigma > sigma_max_){new_sigma = 2.0 * sigma_max_ - new_sigma;}
  if(new_sigma < sigma_min_){new_sigma = 2.0 * sigma_min_ - new_sigma;}
  // calculate new energy
  double new_energy = getEnergy(new_sigma);
  // accept or reject
  double delta = ( new_energy - old_energy ) / kbt_;
  // if delta is negative always accept move
  if( delta <= 0.0 ){
   old_energy = new_energy;
   sigma_ = new_sigma;
   MCaccept_++;
  // otherwise extract random number
  } else {
   double s = static_cast<double>(rand()) / RAND_MAX;
   if( s < exp(-delta) ){
    old_energy = new_energy;
    sigma_ = new_sigma;
    MCaccept_++;
   }
  }
 }
}

void EMOutrestraint::update(){
  // get time step 
  long int step = getStep();
  // do MC stuff at the right time step
  if(step%MCstride_==0) doMonteCarlo();
  // this is needed when restarting simulations
  if(MCfirst_==-1) MCfirst_=step;
  // calculate acceptance
  double MCtrials = std::floor(static_cast<double>(step-MCfirst_) / static_cast<double>(MCstride_))+1.0;
  double accept = static_cast<double>(MCaccept_) / static_cast<double>(MCsteps_) / MCtrials;
  // set value of acceptance
  valueAccept->set(accept);
}

double EMOutrestraint::getEnergy(double sigma)
{
  double ene = 0.0;
  // cycle on arguments
  // count number of non-zero overlaps
  double ndata_zero = 0.0;
  for(unsigned i=0;i<getNumberOfArguments()/2;++i){
    // check for zero overlaps
    double ovmd = getArgument(i);
    if(ovmd > 0.0){
     // individual term
     double ene_tmp = std::log(ovmd/ovdd_[i]);
     // increment energy
     ene += std::log(ene_tmp*ene_tmp+ 2.0*sigma*sigma);
     // increment counter
     ndata_zero += 1.0;
    }
  }
  // add normalization and Jeffrey's prior
  ene += std::log(sigma) - ndata_zero*std::log(sqrt2_div_pi*sigma);

  return ene;
}

void EMOutrestraint::calculate()
{   
  double ene = 0.0;
  unsigned int ndata = getNumberOfArguments()/2;
  
  vector<double> ene_der(ndata);
  
  // cycle on arguments
  // count number of non-zero overlaps
  double ndata_zero = 0.0;
  for(unsigned i=0;i<ndata;++i){
    // check for zero overlaps
    double ovmd = getArgument(i);
    if(ovmd > 0.0){
     // individual term
     ene_der[i] = std::log(ovmd/ovdd_[i]);
     // increment energy
     ene += std::log(ene_der[i] * ene_der[i] + 2.0*sigma_*sigma_);
     // increment counter
     ndata_zero += 1.0;
    }
  };
  
  // get derivatives
  for(unsigned i=0;i<ndata;++i){
    // check for zero overlaps
    double ovmd = getArgument(i);
    if(ovmd > 0.0 && ene > 0.0){
     // calculate derivative
     double der =  1.0 / (ene_der[i] * ene_der[i] + 2.0*sigma_*sigma_) * 2.0 * ene_der[i] / ovmd;
     // set forces
     setOutputForce(i, -kbt_*der);
    }
  }
  
  // add normalization and Jeffrey's prior
  ene += std::log(sigma_) - ndata_zero*std::log(sqrt2_div_pi*sigma_);

  // set value of the bias
  valueBias->set(kbt_ * ene);
  // set value of data uncertainty
  valueSigma->set(sigma_);
}


}
}


