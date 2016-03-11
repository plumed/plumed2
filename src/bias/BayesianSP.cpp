/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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

//+PLUMEDOC BIAS BAYESIANSP
/*
Calculate a Bayesian Score to use with the Chemical Shifts CV \ref CS2BACKBONE.

The functional form of this bias is
\f[
C=k_BT \sum_{i=1}^{N_{arg}} \log{ \left[ \frac{\pi}{\sqrt{2} \sigma} (x_i+2\sigma^2) \right]} + k_BT\log{\sigma}
\f]

where sigma is an uncertainty parameter,
sampled by a MC algorithm in the bounded interval defined by SIGMA_MIN and SIGMA_MAX.
The initial value is set at SIGMA0. The MC move is a random displacement
of maximum value equal to DSIGMA.



\par Examples
The following input tells plumed to use all the HA chemical shifts with the Bayesian Score and
to print the values of the uncertainty parameter, MC acceptance, and Bayesian score.
\verbatim
WHOLEMOLECULES ENTITY0=1-174
cs:  CS2BACKBONE ATOMS=1-174 DATA=data/ FF=a03_gromacs.mdb FLAT=0.0 NRES=13 ENSEMBLE COMPONENTS
csb: BAYESIANSP ARG=cs.ha SIGMA0=1.0 SIGMA_MIN=0.00001 SIGMA_MAX=10.0 DSIGMA=0.1 NDATA=13
PRINT ARG=csb.sigma,csb.accept,csb.bias
\endverbatim
(See also \ref CS2BACKBONE and \ref PRINT).


*/
//+ENDPLUMEDOC


class BayesianSP : public Bias
{
  const double sqrt2_div_pi;
  // sigma is data uncertainty
  double sigma_;
  double sigma_min_;
  double sigma_max_;
  double Dsigma_;
  // sigma_mean is uncertainty in the mean estimate
  double sigma_mean_;
  // temperature in kbt
  double kbt_;
  // number of data points
  unsigned ndata_;
  // number of replicas
  unsigned nrep_;
  // Monte Carlo stuff
  unsigned MCsteps_;
  unsigned MCstride_;
  unsigned MCaccept_;
  long int MCfirst_;
  // output
  Value* valueBias;
  Value* valueSigma;
  Value* valueAccept;
 
  void doMonteCarlo();
  double getEnergy(const double sigma);
  
public:
  BayesianSP(const ActionOptions&);
  void calculate();
  void update();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(BayesianSP,"BAYESIANSP")

void BayesianSP::registerKeywords(Keywords& keys){
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","SIGMA0","initial value of the uncertainty parameter");
  keys.add("compulsory","SIGMA_MIN","minimum value of the uncertainty parameter");
  keys.add("compulsory","SIGMA_MAX","maximum value of the uncertainty parameter");
  keys.add("compulsory","DSIGMA","maximum MC move of the uncertainty parameter");
  keys.add("compulsory","SIGMA_MEAN","uncertainty in the mean estimate");
  keys.add("optional","NDATA","number of experimental data points");
  keys.add("optional","TEMP","the system temperature - this is only needed if code doesnt' pass the temperature to plumed");
  keys.add("optional","MC_STEPS","number of MC steps");
  keys.add("optional","MC_STRIDE","MC stride");
  componentsAreNotOptional(keys); 
  keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");
  keys.addOutputComponent("sigma", "default","uncertainty parameter");
  keys.addOutputComponent("accept","default","MC acceptance");
}

BayesianSP::BayesianSP(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao), 
sqrt2_div_pi(0.45015815807855),
ndata_(0),
MCsteps_(1), 
MCstride_(1), 
MCaccept_(0), 
MCfirst_(-1)
{
  parse("SIGMA0",sigma_);
  parse("SIGMA_MIN",sigma_min_);
  parse("SIGMA_MAX",sigma_max_);
  parse("DSIGMA",Dsigma_);
  parse("SIGMA_MEAN",sigma_mean_);
  parse("NDATA",ndata_);
  parse("MC_STEPS",MCsteps_);
  parse("MC_STRIDE",MCstride_);
  // get temperature
  double temp=0.0;
  parse("TEMP",temp);

  checkRead();

  if(ndata_==0) ndata_=getNumberOfArguments();
  if(temp>0.0) kbt_=plumed.getAtoms().getKBoltzmann()*temp;
  else kbt_=plumed.getAtoms().getKbT();

  // get number of replicas
  if(comm.Get_rank()==0) nrep_ = multi_sim_comm.Get_size();
  else nrep_ = 0;
  comm.Sum(&nrep_,1);

  // divide sigma_mean by the square root of the number of replicas
  sigma_mean_ /= sqrt(static_cast<double>(nrep_));

  // adjust for multiple-time steps
  MCstride_ *= getStride();

  log.printf("  initial data uncertainty %f\n",sigma_);
  log.printf("  minimum data uncertainty %f\n",sigma_min_);
  log.printf("  maximum data uncertainty %f\n",sigma_max_);
  log.printf("  maximum MC move of data uncertainty %f\n",Dsigma_);
  log.printf("  uncertainty in the mean estimate %f\n",sigma_mean_);
  log.printf("  temperature of the system %f\n",kbt_);
  log.printf("  number of experimental data points %u\n",ndata_);
  log.printf("  number of replicas %u\n",nrep_);
  log.printf("  MC steps %u\n",MCsteps_);
  log.printf("  MC stride %u\n",MCstride_);

  addComponent("bias");   componentIsNotPeriodic("bias");
  addComponent("sigma");  componentIsNotPeriodic("sigma");
  addComponent("accept"); componentIsNotPeriodic("accept");
  valueBias=getPntrToComponent("bias");
  valueSigma=getPntrToComponent("sigma");
  valueAccept=getPntrToComponent("accept");

  // initialize random seed
  unsigned iseed;
  if(comm.Get_rank()==0) iseed = time(NULL)+multi_sim_comm.Get_rank();
  else iseed = 0;     
  comm.Sum(&iseed, 1);
  srand(iseed);
}

double BayesianSP::getEnergy(const double sigma){
  // calculate effective sigma
  const double smean2 = sigma_mean_*sigma_mean_;
  const double s = sqrt( sigma*sigma + smean2 );
  // cycle on arguments
  double ene = 0.0;
  for(unsigned i=0;i<getNumberOfArguments();++i){
    // argument
    const double a2 = 0.5 * getArgument(i) + s*s;
    // increment energy
    ene += std::log( 2.0 * a2 / ( 1.0 - exp(- a2 / smean2) ) );
  }
  // add normalization and Jeffrey's prior
  ene += std::log(s) - static_cast<double>(ndata_)*std::log(sqrt2_div_pi*s);
  return kbt_ * ene;
}

void BayesianSP::doMonteCarlo(){
 // store old energy
 double old_energy = getEnergy(sigma_);
 // cycle on MC steps 
 for(unsigned i=0;i<MCsteps_;++i){
  // propose move
  const double r = static_cast<double>(rand()) / RAND_MAX;
  const double ds = -Dsigma_ + r * 2.0 * Dsigma_;
  double new_sigma = sigma_ + ds;
  // check boundaries
  if(new_sigma > sigma_max_){new_sigma = 2.0 * sigma_max_ - new_sigma;}
  if(new_sigma < sigma_min_){new_sigma = 2.0 * sigma_min_ - new_sigma;}
  // calculate new energy
  double new_energy = getEnergy(new_sigma);
  // accept or reject
  const double delta = ( new_energy - old_energy ) / kbt_;
  // if delta is negative always accept move
  if( delta <= 0.0 ){
   old_energy = new_energy;
   sigma_ = new_sigma;
   MCaccept_++;
  // otherwise extract random number
  } else {
   const double s = static_cast<double>(rand()) / RAND_MAX;
   if( s < exp(-delta) ){
    old_energy = new_energy;
    sigma_ = new_sigma;
    MCaccept_++;
   }
  }
 }
}

void BayesianSP::update(){
  // get time step 
  const long int step = getStep();
  // do MC stuff at the right time step
  if(step%MCstride_==0) doMonteCarlo();
  // this is needed when restarting simulations
  if(MCfirst_==-1) MCfirst_=step;
  // calculate acceptance
  const double MCtrials = floor(static_cast<double>(step-MCfirst_) / static_cast<double>(MCstride_))+1.0;
  const double accept = static_cast<double>(MCaccept_) / static_cast<double>(MCsteps_) / MCtrials;
  // set value of acceptance
  valueAccept->set(accept);
}

void BayesianSP::calculate(){

  // calculate local effective sigma
  const double smean2 = sigma_mean_*sigma_mean_; 
  const double s = sqrt( sigma_*sigma_ + smean2 );
  double ene = 0.0;
  vector<double> f(getNumberOfArguments());
  
  // cycle on arguments (experimental data) 
  if(comm.Get_rank()==0){
   for(unsigned i=0;i<getNumberOfArguments();++i){
     // argument
     const double a2 = 0.5 * getArgument(i) + s*s;
     // useful quantity
     const double t = exp(-a2/smean2);
     const double it = 1./(1.-t);
     // increment energy
     ene += std::log(2.*a2*it);
     // store force
     f[i] = - 0.5/a2 + it*t*0.5/smean2;
   }
   // collect contribution to forces and energy from other replicas
   multi_sim_comm.Sum(&f[0],getNumberOfArguments());
   multi_sim_comm.Sum(&ene,1);
   // add normalizations and priors of local replica
   ene += std::log(s) - static_cast<double>(ndata_)*std::log(sqrt2_div_pi*s);
  } else {
    // forces to zero for non master threads
    for(unsigned i=0;i<getNumberOfArguments();++i) f[i]=0.;
  }
  // intra-replica summation
  comm.Sum(&f[0],getNumberOfArguments());
  comm.Sum(&ene,1);
  
  // set value of the bias
  valueBias->set(kbt_*ene);
  // set value of data uncertainty
  valueSigma->set(sigma_);
  // set forces
  for(unsigned i=0; i<f.size(); ++i) setOutputForce(i, kbt_ * f[i]);
}


}
}


