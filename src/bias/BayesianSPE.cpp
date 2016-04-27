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

//+PLUMEDOC BIAS BAYESIANSPE
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
csb: BAYESIANSPE ARG=cs.ha SIGMA0=1.0 SIGMA_MIN=0.00001 SIGMA_MAX=10.0 DSIGMA=0.1 NDATA=13
PRINT ARG=csb.sigma,csb.accept,csb.bias
\endverbatim
(See also \ref CS2BACKBONE and \ref PRINT).


*/
//+ENDPLUMEDOC


class BayesianSPE : public Bias
{
  const double sqrt2_div_pi;
  // experimental values
  vector<double> parameters;
  // scale is data scaling factor
  bool   doscale_;
  double scale_;
  double scale_min_;
  double scale_max_;
  double Dscale_;
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
  // Monte Carlo stuff
  unsigned MCsteps_;
  unsigned MCstride_;
  unsigned MCaccept_;
  long int MCfirst_;
  // output
  Value* valueBias;
  Value* valueSigma;
  Value* valueScale;
  Value* valueAccept;
 
  void doMonteCarlo();
  double getEnergy(const double sigma, const double scale);
  
public:
  BayesianSPE(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(BayesianSPE,"BAYESIANSPE")

void BayesianSPE::registerKeywords(Keywords& keys){
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("optional","PARARG","the input for this action is the scalar output from other actions without derivatives."); 
  keys.add("optional","PARAMETERS","the parameters of the arguments in your function");
  keys.addFlag("SCALEDATA",false,"Set to TRUE if you want to sample a scaling factor common to all values and replicas.");  
  keys.add("compulsory","SCALE0","initial value of the uncertainty parameter");
  keys.add("compulsory","SCALE_MIN","minimum value of the uncertainty parameter");
  keys.add("compulsory","SCALE_MAX","maximum value of the uncertainty parameter");
  keys.add("compulsory","DSCALE","maximum MC move of the uncertainty parameter");
  keys.add("compulsory","SIGMA0","initial value of the uncertainty parameter");
  keys.add("compulsory","SIGMA_MIN","minimum value of the uncertainty parameter");
  keys.add("compulsory","SIGMA_MAX","maximum value of the uncertainty parameter");
  keys.add("compulsory","DSIGMA","maximum MC move of the uncertainty parameter");
  keys.add("compulsory","SIGMA_MEAN","uncertainty in the mean estimate");
  keys.add("optional","TEMP","the system temperature - this is only needed if code doesnt' pass the temperature to plumed");
  keys.add("optional","MC_STEPS","number of MC steps");
  keys.add("optional","MC_STRIDE","MC stride");
  componentsAreNotOptional(keys);
  keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");
  keys.addOutputComponent("sigma", "default","uncertainty parameter");
  keys.addOutputComponent("scale", "default","scale parameter");
  keys.addOutputComponent("accept","default","MC acceptance");
}

BayesianSPE::BayesianSPE(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao), 
sqrt2_div_pi(0.45015815807855),
doscale_(false),
ndata_(getNumberOfArguments()),
MCsteps_(1), 
MCstride_(1), 
MCaccept_(0), 
MCfirst_(-1)
{
  parseVector("PARAMETERS",parameters);
  if(parameters.size()!=static_cast<unsigned>(getNumberOfArguments())&&!parameters.empty())
    error("Size of PARAMETERS array should be either 0 or the same as of the number of arguments in ARG1");

  vector<Value*> arg2;
  parseArgumentList("PARARG",arg2);

  if(!arg2.empty()) {
    if(parameters.size()>0) error("It is not possible to use PARARG and PARAMETERS together");
    if(arg2.size()!=getNumberOfArguments()) error("Size of PARARG array should be the same as number for arguments in ARG");
    for(unsigned i=0;i<arg2.size();i++){
      parameters.push_back(arg2[i]->get()); 
      if(arg2[i]->hasDerivatives()==true) error("PARARG can only accept arguments without derivatives");
    }
  }

  if(parameters.size()!=getNumberOfArguments()) 
    error("PARARG or PARAMETERS arrays should include the same number of elements as the arguments in ARG");

  parseFlag("SCALEDATA", doscale_);
  if(doscale_) {
    parse("SCALE0",scale_);
    parse("SCALE_MIN",scale_min_);
    parse("SCALE_MAX",scale_max_);
    parse("DSCALE",Dscale_);
  } else {
    scale_=1.0;
  }

  parse("SIGMA0",sigma_);
  parse("SIGMA_MIN",sigma_min_);
  parse("SIGMA_MAX",sigma_max_);
  parse("DSIGMA",Dsigma_);
  parse("SIGMA_MEAN",sigma_mean_);
  parse("MC_STEPS",MCsteps_);
  parse("MC_STRIDE",MCstride_);
  // get temperature
  double temp=0.0;
  parse("TEMP",temp);

  checkRead();

  if(temp>0.0) kbt_=plumed.getAtoms().getKBoltzmann()*temp;
  else kbt_=plumed.getAtoms().getKbT();

  // get number of replicas
  unsigned nrep_;
  if(comm.Get_rank()==0) nrep_ = multi_sim_comm.Get_size();
  else nrep_ = 0;
  comm.Sum(&nrep_,1);

  // divide sigma_mean by the square root of the number of replicas
  sigma_mean_ /= sqrt(static_cast<double>(nrep_));

  // adjust for multiple-time steps
  MCstride_ *= getStride();

  if(doscale_) {
    log.printf("  sampling a common scaling factor with:\n");
    log.printf("    initial scale parameter %f\n",scale_);
    log.printf("    minimum scale parameter %f\n",scale_min_);
    log.printf("    maximum scale parameter %f\n",scale_max_);
    log.printf("    maximum MC move of scale parameter %f\n",Dscale_);
  }

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

  addComponent("bias");
  componentIsNotPeriodic("bias");
  valueBias=getPntrToComponent("bias");
  if(doscale_) {
    addComponent("scale");
    componentIsNotPeriodic("scale");
    valueScale=getPntrToComponent("scale");
  }
  addComponent("accept");
  componentIsNotPeriodic("accept");
  valueAccept=getPntrToComponent("accept");
  addComponent("sigma");
  componentIsNotPeriodic("sigma");
  valueSigma=getPntrToComponent("sigma");

  // initialize random seed
  unsigned iseed;
  if(comm.Get_rank()==0) iseed = time(NULL)+multi_sim_comm.Get_rank();
  else iseed = 0;     
  comm.Sum(&iseed, 1);
  srand(iseed);
}

double BayesianSPE::getEnergy(const double sigma, const double scale){
  // calculate effective sigma
  const double smean2 = sigma_mean_*sigma_mean_;
  const double s = sqrt( sigma*sigma + smean2 );
  // cycle on arguments
  double ene = 0.0;
  for(unsigned i=0;i<getNumberOfArguments();++i){
    const double dev = scale*getArgument(i)-parameters[i]; 
    // argument
    const double a2 = 0.5*dev*dev + s*s;
    // increment energy
    ene += std::log( 2.0 * a2 / ( 1.0 - exp(- a2 / smean2) ) );
  }
  // add normalization and Jeffrey's prior
  ene += std::log(s) - static_cast<double>(ndata_)*std::log(sqrt2_div_pi*s);
  return kbt_ * ene;
}

void BayesianSPE::doMonteCarlo(){
 // store old energy
 double old_energy = getEnergy(sigma_, scale_);
 // cycle on MC steps 
 for(unsigned i=0;i<MCsteps_;++i){
  // propose move for scale
  double new_scale = scale_;
  if(doscale_) {
    const double r1 = static_cast<double>(rand()) / RAND_MAX;
    const double ds1 = -Dscale_ + r1 * 2.0 * Dscale_;
    new_scale += ds1;
    // check boundaries
    if(new_scale > scale_max_){new_scale = 2.0 * scale_max_ - new_scale;}
    if(new_scale < scale_min_){new_scale = 2.0 * scale_min_ - new_scale;}
    // the scaling factor should be the same for all the replicas
    if(comm.Get_rank()==0) multi_sim_comm.Bcast(new_scale,0);
    comm.Bcast(new_scale,0);
  }

  // propose move for sigma
  const double r2 = static_cast<double>(rand()) / RAND_MAX;
  const double ds2 = -Dsigma_ + r2 * 2.0 * Dsigma_;
  double new_sigma = sigma_ + ds2;
  // check boundaries
  if(new_sigma > sigma_max_){new_sigma = 2.0 * sigma_max_ - new_sigma;}
  if(new_sigma < sigma_min_){new_sigma = 2.0 * sigma_min_ - new_sigma;}

  // calculate new energy
  const double new_energy = getEnergy(new_sigma,new_scale);
  // accept or reject
  const double delta = ( new_energy - old_energy ) / kbt_;
  // if delta is negative always accept move
  if( delta <= 0.0 ){
   old_energy = new_energy;
   scale_ = new_scale;
   sigma_ = new_sigma;
   MCaccept_++;
  // otherwise extract random number
  } else {
   const double s = static_cast<double>(rand()) / RAND_MAX;
   if( s < exp(-delta) ){
    old_energy = new_energy;
    scale_ = new_scale;
    sigma_ = new_sigma;
    MCaccept_++;
   }
  }

  if(doscale_) {
    // the scaling factor should be the same for all the replicas
    if(comm.Get_rank()==0) multi_sim_comm.Bcast(scale_,0);
    comm.Bcast(scale_,0);
  }
 }
}

void BayesianSPE::calculate(){
  // get time step 
  const long int step = getStep();
  // do MC stuff at the right time step
  if(step%MCstride_==0&&!getExchangeStep()) doMonteCarlo();
  // this is needed when restarting simulations
  if(MCfirst_==-1) MCfirst_=step;
  // calculate acceptance
  const double MCtrials = floor(static_cast<double>(step-MCfirst_) / static_cast<double>(MCstride_))+1.0;
  const double accept = static_cast<double>(MCaccept_) / static_cast<double>(MCsteps_) / MCtrials;
  // set value of acceptance
  valueAccept->set(accept);

  // calculate local effective sigma
  const double smean2 = sigma_mean_*sigma_mean_; 
  const double s = sqrt( sigma_*sigma_ + smean2 );
  double ene = 0.0;
  vector<double> f(getNumberOfArguments());
  
  // cycle on arguments (experimental data) 
  if(comm.Get_rank()==0){
   for(unsigned i=0;i<getNumberOfArguments();++i){
     const double dev = scale_*getArgument(i)-parameters[i]; 
     // argument
     const double a2 = 0.5*dev*dev + s*s;
     // useful quantity
     const double t = exp(-a2/smean2);
     const double dt = 1./t;
     const double it = 1./(1.-t);
     const double dit = 1./(1.-dt);
     // increment energy
     ene += std::log(2.*a2*it);
     // store force
     f[i] = -scale_*dev*(dit/smean2 + 1./a2);
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
  // set value of scale parameter 
  if(doscale_) valueScale->set(scale_);
  // set forces
  for(unsigned i=0; i<f.size(); ++i) setOutputForce(i, kbt_ * f[i]);
}


}
}


