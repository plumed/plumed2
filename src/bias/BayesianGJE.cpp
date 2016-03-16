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

//+PLUMEDOC BIAS BAYESIANGJE
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
csb: BAYESIANGJE ARG=cs.ha SIGMA0=1.0 SIGMA_MIN=0.00001 SIGMA_MAX=10.0 DSIGMA=0.1 NDATA=13
PRINT ARG=csb.sigma,csb.accept,csb.bias
\endverbatim
(See also \ref CS2BACKBONE and \ref PRINT).


*/
//+ENDPLUMEDOC


class BayesianGJE : public Bias
{
  const double sqrt2pi;
  // experimental values
  vector<double> parameters;
  // scale is data scaling factor
  double scale_;
  double scale_min_;
  double scale_max_;
  double Dscale_;
  // sigma is data uncertainty
  vector<double> sigma_;
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
  Value* valueScale;
  Value* valueAccept;
  vector<Value*> valueSigma;
  vector<Value*> valueKappa;
 
  void doMonteCarlo();
  double getEnergy(const vector<double> &sigma, const double scale);
  
public:
  BayesianGJE(const ActionOptions&);
  void calculate();
  void update();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(BayesianGJE,"BAYESIANGJE")

void BayesianGJE::registerKeywords(Keywords& keys){
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("optional","PARARG","the input for this action is the scalar output from other actions without derivatives."); 
  keys.add("optional","PARAMETERS","the parameters of the arguments in your function");
  keys.add("compulsory","SCALE0","initial value of the uncertainty parameter");
  keys.add("compulsory","SCALE_MIN","minimum value of the uncertainty parameter");
  keys.add("compulsory","SCALE_MAX","maximum value of the uncertainty parameter");
  keys.add("compulsory","DSCALE","maximum MC move of the uncertainty parameter");
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
  useCustomisableComponents(keys);
  keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");
  keys.addOutputComponent("sigma", "default","uncertainty parameter");
  keys.addOutputComponent("scale", "default","scale parameter");
  keys.addOutputComponent("kappa", "default","intensity of the harmonic restraint");
  keys.addOutputComponent("accept","default","MC acceptance");
}

BayesianGJE::BayesianGJE(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao), 
sqrt2pi(2.506628274631001),
ndata_(0),
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

  parse("SCALE0",scale_);
  parse("SCALE_MIN",scale_min_);
  parse("SCALE_MAX",scale_max_);
  parse("DSCALE",Dscale_);
  vector<double> readsigma;
  parseVector("SIGMA0",readsigma);
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

  if(readsigma.size()==ndata_) {
    sigma_.resize(ndata_);
    sigma_=readsigma;
  } else if(readsigma.size()==1) {
    sigma_.resize(ndata_,readsigma[0]);
  } else {
    plumed_merror("SIGMA0 can accept either one single value or as many values as the number of arguments");
  } 

  // get number of replicas
  unsigned nrep_;
  if(comm.Get_rank()==0) nrep_ = multi_sim_comm.Get_size();
  else nrep_ = 0;
  comm.Sum(&nrep_,1);

  // divide sigma_mean by the square root of the number of replicas
  sigma_mean_ /= sqrt(static_cast<double>(nrep_));

  // adjust for multiple-time steps
  MCstride_ *= getStride();

  log.printf("  initial scale parameter %f\n",scale_);
  log.printf("  minimum scale parameter %f\n",scale_min_);
  log.printf("  maximum scale parameter %f\n",scale_max_);
  log.printf("  maximum MC move of scale parameter %f\n",Dscale_);
  if(readsigma.size()==1) log.printf("  initial data uncertainty %f\n",sigma_[0]);
  else {
    log.printf("  initial data uncertainty");
    for(unsigned i=0;i<sigma_.size();++i) log.printf(" %f", sigma_[i]);
    log.printf("\n");
  }
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
  addComponent("scale");  componentIsNotPeriodic("scale");
  addComponent("accept"); componentIsNotPeriodic("accept");
  valueBias=getPntrToComponent("bias");
  valueScale=getPntrToComponent("scale");
  valueAccept=getPntrToComponent("accept");
  for(unsigned i=0;i<sigma_.size();++i){
    std::string num; Tools::convert(i,num);
    addComponent("sigma_"+num); componentIsNotPeriodic("sigma_"+num);
    addComponent("kappa_"+num); componentIsNotPeriodic("kappa_"+num);
    valueSigma.push_back(getPntrToComponent("sigma_"+num));
    valueKappa.push_back(getPntrToComponent("kappa_"+num));
  }
  // initialize random seed
  unsigned iseed;
  if(comm.Get_rank()==0) iseed = time(NULL)+multi_sim_comm.Get_rank();
  else iseed = 0;     
  comm.Sum(&iseed, 1);
  srand(iseed);
}

double BayesianGJE::getEnergy(const vector<double> &sigma, const double scale){
  // cycle on arguments
  double ene = 0.0;
  const double smean2 = sigma_mean_*sigma_mean_;
  for(unsigned i=0;i<getNumberOfArguments();++i){
    const double ss = sigma[i]*sigma[i] + smean2; 
    const double dev = scale*getArgument(i)-parameters[i]; 
    ene += 0.5*dev*dev/ss + std::log(ss*sqrt2pi);
  }
  return kbt_ * ene;
}

void BayesianGJE::doMonteCarlo(){
 // store old energy
 double old_energy = getEnergy(sigma_, scale_);
 // cycle on MC steps 
 for(unsigned i=0;i<MCsteps_;++i){
   // propose move for scale
   const double r1 = static_cast<double>(rand()) / RAND_MAX;
   const double ds1 = -Dscale_ + r1 * 2.0 * Dscale_;
   double new_scale = scale_ + ds1;
   // check boundaries
   if(new_scale > scale_max_){new_scale = 2.0 * scale_max_ - new_scale;}
   if(new_scale < scale_min_){new_scale = 2.0 * scale_min_ - new_scale;}
   // the scaling factor should be the same for all the replicas
   if(comm.Get_rank()==0) multi_sim_comm.Bcast(new_scale,0);
   comm.Bcast(new_scale,0);
   // propose move for sigma
   vector<double> new_sigma(sigma_.size());
   for(unsigned j=0;j<sigma_.size();++j) {
     const double r2 = static_cast<double>(rand()) / RAND_MAX;
     const double ds2 = -Dsigma_ + r2 * 2.0 * Dsigma_;
     new_sigma[j] = sigma_[j] + ds2;
     // check boundaries
     if(new_sigma[j] > sigma_max_){new_sigma[j] = 2.0 * sigma_max_ - new_sigma[j];}
     if(new_sigma[j] < sigma_min_){new_sigma[j] = 2.0 * sigma_min_ - new_sigma[j];}
   }
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
   // the scaling factor should be the same for all the replicas
   if(comm.Get_rank()==0) multi_sim_comm.Bcast(scale_,0);
   comm.Bcast(scale_,0);
 }
}

void BayesianGJE::update(){
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

void BayesianGJE::calculate(){
  const unsigned ssize = sigma_.size();
  // calculate effective sigma
  vector<double> ss(ssize);
  vector<double> inv_s2(ssize);
  const double smean2 = sigma_mean_*sigma_mean_;

  if(comm.Get_rank()==0) {
    for(unsigned i=0;i<ssize; ++i) {
      ss[i] = sigma_[i]*sigma_[i] + smean2;
      inv_s2[i] = 1.0/ss[i];
    }
  } else {
    for(unsigned i=0;i<ssize; ++i) {
      ss[i] = sigma_[i]*sigma_[i] + smean2;
      inv_s2[i] = 0.;
    }
  }
 
  // inter-replicas summation
  if(comm.Get_rank()==0) multi_sim_comm.Sum(&inv_s2[0],ssize); 
  // intra-replica summation
  comm.Sum(&inv_s2[0],ssize);  
  
  // cycle on arguments
  const unsigned narg=getNumberOfArguments();
  double ene = 0.0;
  for(unsigned i=0;i<narg;++i){
    const double dev = scale_*getArgument(i)-parameters[i]; 
    // increment energy
    ene += 0.5*dev*dev*inv_s2[i] + std::log(ss[i]*sqrt2pi); 
    // set derivatives
    setOutputForce(i, -kbt_*dev*scale_*inv_s2[i]);
    // set value of data uncertainty
    valueSigma[i]->set(sigma_[i]);
    // and of the harmonic restraint
    valueKappa[i]->set(kbt_*inv_s2[i]);
  }

  // set value of the bias
  valueBias->set(kbt_*ene);
  // set value of scale parameter 
  valueScale->set(scale_);
}

}
}


