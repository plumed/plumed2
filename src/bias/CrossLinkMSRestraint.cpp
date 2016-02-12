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

#include <math.h>

using namespace std;

namespace PLMD{
namespace bias{

//+PLUMEDOC BIAS CROSSLINKMSRESTRAINT
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
csb: BAYESIANSP ARG=cs.ha SIGMA0=1.0 SIGMA_MIN=0.00001 SIGMA_MAX=10.0 DSIGMA=0.1 KBT=2.494 MC_STEPS=1 MC_STRIDE=1 MC_SEED=1234
PRINT ARG=csb.sigma,csb.accept,csb.bias
\endverbatim
(See also \ref CS2BACKBONE and \ref PRINT).


*/
//+ENDPLUMEDOC


class CrossLinkMSRestraint : public Bias
{
// constants
  double sq2;
  double sqPi;
  double sq2Pi;
// uncertainty parameter
  double sigma_;
  double sigma_min_;
  double sigma_max_;
  double Dsigma_;
// psi parameter
  std::vector<double> psi_;
  double psi_min_;
  double psi_max_;
  double Dpsi_;
// XL length 
  double length_;
// slope
  double slope_;
  bool entropic_;
// data classes
  std::vector<int> ndata_;
// temperature in kbt
  double kbt_;
// Monte Carlo stuff
  int MCsteps_;
  int MCstride_;
  unsigned int MCaccsig_;
  std::vector<unsigned int> MCaccpsi_;
  long int MCfirst_;
 
  void doMonteCarlo();
  double getEnergy(double sigma, double psi, int i0, int i1);
  double getEnergy(double sigma, std::vector<double> psi);
  double proposeMove(double x, double xmin, double xmax, double dxmax);
  bool doAccept(double oldE, double newE);

public:
  CrossLinkMSRestraint(const ActionOptions&);
  void calculate();
  void update();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(CrossLinkMSRestraint,"CROSSLINKMSRESTRAINT")

void CrossLinkMSRestraint::registerKeywords(Keywords& keys){
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","SIGMA0","initial value of the uncertainty parameter");
  keys.add("compulsory","SIGMA_MIN","minimum value of the uncertainty parameter");
  keys.add("compulsory","SIGMA_MAX","maximum value of the uncertainty parameter");
  keys.add("compulsory","DSIGMA","maximum MC move of the uncertainty parameter");
  keys.add("compulsory","PSI0","initial value of the psi parameter");
  keys.add("compulsory","PSI_MIN","minimum value of the psi parameter");
  keys.add("compulsory","PSI_MAX","maximum value of the psi parameter");
  keys.add("compulsory","DPSI","maximum MC move of the psi parameter");
  keys.add("compulsory","LENGTH","length of XL");
  keys.add("compulsory","KBT","temperature");
  keys.add("compulsory","NDATA","number of data points per class");
  keys.add("optional","SLOPE","slope");
  keys.addFlag("ENTROPY",false,"add entropic contribution");
  keys.add("optional","MC_STEPS","number of MC steps");
  keys.add("optional","MC_STRIDE","MC stride");
  componentsAreNotOptional(keys); 
  keys.addOutputComponent("bias",   "default","the instantaneous value of the bias potential");
  keys.addOutputComponent("sigma",  "default","uncertainty parameter");
  keys.addOutputComponent("accsig", "default","MC acceptance sigma");
  keys.addOutputComponent("psi",    "COMPONENTS","psi parameter");
  keys.addOutputComponent("accpsi", "COMPONENTS","MC acceptance psi");
}

CrossLinkMSRestraint::CrossLinkMSRestraint(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao),
MCsteps_(1), MCstride_(1), MCaccsig_(0), MCfirst_(-1), slope_(0.0),
entropic_(false)
{
  double psi0;
  parse("SIGMA0",   sigma_);
  parse("SIGMA_MIN",sigma_min_);
  parse("SIGMA_MAX",sigma_max_);
  parse("DSIGMA",   Dsigma_);
  parse("PSI0",     psi0);
  parse("PSI_MIN",  psi_min_);
  parse("PSI_MAX",  psi_max_);
  parse("DPSI",     Dpsi_);
  parse("LENGTH",   length_);
  parse("SLOPE",    slope_);
  parse("KBT",      kbt_);
  parseVector("NDATA", ndata_);
  parse("MC_STEPS", MCsteps_);
  parse("MC_STRIDE",MCstride_);

  parseFlag("ENTROPY",entropic_);

  checkRead();

  // adjust for multiple-time steps
  MCstride_ *= getStride();

  log.printf("  initial value of uncertainty %f\n",sigma_);
  log.printf("  minimum value of uncertainty %f\n",sigma_min_);
  log.printf("  maximum value of uncertainty %f\n",sigma_max_);
  log.printf("  maximum MC move of the uncertainty parameter %f\n",Dsigma_);
  log.printf("  initial value of psi %f\n",psi0);
  log.printf("  minimum value of psi %f\n",psi_min_);
  log.printf("  maximum value of psi %f\n",psi_max_);
  log.printf("  maximum MC move of the psi parameter %f\n",Dpsi_);
  log.printf("  XL length %f\n", length_);
  log.printf("  slope %f\n", slope_);
  log.printf("  temperature of the system %f\n",kbt_);
  log.printf("  number of data points %d\n",getNumberOfArguments());
  log.printf("  number of MC steps %d\n",MCsteps_);
  log.printf("  do MC every %d steps\n", MCstride_);

  addComponent("bias");   componentIsNotPeriodic("bias");
  addComponent("sigma");  componentIsNotPeriodic("sigma");
  addComponent("accsig"); componentIsNotPeriodic("accsig");
  for(unsigned i=0;i<ndata_.size();i++) {
      std::string num; Tools::convert(i,num);
      addComponent("psi"+num);    componentIsNotPeriodic("psi"+num);
      addComponent("accpsi"+num); componentIsNotPeriodic("accpsi"+num);
      psi_.push_back(psi0);
      MCaccpsi_.push_back(0);      
  }

  // initialize random seed
  srand (time(NULL));

  // initialize constants
  sq2   = 1.4142135623730951;
  sqPi  = 1.7724538509055159;
  sq2Pi = 2.5066282746310002;

}

// used to update a single psi
double CrossLinkMSRestraint::getEnergy(double sigma, double psi, int i0, int i1)
{
  // calculate stuff
  double sig  = sq2 * sigma;
  double sig2 = sig * sig;
  // energy
  double ene = 0.0;
  // cycle on arguments 
  for(unsigned i=i0;i<i1;++i){
    // get distance
    double dist = getArgument(i);
    // get various useful stuff
    double log_eLpR_2  = -(dist+length_)*(dist+length_)/(2.0*sig2);
    double eLpR_2      = std::exp(log_eLpR_2);
    double log_e2LR    = 2.0 * length_ * dist / sig2;
    double eLpR_2_e2LR = std::exp(log_eLpR_2+log_e2LR);
    double erfLmR = erf((length_-dist)/sq2/sig);
    double erfLpR = erf((length_+dist)/sq2/sig);
    // calculate posterior
    double post = -sig/(sq2Pi*dist) * (eLpR_2_e2LR - eLpR_2) + 0.5 * (erfLmR + erfLpR);
    // and add to energy
    ene += -kbt_ * std::log(psi*(1.0-post)+(1.0-psi)*post);
  }
  // add Jeffrey's prior on psi
  ene += kbt_ * std::log(psi);
  return ene;
}

// used to update sigma 
double CrossLinkMSRestraint::getEnergy(double sigma, std::vector<double> psi)
{
  // calculate stuff
  double sig  = sq2 * sigma;
  double sig2 = sig * sig;
  // counter for ndata
  int imax = ndata_[0];
  int ipsi = 0;
  // energy
  double ene = 0.0;
  // cycle on arguments 
  for(unsigned i=0;i<getNumberOfArguments();++i){
    // get distance
    double dist = getArgument(i);
    // get various useful stuff
    double log_eLpR_2  = -(dist+length_)*(dist+length_)/(2.0*sig2);
    double eLpR_2      = std::exp(log_eLpR_2);
    double log_e2LR    = 2.0 * length_ * dist / sig2;
    double eLpR_2_e2LR = std::exp(log_eLpR_2+log_e2LR);
    double erfLmR = erf((length_-dist)/sq2/sig);
    double erfLpR = erf((length_+dist)/sq2/sig);
    // calculate posterior
    double post = -sig/(sq2Pi*dist) * (eLpR_2_e2LR - eLpR_2) + 0.5 * (erfLmR + erfLpR);
    // and add to energy
    if(i >= imax){
     ipsi += 1;
     imax += ndata_[ipsi];
    }
    ene += -kbt_ * std::log(psi[ipsi]*(1.0-post)+(1.0-psi[ipsi])*post);
  }
  return ene;
}

double CrossLinkMSRestraint::proposeMove(double x, double xmin, double xmax, double dxmax)
{
 double r = static_cast<double>(rand()) / RAND_MAX;
 double dx = -dxmax + r * 2.0 * dxmax;
 double x_new = x + dx;
 // check boundaries
 if(x_new > xmax){x_new = 2.0 * xmax - x_new;}
 if(x_new < xmin){x_new = 2.0 * xmin - x_new;}
 return x_new;
}

bool CrossLinkMSRestraint::doAccept(double oldE, double newE){
  bool accept = false;
  // calculate delta energy 
  double delta = ( newE - oldE ) / kbt_;
  // if delta is negative always accept move
  if( delta < 0.0 ){ 
   accept = true;
  }else{
   // otherwise extract random number   
   double s = static_cast<double>(rand()) / RAND_MAX;
   if( s < exp(-delta) ) { accept = true; }
  }
  return accept;
}
    
void CrossLinkMSRestraint::doMonteCarlo(){
 double oldE, newE;
 bool accept;
 // cycle on MC steps 
 for(unsigned i=0;i<MCsteps_;++i){
  // get old energy
  oldE = getEnergy(sigma_, psi_);
  // propose move in sigma
  double new_sigma = proposeMove(sigma_,sigma_min_,sigma_max_,Dsigma_);
  // calculate new energy
  newE = getEnergy(new_sigma, psi_);
  // accept or reject
  accept = doAccept(oldE, newE);
  if(accept){
   sigma_ = new_sigma;
   MCaccsig_++;
  }
  // propose moves in psi
  int imin = 0; 
  for(unsigned j=0; j<psi_.size(); ++j){
      // get old energy
      oldE = getEnergy(sigma_, psi_[j], imin, imin+ndata_[j]);
      double new_psi = proposeMove(psi_[j],psi_min_,psi_max_,Dpsi_);
      // calculate new energy
      newE = getEnergy(sigma_, new_psi, imin, imin+ndata_[j]);
      // accept or reject
      accept = doAccept(oldE, newE);
      if(accept){
       psi_[j] = new_psi;
       MCaccpsi_[j]++;
      }
      // increment imin
      imin += ndata_[j];
  }
 }
}

void CrossLinkMSRestraint::update(){
  // get time step 
  long int step = getStep();
  // do MC stuff at the right time step
  if(step%MCstride_==0) doMonteCarlo();
  // this is needed when restarting simulations
  if(MCfirst_==-1) MCfirst_=step;
  // calculate acceptance
  double MCtrials = std::floor(static_cast<double>(step-MCfirst_) / static_cast<double>(MCstride_))+1.0;
  // sigma acceptance
  double accsig = static_cast<double>(MCaccsig_) / static_cast<double>(MCsteps_) / MCtrials;
  getPntrToComponent("accsig")->set(accsig);
  // psi acceptances
  for(unsigned i=0; i<MCaccpsi_.size(); ++i){
   double accpsi = static_cast<double>(MCaccpsi_[i]) / static_cast<double>(MCsteps_) / MCtrials;
   std::string num; Tools::convert(i,num);
   getPntrToComponent("accpsi"+num)->set(accpsi);
  }
}

void CrossLinkMSRestraint::calculate(){
  // calculate stuff
  double sig  = sq2 * sigma_;
  double sig2 = sig * sig;
  // counter for ndata
  int imax = ndata_[0];
  int ipsi = 0;
  // energy
  double ene = 0.0;
  // cycle on arguments
  for(unsigned i=0;i<getNumberOfArguments();++i){
    // get distance
    double dist = getArgument(i);
    // get various useful stuff
    double log_eLpR_2  = -(dist+length_)*(dist+length_)/(2.0*sig2);
    double eLpR_2      = std::exp(log_eLpR_2);
    double log_e2LR    = 2.0 * length_ * dist / sig2;
    double eLpR_2_e2LR = std::exp(log_eLpR_2+log_e2LR);
    double erfLmR = erf((length_-dist)/sq2/sig);
    double erfLpR = erf((length_+dist)/sq2/sig);
    // calculate posterior
    double post = -sig/(sq2Pi*dist) * (eLpR_2_e2LR - eLpR_2) + 0.5 * (erfLmR + erfLpR);
    // and add to energy
    if(i >= imax){
     ipsi += 1;
     imax += ndata_[ipsi];
    }
    ene += -kbt_ * std::log(psi_[ipsi]*(1.0-post)+(1.0-psi_[ipsi])*post) + kbt_ * slope_ * dist;
    // add entropic contribution
    if(entropic_) ene += 2.0*std::log(dist);
    // set force on the i-th component
    double dlog_eLpR_2 = -2.0 * (dist+length_) / (2.0*sig2);
    double dlog_e2LR = 2.0 * length_ / sig2;
    double deLpR_2_e2LR = eLpR_2_e2LR * (dlog_eLpR_2 + dlog_e2LR);
    double deLpR_2 = eLpR_2 * dlog_eLpR_2;
    double derfLmR = -2.0 / sqPi * exp(-(length_-dist)*(length_-dist)/2.0/sig2) / sq2 / sig;
    double derfLpR =  2.0 / sqPi * exp(-(length_+dist)*(length_+dist)/2.0/sig2) / sq2 / sig;
    double dpost = sig/sq2Pi/dist/dist * (eLpR_2_e2LR - eLpR_2) -sig/(sq2Pi*dist) * (deLpR_2_e2LR - deLpR_2) + 0.5 * (derfLmR + derfLpR);
    double force = kbt_ * 1.0 / (psi_[ipsi]*(1.0-post)+(1.0-psi_[ipsi])*post) * (1.0-2.0*psi_[ipsi]) * dpost - kbt_ * slope_;
    if(entropic_) force += -2.0 / dist;
    setOutputForce(i, force);
  };
  // add Jeffrey's priors on psi_
  for(unsigned i=0;i<psi_.size();++i) ene += kbt_ * std::log(psi_[i]);
  // set value of the bias
  getPntrToComponent("bias")->set(ene);
  // set value of uncertainty
  getPntrToComponent("sigma")->set(sigma_);
  // set value of psi 
  for(unsigned i=0;i<psi_.size();++i){ 
     std::string num; Tools::convert(i,num);
     getPntrToComponent("psi"+num)->set(psi_[i]);
  }
}


}
}


