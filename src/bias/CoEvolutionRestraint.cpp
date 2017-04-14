/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017 The plumed team
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
#include <math.h>

using namespace std;

namespace PLMD{
namespace bias{

//+PLUMEDOC BIAS COEVOLUTION
/*


*/
//+ENDPLUMEDOC


class CoEvolutionRestraint : public Bias
{
  // psi parameter - true positive discovery rate
  double psi_;
  double Dpsi_;
  // phi parameter - true negative discovery rate
  double phi_;
  double Dphi_;
  // number of positives/negatives
  unsigned npos_;
  unsigned nneg_;
  // psi/phi prior parameters
  double psi_mean_;
  double psi_sig_;
  double phi_mean_;
  double phi_sig_; 
  // cutoff and slope parameters;
  double R0_;
  double alpha_;
  // temperature in kbt
  double kbt_;
  // Monte Carlo stuff
  int MCsteps_;
  int MCstride_;
  unsigned int MCaccpsi_;
  unsigned int MCaccphi_;
  long int MCfirst_;
  // parallel stuff
  unsigned rank_;
  unsigned nrep_;
  
  double getPrior(double p, double pmean, double psig);
  void   doMonteCarlo(double oldE, long int step);
  double getEnergy(double psi, double phi);
  double proposeMove(double x, double xmin, double xmax, double dxmax);
  bool   doAccept(double oldE, double newE);

public:
  CoEvolutionRestraint(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(CoEvolutionRestraint,"COEVOLUTION")

void CoEvolutionRestraint::registerKeywords(Keywords& keys){
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","PSI0","initial value of the psi parameter");
  keys.add("compulsory","PHI0","initial value of the phi parameter");
  keys.add("compulsory","DPSI","maximum MC move of the psi parameter");
  keys.add("compulsory","DPHI","maximum MC move of the phi parameter");
  keys.add("compulsory","PSI_MEAN","psi prior parameter - average value");
  keys.add("compulsory","PSI_SIGMA","psi prior parameter - standard deviation");
  keys.add("compulsory","PHI_MEAN","phi prior parameter - average value");
  keys.add("compulsory","PHI_SIGMA","phi prior parameter - standard deviation");
  keys.add("compulsory","NPOS","number of positives");
  keys.add("compulsory","NNEG","number of negatives");
  keys.add("compulsory","R0","Value of the R0 parameter");
  keys.add("optional","ALPHA","Value of the alpha parameter");
  keys.add("optional","TEMP","temperature in energy units");
  keys.add("optional","MC_STEPS","number of MC steps");
  keys.add("optional","MC_STRIDE","MC stride");
  componentsAreNotOptional(keys);
  useCustomisableComponents(keys);
  keys.addOutputComponent("psi",   "default","psi parameter");
  keys.addOutputComponent("phi",   "default","phi parameter");
  keys.addOutputComponent("accpsi","default","MC acceptance psi");
  keys.addOutputComponent("accphi","default","MC acceptance phi");
}

CoEvolutionRestraint::CoEvolutionRestraint(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao), alpha_(1.0),
MCsteps_(1), MCstride_(1), MCaccpsi_(0),
MCaccphi_(0), MCfirst_(-1)
{
  // parse positive arguments
  vector<Value*> parg;
  parseArgumentList("ARG",1,parg);
  // parse negative arguments
  vector<Value*> narg;
  parseArgumentList("ARG",2,narg);
  // merge lists into global list
  vector<Value*> arg;
  for(unsigned i=0; i<parg.size(); ++i) arg.push_back(parg[i]);
  for(unsigned i=0; i<narg.size(); ++i) arg.push_back(narg[i]);
  
  // psi stuff
  parse("PSI0", psi_);
  parse("DPSI", Dpsi_);
  // phi stuff
  parse("PHI0", phi_);
  parse("DPHI", Dphi_);
  // priors parameters
  parse("PSI_MEAN",  psi_mean_);
  parse("PSI_SIGMA", psi_sig_);
  parse("PHI_MEAN",  phi_mean_);
  parse("PHI_SIGMA", phi_sig_);
  // number of positives and negatives
  parse("NPOS", npos_);
  if(npos_<=0) error("NPOS should be strictly positive");
  parse("NNEG", nneg_);
  if(nneg_<=0) error("NNEG should be strictly positive");
  // R0 stuff
  parse("R0", R0_);
  // alpha parameter
  parse("ALPHA", alpha_);
  // temperature
  double temp=0.0;
  parse("TEMP",temp);
  // convert temp to kbt
  if(temp>0.0) kbt_=plumed.getAtoms().getKBoltzmann()*temp;
  else kbt_=plumed.getAtoms().getKbT();
  // MC stuff
  parse("MC_STEPS", MCsteps_);
  parse("MC_STRIDE",MCstride_);

  checkRead();
  
  // control number of positive and negative arguments
  if(parg.size()!=npos_) error("The number of arguments in ARG1 should be equal to NPOS");
  if(narg.size()!=0 && narg.size()!=nneg_)
   error("The number of arguments in ARG2 should be equal to either zero or NNEG");
  
  // ask for arguments
  requestArguments(arg);

  // adjust for multiple-time steps
  MCstride_ *= getStride();

  log.printf("  initial value of psi %f\n",psi_);
  log.printf("  initial value of phi %f\n",psi_);
  log.printf("  maximum MC move of the psi parameter %f\n",Dpsi_);
  log.printf("  maximum MC move of the phi parameter %f\n",Dphi_);
  log.printf("  psi prior average %f\n",psi_mean_);
  log.printf("  psi prior standard deviation %f\n",psi_sig_);
  log.printf("  phi prior average %f\n",phi_mean_);
  log.printf("  phi prior standard deviation %f\n",phi_sig_);
  log.printf("  number of positive data points %d\n",npos_);
  log.printf("  number of negative data points %d\n",nneg_);
  log.printf("  value of R0 parameter %f\n",R0_);
  log.printf("  value of alpha parameter %f\n", alpha_);
  log.printf("  temperature of the system in energy unit %f\n",kbt_);
  log.printf("  number of MC steps %d\n",MCsteps_);
  log.printf("  do MC every %d steps\n", MCstride_);

  addComponent("psi");    componentIsNotPeriodic("psi");
  addComponent("accpsi"); componentIsNotPeriodic("accpsi");
  addComponent("phi");    componentIsNotPeriodic("phi");
  addComponent("accphi"); componentIsNotPeriodic("accphi");
  
  // initialize parallel stuff
  rank_ = comm.Get_rank();
  nrep_ = comm.Get_size();
  
  // initialize random seed
  unsigned iseed;
  if(rank_ == 0) iseed = time(NULL);
  else           iseed = 0;
  comm.Sum(&iseed, 1);
  // initialize random generator
  srand (iseed);
  
}

// get truncated Gaussian prior
double CoEvolutionRestraint::getPrior(double p, double pmean, double psig)
{
  double sqrt2 = sqrt(2.0);
  // calculate normalization
  double phi_B = 0.5 * ( 1.0 + erf ( (1.0 - pmean) / psig / sqrt2 ) );
  double phi_A = 0.5 * ( 1.0 - erf ( pmean / psig / sqrt2 ) );
  double norm = psig * ( phi_B - phi_A );
  // calculate prior
  double eps = ( p - pmean ) / psig;
  double prior = 0.5 * eps * eps + std::log(norm);

  return kbt_ * prior;
}

// used to update Bayesian nuisance parameters
double CoEvolutionRestraint::getEnergy(double psi, double phi)
{
  // ratio of negative and positives data points
  double ratio = static_cast<double>(nneg_)/static_cast<double>(npos_);
  // calculate TPR
  double tpr = 1.0 / ( 1.0 + ratio * ( 1.0 - phi ) / psi );
  // calculate FPR
  double fpr  = 1.0 / ( 1.0 + ratio * phi / (1.0 - psi ) );
  
    // calculate energy
  double ene = 0.0;
  // cycle on positive arguments
  for(unsigned i=rank_;i<npos_;i=i+nrep_){
    // get distance
    double dist = getArgument(i);
    // calculate forward model
    double tmp = exp(-alpha_*(dist-R0_));
    double p = 1.0 - 1.0 / (1.0+tmp);
    // calculate data likelihood
    double like = tpr * p + fpr * ( 1.0 - p );
    // add to energy
    ene += -kbt_ * std::log(like);
  }
  // number of negative arguments (either zero or nneg_)
  unsigned nneg = getNumberOfArguments() - npos_;
  // cycle on negative arguments
  for(unsigned i=rank_;i<nneg;i=i+nrep_){
    // get distance
    double dist = getArgument(i+npos_);
    // calculate forward model
    double tmp = exp(-alpha_*(dist-R0_));
    double p = 1.0 - 1.0 / (1.0+tmp);
    // calculate data likelihood
    double like = ( 1.0 - tpr ) * p + ( 1.0 - fpr ) * ( 1.0 - p );
    // add to energy
    ene += -kbt_ * std::log(like);
  }
  
  // sum energy
  comm.Sum(&ene, 1);
  
  // add prior on psi
  ene += getPrior(psi, psi_mean_, psi_sig_);
  // add prior on phi
  ene += getPrior(phi, phi_mean_, phi_sig_); 

  return ene;
}

double CoEvolutionRestraint::proposeMove(double x, double xmin, double xmax, double dxmax)
{
 double r = static_cast<double>(rand()) / RAND_MAX;
 double dx = -dxmax + r * 2.0 * dxmax;
 double x_new = x + dx;
 // check boundaries
 if(x_new > xmax){x_new = 2.0 * xmax - x_new;}
 if(x_new < xmin){x_new = 2.0 * xmin - x_new;}
 return x_new;
}

bool CoEvolutionRestraint::doAccept(double oldE, double newE){
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
    
void CoEvolutionRestraint::doMonteCarlo(double oldE, long int step)
{ 
 // cycle on MC steps 
 for(unsigned i=0;i<MCsteps_;++i){
   // propose move in psi and phi
   double new_psi = proposeMove(psi_, 0.0, 1.0, Dpsi_);
   double new_phi = proposeMove(phi_, 0.0, 1.0, Dphi_);
   // calculate new energy
   double newE = getEnergy(new_psi, new_phi);
   // accept or reject
   bool accept = doAccept(oldE, newE);
   if(accept){
    psi_ = new_psi;
    phi_ = new_phi;
    MCaccpsi_++;
    MCaccphi_++;
    oldE = newE;
   }
 }
 // this is needed when restarting simulations
 if(MCfirst_==-1) MCfirst_=step;
 // calculate number of trials
 double MCtrials = std::floor(static_cast<double>(step-MCfirst_) / static_cast<double>(MCstride_))+1.0;
 // psi acceptance
 double accpsi = static_cast<double>(MCaccpsi_) / static_cast<double>(MCsteps_) / MCtrials;
 getPntrToComponent("accpsi")->set(accpsi);
 // phi acceptance
 double accphi = static_cast<double>(MCaccphi_) / static_cast<double>(MCsteps_) / MCtrials;
 getPntrToComponent("accphi")->set(accphi);
}

void CoEvolutionRestraint::calculate()
{
  // allocate force vector
  vector<double> force(getNumberOfArguments(), 0.0);

  // ratio of negative and positives data points
  double ratio = static_cast<double>(nneg_)/static_cast<double>(npos_);
  // calculate TPR
  double tpr = 1.0 / ( 1.0 + ratio * ( 1.0 - phi_ ) / psi_ );
  // calculate FPR
  double fpr  = 1.0 / ( 1.0 + ratio * phi_ / (1.0 - psi_ ) );
  
  // calculate energy
  double ene = 0.0;
  // cycle on positive arguments
  for(unsigned i=rank_;i<npos_;i=i+nrep_){
    // get distance
    double dist = getArgument(i);
    // calculate forward model
    double tmp = exp(-alpha_*(dist-R0_));
    double p = 1.0 - 1.0 / (1.0+tmp);
    // calculate data likelihood
    double like = tpr * p + fpr * ( 1.0 - p );
    // add to energy
    ene += -kbt_ * std::log(like);
    // calculate force
    double dene_dlike = -kbt_ / like;
    double dlike_dp   = tpr - fpr;
    double dp_ddist   = -1.0 / (1.0+tmp) / (1.0+tmp) * tmp * alpha_;
    // apply chain rule
    force[i] = -dene_dlike * dlike_dp * dp_ddist;
  }
  // number of negative arguments (either zero or nneg_)
  unsigned nneg = getNumberOfArguments() - npos_;
  // cycle on negative arguments
  for(unsigned i=rank_;i<nneg;i=i+nrep_){
    // get distance
    double dist = getArgument(i+npos_);
    // calculate forward model
    double tmp = exp(-alpha_*(dist-R0_));
    double p = 1.0 - 1.0 / (1.0+tmp);
    // calculate data likelihood
    double like = ( 1.0 - tpr ) * p + ( 1.0 - fpr ) * ( 1.0 - p );
    // add to energy
    ene += -kbt_ * std::log(like);
    // calculate force
    double dene_dlike = -kbt_ / like;
    double dlike_dp   = - tpr + fpr;
    double dp_ddist   = -1.0 / (1.0+tmp) / (1.0+tmp) * tmp * alpha_;
    // apply chain rule
    force[i+npos_]  = -dene_dlike * dlike_dp * dp_ddist;
  }
  
  // sum energy and derivatives
  comm.Sum(&force[0], force.size());
  comm.Sum(&ene, 1);
  
  // apply forces
  for(unsigned i=0; i<force.size(); ++i) setOutputForce(i, force[i]);
  
  // add prior on psi
  ene += getPrior(psi_, psi_mean_, psi_sig_);
  // add prior on phi
  ene += getPrior(phi_, phi_mean_, phi_sig_); 
   
  // set value of the bias
  setBias(ene);
  // set values of psi
  getPntrToComponent("psi")->set(psi_);
  // set values of phi
  getPntrToComponent("phi")->set(phi_);

  // get time step 
  long int step = getStep();
  // do MC stuff at the right time step
  if(step%MCstride_==0&&!getExchangeStep()) doMonteCarlo(ene, step);
   
}


}
}


