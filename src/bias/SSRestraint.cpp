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
#include <math.h>
#include <string.h>


using namespace std;

namespace PLMD{
namespace bias{

//+PLUMEDOC BIAS SS
/*


*/
//+ENDPLUMEDOC


class SSRestraint : public Bias
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
  // reference dihedrals and list
  double phi0_;
  double psi0_;
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

  void set_reference_dihedrals(std::string ss_type);
  double getPrior(double p, double pmean, double psig);
  void   doMonteCarlo(double oldE, long int step);
  double getEnergy(double psi, double phi);
  double proposeMove(double x, double xmin, double xmax, double dxmax);
  bool   doAccept(double oldE, double newE);

public:
  SSRestraint(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(SSRestraint,"SS")

void SSRestraint::registerKeywords(Keywords& keys){
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
  keys.add("compulsory","NRES","number of residues");
  keys.add("compulsory","TYPE","type of secondary structure (currently H, E, or P)");
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

SSRestraint::SSRestraint(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao),
MCsteps_(1), MCstride_(1), MCaccpsi_(0),
MCaccphi_(0), MCfirst_(-1)
{
  // parse phi dihedrals
  vector<Value*> phi_arg;
  parseArgumentList("ARG",1,phi_arg);
  // parse psi dihedrals
  vector<Value*> psi_arg;
  parseArgumentList("ARG",2,psi_arg);
  // check length
  if(phi_arg.size()!=psi_arg.size()) error("The number of arguments in ARG1 and ARG2 should be the same");
  // and set number of positives
  npos_ = phi_arg.size();
  // merge lists into global list
  vector<Value*> arg;
  for(unsigned i=0; i<phi_arg.size(); ++i) arg.push_back(phi_arg[i]);
  for(unsigned i=0; i<psi_arg.size(); ++i) arg.push_back(psi_arg[i]);

  // number of residues
  unsigned nres;
  parse("NRES", nres);
  if(nres<=0) error("NRES should be strictly positive");
  // and number of negatives
  nneg_ = nres - npos_;
  if(nneg_<=0) error("The number of negatives should be strictly positive");
  // secondary structure type
  string ss_type;
  parse("TYPE", ss_type);
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
  
  // ask for arguments
  requestArguments(arg);

  // set reference dihedrals
  set_reference_dihedrals(ss_type);

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

// set the reference dihedrals from SS prediction
void SSRestraint::set_reference_dihedrals(std::string ss_type)
{
  // assign dihedrals
  // right-handed HELIX
  if(ss_type.compare("H")==0){
      phi0_= -1.05;
      psi0_= -0.79;
  } else if(ss_type.compare("E")==0) {
  // beta strand
      phi0_ = -2.36;
      psi0_ =  2.36;
  } else if(ss_type.compare("P")==0){
  // polyproline
      phi0_ = -1.31;
      psi0_ =  2.71 ;
  } else {
      error("Unsupported secondary structure type");
  }
}

// get truncated Gaussian prior
double SSRestraint::getPrior(double p, double pmean, double psig)
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
double SSRestraint::getEnergy(double psi, double phi)
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
    // get dihedrals
    double phid = getArgument(i);
    double psid = getArgument(i+npos_);
    // calculate forward model
    double p = 0.25*(1.0+cos(phid-phi0_))*(1.0+cos(psid-psi0_));
    // calculate data likelihood
    double like = tpr * p + fpr * ( 1.0 - p );
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

double SSRestraint::proposeMove(double x, double xmin, double xmax, double dxmax)
{
 double r = static_cast<double>(rand()) / RAND_MAX;
 double dx = -dxmax + r * 2.0 * dxmax;
 double x_new = x + dx;
 // check boundaries
 if(x_new > xmax){x_new = 2.0 * xmax - x_new;}
 if(x_new < xmin){x_new = 2.0 * xmin - x_new;}
 return x_new;
}

bool SSRestraint::doAccept(double oldE, double newE){
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
    
void SSRestraint::doMonteCarlo(double oldE, long int step)
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
  
void SSRestraint::calculate()
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
  // cycle on dihedrals
  for(unsigned i=rank_;i<npos_;i=i+nrep_){
    // get dihedrals
    double phid = getArgument(i);
    double psid = getArgument(i+npos_);
    // calculate forward model
    double p = 0.25*(1.0+cos(phid-phi0_))*(1.0+cos(psid-psi0_));
    // calculate data likelihood
    double like = tpr * p + fpr * ( 1.0 - p );
    // add to energy
    ene += -kbt_ * std::log(like);
    // calculate force
    double dene_dlike = -kbt_ / like;
    double dlike_dp   = tpr - fpr;
    double dp_dphid   = -0.25 * sin(phid-phi0_) * (1.0+cos(psid-psi0_));
    double dp_dpsid   = -0.25 * (1.0+cos(phid-phi0_)) * sin(psid-psi0_);
    // apply chain rule
    force[i]       = -dene_dlike * dlike_dp * dp_dphid;
    force[i+npos_] = -dene_dlike * dlike_dp * dp_dpsid;
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


