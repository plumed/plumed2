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
#include "bias/Bias.h"
#include "bias/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include <math.h>

using namespace std;

namespace PLMD {
namespace isdb {

//+PLUMEDOC BIAS COEVOLUTION
/*


*/
//+ENDPLUMEDOC


class CoEvolutionRestraint : public bias::Bias
{
  // alpha parameters
  double alpha_p_;
  double alpha_n_;
  double Dalpha_;
  // number of positives/negatives
  unsigned npos_;
  unsigned nneg_;
  // alpha prior parameters
  double alpha_p_mean_;
  double alpha_p_sig_;
  double alpha_n_mean_;
  double alpha_n_sig_;
  // "cutoff", gamma, and P0 parameters;
  double R0_;
  double gamma_;
  double P0_;
  // temperature in kbt
  double kbt_;
  // Monte Carlo stuff
  int MCsteps_;
  int MCstride_;
  unsigned int MCaccalpha_p_;
  unsigned int MCaccalpha_n_;
  long int MCfirst_;
  // parallel stuff
  unsigned rank_;
  unsigned nrep_;

  void   setup_restraint();
  double getPrior(double a, double amean, double asig);
  void   doMonteCarlo(double oldE, const vector<double> &fmod, long int step);
  double getEnergy(double alpha_p, double alpha_n, const vector<double> &fmod);
  double proposeMove(double x, double xmin, double xmax, double dxmax);
  bool   doAccept(double oldE, double newE);

public:
  CoEvolutionRestraint(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(CoEvolutionRestraint,"COEVOLUTION")

void CoEvolutionRestraint::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","ALPHA_P0","initial value of the alpha positive parameter");
  keys.add("compulsory","ALPHA_N0","initial value of the alpha negative parameter");
  keys.add("compulsory","DALPHA","maximum MC move of the alpha parameter");
  keys.add("compulsory","NPOS","number of positives");
  keys.add("compulsory","NNEG","number of negatives");
  keys.add("optional","TEMP","temperature in energy units");
  keys.add("optional","MC_STEPS","number of MC steps");
  keys.add("optional","MC_STRIDE","MC stride");
  componentsAreNotOptional(keys);
  useCustomisableComponents(keys);
  keys.addOutputComponent("alpha_p",   "default","alpha positive parameter");
  keys.addOutputComponent("alpha_n",   "default","alpha negative parameter");
  keys.addOutputComponent("accalpha_p","default","MC acceptance alpha positive");
  keys.addOutputComponent("accalpha_n","default","MC acceptance alpha negative");
}

CoEvolutionRestraint::CoEvolutionRestraint(const ActionOptions&ao):
  PLUMED_BIAS_INIT(ao),
  MCsteps_(1), MCstride_(1), MCaccalpha_p_(0),
  MCaccalpha_n_(0), MCfirst_(-1)
{
  // alpha stuff
  parse("ALPHA_P0", alpha_p_);
  parse("ALPHA_N0", alpha_n_);
  parse("DALPHA",   Dalpha_);
  // number of positives and negatives
  parse("NPOS", npos_);
  if(npos_<=0) error("NPOS should be strictly positive");
  parse("NNEG", nneg_);
  if(nneg_<=0) error("NNEG should be strictly positive");
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

  // check number of arguments
  if(getNumberOfArguments()!=(npos_+nneg_)) error("The number of arguments should be equal to NPOS + NNEG");

  // prepare stuff
  setup_restraint();

  // adjust for multiple-time steps
  MCstride_ *= getStride();

  log.printf("  initial value of alpha_p %f\n",alpha_p_);
  log.printf("  initial value of alpha_n %f\n",alpha_n_);
  log.printf("  maximum MC move of the alpha parameter %f\n",Dalpha_);
  log.printf("  number of positive data points %d\n",npos_);
  log.printf("  number of negative data points %d\n",nneg_);
  log.printf("  temperature of the system in energy unit %f\n",kbt_);
  log.printf("  number of MC steps %d\n",MCsteps_);
  log.printf("  do MC every %d steps\n", MCstride_);

  addComponent("alpha_p");    componentIsNotPeriodic("alpha_p");
  addComponent("accalpha_p"); componentIsNotPeriodic("accalpha_p");
  addComponent("alpha_n");    componentIsNotPeriodic("alpha_n");
  addComponent("accalpha_n"); componentIsNotPeriodic("accalpha_n");

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

// setup the restraint
void CoEvolutionRestraint::setup_restraint()
{
  // set up parameters for the forward models
  R0_ = 0.0;
  P0_ = 0.0;
  gamma_ = 0.0;
  // alpha prior parameters
  alpha_p_mean_ = 0.0;
  alpha_p_sig_ = 0.0;
  alpha_n_mean_ = 0.0;
  alpha_n_sig_ = 0.0;
}

// calculate Gaussian prior for a single alpha
double CoEvolutionRestraint::getPrior(double a, double amean, double asig)
{
  // calculate sigma
  double s = sqrt(a*(1.0-a)/asig);
  // calculate prior - excluding 2pi which is constant
  double prior = 0.5 * ( a - amean ) * ( a - amean ) / s / s + std::log(s);

  return kbt_ * prior;
}

// used to update Bayesian nuisance parameters
double CoEvolutionRestraint::getEnergy
(double alpha_p, double alpha_n, const vector<double> &fmod)
{
  // calculate energy
  double ene = 0.0;
  // cycle on positive arguments
  for(unsigned i=rank_; i<npos_; i=i+nrep_) {
    // calculate data likelihood
    double like = alpha_p * fmod[i] + alpha_n * ( 1.0 - fmod[i] );
    // add to energy
    ene += -kbt_ * std::log(like);
  }
  // cycle on negative arguments
  for(unsigned i=rank_; i<nneg_; i=i+nrep_) {
    // calculate data likelihood
    double like = ( 1.0 - alpha_p ) * fmod[i+npos_] + ( 1.0 - alpha_n ) * ( 1.0 - fmod[i+npos_] );
    // add to energy
    ene += -kbt_ * std::log(like);
  }

  // sum energy
  comm.Sum(&ene, 1);

  // add prior on alpha_p
  ene += getPrior(alpha_p, alpha_p_mean_, alpha_p_sig_);
  // add prior on alpha_n
  ene += getPrior(alpha_n, alpha_n_mean_, alpha_n_sig_);

  return ene;
}

double CoEvolutionRestraint::proposeMove(double x, double xmin, double xmax, double dxmax)
{
  double r = static_cast<double>(rand()) / RAND_MAX;
  double dx = -dxmax + r * 2.0 * dxmax;
  double x_new = x + dx;
// check boundaries
  if(x_new > xmax) {x_new = 2.0 * xmax - x_new;}
  if(x_new < xmin) {x_new = 2.0 * xmin - x_new;}
  return x_new;
}

bool CoEvolutionRestraint::doAccept(double oldE, double newE) {
  bool accept = false;
  // calculate delta energy
  double delta = ( newE - oldE ) / kbt_;
  // if delta is negative always accept move
  if( delta < 0.0 ) {
    accept = true;
  } else {
    // otherwise extract random number
    double s = static_cast<double>(rand()) / RAND_MAX;
    if( s < exp(-delta) ) { accept = true; }
  }
  return accept;
}

void CoEvolutionRestraint::doMonteCarlo(double oldE, const vector<double> &fmod, long int step)
{
// cycle on MC steps
  for(unsigned i=0; i<MCsteps_; ++i) {
    // 1) propose move in alpha_p
    double new_alpha_p = proposeMove(alpha_p_, 0.0, 1.0, Dalpha_);
    // calculate new energy
    double newE = getEnergy(new_alpha_p, alpha_n_, fmod);
    // accept or reject
    bool accept = doAccept(oldE, newE);
    if(accept) {
      alpha_p_ = new_alpha_p;
      MCaccalpha_p_++;
      oldE = newE;
    }
    // 2) propose move in alpha_n
    double new_alpha_n = proposeMove(alpha_n_, 0.0, 1.0, Dalpha_);
    // calculate new energy
    newE = getEnergy(alpha_p_, new_alpha_n, fmod);
    // accept or reject
    accept = doAccept(oldE, newE);
    if(accept) {
      alpha_n_ = new_alpha_n;
      MCaccalpha_n_++;
      oldE = newE;
    }

  }
// this is needed when restarting simulations
  if(MCfirst_==-1) MCfirst_=step;
// calculate number of trials
  double MCtrials = std::floor(static_cast<double>(step-MCfirst_) / static_cast<double>(MCstride_))+1.0;
// alpha_p acceptance
  double accalpha_p = static_cast<double>(MCaccalpha_p_) / static_cast<double>(MCsteps_) / MCtrials;
  getPntrToComponent("accalpha_p")->set(accalpha_p);
  // alpha_n acceptance
  double accalpha_n = static_cast<double>(MCaccalpha_n_) / static_cast<double>(MCsteps_) / MCtrials;
  getPntrToComponent("accalpha_n")->set(accalpha_n);
}

void CoEvolutionRestraint::calculate()
{
  // allocate force vector
  vector<double> force(getNumberOfArguments(), 0.0);
  // and forward model vector
  vector<double> fmod(getNumberOfArguments(), 0.0);

  // calculate energy
  double ene = 0.0;
  // cycle on positive arguments
  for(unsigned i=rank_; i<npos_; i=i+nrep_) {
    // get distance
    double dist = getArgument(i);
    // calculate forward model
    double tmp = exp(-gamma_*(dist-R0_));
    double p = P0_ * ( 1.0 - 1.0 / (1.0+tmp) );
    // store forward model
    fmod[i] = p;
    // calculate data likelihood
    double like = alpha_p_ * p + alpha_n_ * ( 1.0 - p );
    // add to energy
    ene += -kbt_ * std::log(like);
    // calculate force
    double dene_dlike = -kbt_ / like;
    double dlike_dp   = alpha_p_ - alpha_n_;
    double dp_ddist   = -P0_ / (1.0+tmp) / (1.0+tmp) * tmp * gamma_;
    // apply chain rule
    force[i] = -dene_dlike * dlike_dp * dp_ddist;
  }
  // cycle on negative arguments
  for(unsigned i=rank_; i<nneg_; i=i+nrep_) {
    // get distance
    double dist = getArgument(i+npos_);
    // calculate forward model
    double tmp = exp(-gamma_*(dist-R0_));
    double p = P0_ * ( 1.0 - 1.0 / (1.0+tmp) );
    // store forward model
    fmod[i+npos_] = p;
    // calculate data likelihood
    double like = ( 1.0 - alpha_p_ ) * p + ( 1.0 - alpha_n_ ) * ( 1.0 - p );
    // add to energy
    ene += -kbt_ * std::log(like);
    // calculate force
    double dene_dlike = -kbt_ / like;
    double dlike_dp   = - alpha_p_ + alpha_n_;
    double dp_ddist   = -P0_ / (1.0+tmp) / (1.0+tmp) * tmp * gamma_;
    // apply chain rule
    force[i+npos_]  = -dene_dlike * dlike_dp * dp_ddist;
  }

  // sum energy, fmod, and derivatives
  comm.Sum(&force[0], force.size());
  comm.Sum(&fmod[0], fmod.size());
  comm.Sum(&ene, 1);

  // apply forces
  for(unsigned i=0; i<force.size(); ++i) setOutputForce(i, force[i]);

// add prior on alpha_p
  ene += getPrior(alpha_p_, alpha_p_mean_, alpha_p_sig_);
  // add prior on alpha_n
  ene += getPrior(alpha_n_, alpha_n_mean_, alpha_n_sig_);

  // set value of the bias
  setBias(ene);
  // set values of alpha_p
  getPntrToComponent("alpha_p")->set(alpha_p_);
  // set values of alpha_n
  getPntrToComponent("alpha_n")->set(alpha_n_);

  // get time step
  long int step = getStep();
  // do MC stuff at the right time step
  if(step%MCstride_==0&&!getExchangeStep()) doMonteCarlo(ene, fmod, step);

}


}
}


