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

using namespace std;

namespace PLMD{
namespace bias{

//+PLUMEDOC BIAS COEVOLUTION
/*


*/
//+ENDPLUMEDOC


class CoEvolutionRestraint : public Bias
{
  // psi parameter - for non marginal version
  double psi_;
  double psi_min_;
  double psi_max_;
  double Dpsi_;
  // cutoff and slope parameters;
  double R0_;
  double alpha_;
  // temperature in kbt
  double kbt_;
  // Monte Carlo stuff
  int MCsteps_;
  int MCstride_;
  unsigned int MCaccpsi_;
  long int MCfirst_;
  // marginal posterior
  bool marginal_;
  
  void doMonteCarlo(long int step);
  double getEnergy(double psi);
  double proposeMove(double x, double xmin, double xmax, double dxmax);
  bool doAccept(double oldE, double newE);

public:
  CoEvolutionRestraint(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(CoEvolutionRestraint,"COEVOLUTION")

void CoEvolutionRestraint::registerKeywords(Keywords& keys){
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("optional","PSI0","initial value of the psi parameter");
  keys.add("optional","PSI_MIN","minimum value of the psi parameter");
  keys.add("optional","PSI_MAX","maximum value of the psi parameter");
  keys.add("optional","DPSI","maximum MC move of the psi parameter");
  keys.add("compulsory","R0","Value of the R0 parameter");
  keys.add("optional","ALPHA","Value of the alpha parameter");
  keys.add("compulsory","TEMP","temperature in energy units");
  keys.add("optional","MC_STEPS","number of MC steps");
  keys.add("optional","MC_STRIDE","MC stride");
  keys.addFlag("MARGINAL",false,"use marginal version of posterior");
  componentsAreNotOptional(keys);
  keys.addOutputComponent("psi",   "MARGINAL","psi parameter");
  keys.addOutputComponent("accpsi","MARGINAL","MC acceptance psi");
}

CoEvolutionRestraint::CoEvolutionRestraint(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao), alpha_(1.0),
MCsteps_(1), MCstride_(1),
MCaccpsi_(0), MCfirst_(-1),
marginal_(false)
{
  // psi stuff
  parse("PSI0",     psi_);
  parse("PSI_MIN",  psi_min_);
  parse("PSI_MAX",  psi_max_);
  parse("DPSI",     Dpsi_);
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
  // marginal version
  parseFlag("MARGINAL",marginal_);

  checkRead();

  // adjust for multiple-time steps
  MCstride_ *= getStride();

  if(!marginal_){
   log.printf("  initial value of psi %f\n",psi_);
   log.printf("  minimum value of psi %f\n",psi_min_);
   log.printf("  maximum value of psi %f\n",psi_max_);
   log.printf("  maximum MC move of the psi parameter %f\n",Dpsi_);
  } else {
   log.printf("  using marginal version\n"); 
  }
  log.printf("  value of R0 parameter %f\n",R0_);
  log.printf("  value of alpha parameter %f\n", alpha_);
  log.printf("  temperature of the system in energy unit %f\n",kbt_);
  log.printf("  number of MC steps %d\n",MCsteps_);
  log.printf("  do MC every %d steps\n", MCstride_);

  if(!marginal_){
   addComponent("psi");    componentIsNotPeriodic("psi");
   addComponent("accpsi"); componentIsNotPeriodic("accpsi");
  }
  
  // initialize random seed
  srand (time(NULL));

}

// used to update Bayesian parameters - non marginal version
double CoEvolutionRestraint::getEnergy(double psi)
{
  // calculate energy
  double ene = 0.0;
  for(unsigned i=0;i<getNumberOfArguments();++i){
    // get distance
    double dist = getArgument(i);
    // calculate probability
    double p = 1.0 - 1.0 / (1.0+exp(-alpha_*(dist-R0_)));
    // add to energy
    ene += -kbt_ * std::log( 0.5*p*(1.0-psi*psi)+psi*(1.0-p)*(1.0-psi) );
  }
  // add prior
  ene += kbt_ * 0.5 * std::log(psi);
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
    
void CoEvolutionRestraint::doMonteCarlo(long int step)
{
 double oldE, newE;
 bool accept;
 // this is needed when restarting simulations
 if(MCfirst_==-1) MCfirst_=step;
 // calculate acceptance
 double MCtrials = std::floor(static_cast<double>(step-MCfirst_) / static_cast<double>(MCstride_))+1.0;
 // store old energy
 oldE = getEnergy(psi_);
 // cycle on MC steps 
 for(unsigned i=0;i<MCsteps_;++i){
   // propose move in psi
   double new_psi = proposeMove(psi_,psi_min_,psi_max_,Dpsi_);
   // calculate new energy
   newE = getEnergy(new_psi);
   // accept or reject
   accept = doAccept(oldE, newE);
   if(accept){
    psi_ = new_psi;
    MCaccpsi_++;
    oldE = newE;
   }
 }
 // set values of psi
 getPntrToComponent("psi")->set(psi_);
 // psi acceptance
 double accpsi = static_cast<double>(MCaccpsi_) / static_cast<double>(MCsteps_) / MCtrials;
 getPntrToComponent("accpsi")->set(accpsi);
}

void CoEvolutionRestraint::calculate(){

  // get time step 
  long int step = getStep();
  // do MC stuff at the right time step
  if(step%MCstride_==0&&!getExchangeStep()&&!marginal_) doMonteCarlo(step);
  
  // energy
  double ene = 0.0;
  // cycle on arguments
  // non-marginal version
  if(!marginal_){
   for(unsigned i=0;i<getNumberOfArguments();++i){
    // get distance
    double dist = getArgument(i);
    // calculate probability
    double tmp = exp(-alpha_*(dist-R0_));
    double p = 1.0 - 1.0 / (1.0+tmp);
    double post = 0.5*p*(1.0-psi_*psi_)+psi_*(1.0-p)*(1.0-psi_);
    // add to energy
    ene += -kbt_ * std::log(post);
    // calculate force
    double dene_dpost = -kbt_ / post;
    double dpost_dp   = 0.5*(1.0-psi_*psi_)-psi_*(1.0-psi_);
    double dp_ddist   = -1.0 / (1.0+tmp) / (1.0+tmp) * tmp * alpha_;
    double force = -dene_dpost * dpost_dp * dp_ddist;
    setOutputForce(i, force);
   }
   // add prior
   ene += kbt_ * 0.5 * std::log(psi_);
  // marginal version 
  } else {
   for(unsigned i=0;i<getNumberOfArguments();++i){
    // get distance
    double dist = getArgument(i);
    // calculate probability
    double tmp = exp(-alpha_*(dist-R0_));
    double p = 1.0 - 1.0 / (1.0+tmp);
    // add to energy
    ene += -kbt_ * std::log(1.0 + 2.0*p);
    // calculate force
    double dene_dp  = -kbt_ / (1.0 + 2.0*p) * 2.0;
    double dp_ddist = -1.0 / (1.0+tmp) / (1.0+tmp) * tmp * alpha_;
    double force = -dene_dp * dp_ddist;
    setOutputForce(i, force);
   }
  }
  // set value of the bias
  setBias(ene);
}


}
}


