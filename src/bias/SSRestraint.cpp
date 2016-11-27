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
#include <iostream>


using namespace std;

namespace PLMD{
namespace bias{

//+PLUMEDOC BIAS SS
/*


*/
//+ENDPLUMEDOC


class SSRestraint : public Bias
{
  // psi parameter
  double psi_;
  double psi_min_;
  double psi_max_;
  double Dpsi_;
  // temperature in kbt
  double kbt_;
  // reference dihedrals and list
  std::vector<double> phi0_;
  std::vector<double> psi0_;
  std::vector<unsigned> dlist_;
  // Monte Carlo stuff
  int MCsteps_;
  int MCstride_;
  unsigned int MCaccpsi_;
  long int MCfirst_;

  void set_reference_dihedrals(std::string Spred);
  void doMonteCarlo(long int step);
  double getEnergy(double psi);
  double proposeMove(double x, double xmin, double xmax, double dxmax);
  bool doAccept(double oldE, double newE);

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
  keys.add("compulsory","PSI_MIN","minimum value of the psi parameter");
  keys.add("compulsory","PSI_MAX","maximum value of the psi parameter");
  keys.add("compulsory","DPSI","maximum MC move of the psi parameter");
  keys.add("compulsory","SPRED","secondary structure prediction");
  keys.add("compulsory","TEMP","temperature in energy units");
  keys.add("optional","MC_STEPS","number of MC steps");
  keys.add("optional","MC_STRIDE","MC stride");
  componentsAreNotOptional(keys);
  keys.addOutputComponent("psi",   "default","psi parameter");
  keys.addOutputComponent("accpsi","default","MC acceptance psi");
}

SSRestraint::SSRestraint(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao),
MCsteps_(1), MCstride_(1),
MCaccpsi_(0), MCfirst_(-1)
{
  // psi stuff
  parse("PSI0",     psi_);
  parse("PSI_MIN",  psi_min_);
  parse("PSI_MAX",  psi_max_);
  parse("DPSI",     Dpsi_);
  // secondary structure
  std::string Spred;
  parse("SPRED", Spred);
  // compare length of string and number of arguments
  if(Spred.length() != getNumberOfArguments()/2)
    error("you must provide a SS prediction per residue");
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

  // set reference dihedrals
  set_reference_dihedrals(Spred);

  // adjust for multiple-time steps
  MCstride_ *= getStride();
  
  log.printf("  initial value of psi %f\n",psi_);
  log.printf("  minimum value of psi %f\n",psi_min_);
  log.printf("  maximum value of psi %f\n",psi_max_);
  log.printf("  maximum MC move of the psi parameter %f\n",Dpsi_);
  log.printf("  temperature of the system in energy unit %f\n",kbt_);
  log.printf("  number of MC steps %d\n",MCsteps_);
  log.printf("  do MC every %d steps\n", MCstride_);

  addComponent("psi");    componentIsNotPeriodic("psi");
  addComponent("accpsi"); componentIsNotPeriodic("accpsi");
  
  // initialize random seed
  srand (time(NULL));

}

// set the reference dihedrals from SS prediction
void SSRestraint::set_reference_dihedrals(std::string Spred)
{
  char c;
  // split string into characters
  for(unsigned i=0; i<Spred.length(); ++i){
    c = Spred[i];
    // assign dihedrals
    if(c == 'H'){
      phi0_.push_back(-1.05);
      psi0_.push_back(-0.79);
      dlist_.push_back(i);
    }
    if(c == 'E'){
      phi0_.push_back(-2.36);
      psi0_.push_back(2.36);
      dlist_.push_back(i);
    }
  }
}

// used to update Bayesian parameters
double SSRestraint::getEnergy(double psi)
{
  // calculate energy
  double ene = 0.0;
  for(unsigned i=0;i<dlist_.size();++i){
    // index of residue
    unsigned index = dlist_[i];
    // get dihedrals
    double phid = getArgument(2*index);
    double psid = getArgument(2*index+1);
    // calculate probability
    double p = 0.25*(1.0+cos(phid-phi0_[i]))*(1.0+cos(psid-psi0_[i]));
    // add to energy
    ene += -kbt_ * std::log( 0.5*p*(1.0-psi*psi)+psi*(1.0-p)*(1.0-psi) );
  }
  // add prior
  ene += kbt_ * 0.5 * std::log(psi);
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
    
void SSRestraint::doMonteCarlo(long int step)
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

void SSRestraint::calculate(){

  // get time step 
  long int step = getStep();
  // do MC stuff at the right time step
  if(step%MCstride_==0&&!getExchangeStep()) doMonteCarlo(step);
  
  // energy
  double ene = 0.0;
  // cycle on arguments
  for(unsigned i=0;i<dlist_.size();++i){
    // index of residue
    unsigned index = dlist_[i];
    // get dihedrals
    double phid = getArgument(2*index);
    double psid = getArgument(2*index+1);
    // calculate probability
    double p = 0.25*(1.0+cos(phid-phi0_[i]))*(1.0+cos(psid-psi0_[i]));
    // and posterior
    double post = 0.5*p*(1.0-psi_*psi_)+psi_*(1.0-p)*(1.0-psi_);
    // add to energy
    ene += -kbt_ * std::log(post);
    // calculate force
    double dene_dpost = -kbt_ / post;
    double dpost_dp   = 0.5*(1.0-psi_*psi_)-psi_*(1.0-psi_);
    double dp_darg0   = -0.25 * sin(phid-phi0_[i]) * (1.0+cos(psid-psi0_[i]));
    double dp_darg1   = -0.25 * (1.0+cos(phid-phi0_[i])) * sin(psid-psi0_[i]);
    double force0 = -dene_dpost * dpost_dp * dp_darg0;
    double force1 = -dene_dpost * dpost_dp * dp_darg1;
    setOutputForce(2*index,   force0);
    setOutputForce(2*index+1, force1);
  }
  // add prior
  ene += kbt_ * 0.5 * std::log(psi_);
  
  // set value of the bias
  setBias(ene);
}


}
}


