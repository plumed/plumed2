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
#include "bias/Bias.h"
#include "bias/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include <math.h>
#include <string.h>
#include <iostream>

using namespace std;

namespace PLMD {
namespace isdb {

//+PLUMEDOC BIAS SS
/*


*/
//+ENDPLUMEDOC


class SSRestraint : public bias::Bias
{
  // phi parameters
  map <string, double> phi_;
  double Dphi_;
  // phi prior parameters
  map <string, double> phi_mean_;
  map <string, double> phi_sig_;
  // and labels
  const vector<string> phi_label_ = {"HH","HE","HC","EH","EE","EC"};
  // fmod parameters
  map <string, double> phi_m_;
  map <string, double> psi_m_;
  map <string, double> phi_s_;
  map <string, double> psi_s_;
  map <string, double> fmod_a_;
  map <string, double> fmod_b_;
  // secondary structure prediction
  vector<string> ss_pred_;
  // temperature in kbt
  double kbt_;
  // Monte Carlo stuff
  int MCsteps_;
  int MCstride_;
  map <string, unsigned> MCaccphi_;
  long int MCfirst_;
  // parallel stuff
  unsigned rank_;
  unsigned nrep_;

  void   setup_restraint();
  double getPrior(double p, double pmean, double psig);
  double getPriors(map <string, double> ps, map <string, double> pms, map <string, double> pss);
  void   doMonteCarlo(double oldE, long int step, const vector<double> &pHl, const vector<double> &pEl, const vector<double> &pCl);
  double getEnergy(map <string, double> phi, const vector<double> &pHl, const vector<double> &pEl, const vector<double> &pCl);
  double proposeMove(double x, double xmin, double xmax, double dxmax);
  void   proposeMoveCouple(double &x1, double &x2, double dx);
  bool   doAccept(double oldE, double newE);
  map<string, double> get_p_coefficients(map <string, double> phi);

public:
  SSRestraint(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(SSRestraint,"SS")

void SSRestraint::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","PHI_HH0","initial value of the phi_HH parameter");
  keys.add("compulsory","PHI_HE0","initial value of the phi_HE parameter");
  keys.add("compulsory","PHI_HC0","initial value of the phi_HC parameter");
  keys.add("compulsory","PHI_EH0","initial value of the phi_EH parameter");
  keys.add("compulsory","PHI_EE0","initial value of the phi_EE parameter");
  keys.add("compulsory","PHI_EC0","initial value of the phi_EC parameter");
  keys.add("compulsory","DPHI","maximum MC move of the phi parameters");
  keys.add("compulsory","PSIPRED","secondary structure prediction");
  keys.add("optional","TEMP","temperature in energy units");
  keys.add("optional","MC_STEPS","number of MC steps");
  keys.add("optional","MC_STRIDE","MC stride");
  componentsAreNotOptional(keys);
  useCustomisableComponents(keys);
  keys.addOutputComponent("phi",   "default","phi parameter");
  keys.addOutputComponent("accphi","default","MC acceptance phi");
}

SSRestraint::SSRestraint(const ActionOptions&ao):
  PLUMED_BIAS_INIT(ao),
  MCsteps_(1), MCstride_(1), MCfirst_(-1)
{
  // read initial values of the phi parameters
  for(unsigned i=0; i<phi_label_.size(); ++i) parse("PHI_"+phi_label_[i]+"0", phi_[phi_label_[i]]);

  // check consistency
  if((phi_["HH"]+phi_["EH"])>1.0) error("PHI_HH0+PHI_EH0 should be lower than 1");
  if((phi_["HE"]+phi_["EE"])>1.0) error("PHI_HE0+PHI_EE0 should be lower than 1");
  if((phi_["HC"]+phi_["EC"])>1.0) error("PHI_HC0+PHI_EC0 should be lower than 1");

  // read maximum displecement
  parse("DPHI", Dphi_);

  // secondary structure prediction
  string ss;
  parse("PSIPRED", ss);
  // check length
  if(ss.size()!= getNumberOfArguments()/2) error("Length of prediction should be equal to half the number of arguments in ARG");
  // split ss string into vector of strings - ss_pred_
  for(unsigned i=0; i<ss.size(); ++i) {std::string s(1, ss[i]); ss_pred_.push_back(s);}

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

  // prepare stuff
  setup_restraint();

  // adjust for multiple-time steps
  MCstride_ *= getStride();

  for(unsigned i=0; i<phi_label_.size(); ++i)
    log.printf("  initial value of phi_%s parameter %f\n", phi_label_[i].c_str(), phi_[phi_label_[i]]);
  log.printf("  maximum MC move of the phi parameters %f\n", Dphi_);
  log.printf("  secondary structure prediction %s\n", ss.c_str());
  log.printf("  temperature of the system in energy unit %f\n",kbt_);
  log.printf("  number of MC steps %d\n",MCsteps_);
  log.printf("  do MC every %d steps\n", MCstride_);

  for(unsigned i=0; i<phi_label_.size(); ++i) {
    addComponent("phi_"+phi_label_[i]);    componentIsNotPeriodic("phi_"+phi_label_[i]);
    addComponent("accphi_"+phi_label_[i]); componentIsNotPeriodic("accphi_"+phi_label_[i]);
  }

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
void SSRestraint::setup_restraint()
{
  // set up parameters for the forward models
  // ALPHA-HELIX          BETA SHEETS
  phi_m_["H"] = -1.12;    phi_m_["E"] = -2.33;
  psi_m_["H"] = -0.84;    psi_m_["E"] =  2.26;
  phi_s_["H"] =  0.56;    phi_s_["E"] =  0.68;
  psi_s_["H"] =  0.56;    psi_s_["E"] =  0.58;
  // auxiliary variables, two for each type
  const vector<string> label = {"H","E"};
  for(unsigned i=0; i<label.size(); ++i) {
    string l = label[i];
    fmod_a_[l]= exp(-1.0/phi_s_[l]/phi_s_[l]-1.0/psi_s_[l]/psi_s_[l]);
    fmod_b_[l]= exp(+1.0/phi_s_[l]/phi_s_[l]+1.0/psi_s_[l]/psi_s_[l]) - fmod_a_[l];
  }

  // set prior parameters - mean and sigma prefactor
  // 1) estimate of proportion from benchmark
  phi_mean_["HH"] = 0.859981409;
  phi_mean_["HE"] = 0.019667669;
  phi_mean_["HC"] = 0.073509029;
  phi_mean_["EH"] = 0.007829839;
  phi_mean_["EE"] = 0.742872412;
  phi_mean_["EC"] = 0.062146985;
  // 2) prefactor for sigma calculation
  phi_sig_["HH"] = 1268353.0;
  phi_sig_["HE"] = 759724.0;
  phi_sig_["HC"] = 1220027.0;
  phi_sig_["EH"] = 1268353.0;
  phi_sig_["EE"] = 759724.0;
  phi_sig_["EC"] = 1220027.0;

  // reset acceptances
  for(unsigned i=0; i<phi_label_.size(); ++i) MCaccphi_[phi_label_[i]] = 0;
}

// get Gaussian prior for a single phi
double SSRestraint::getPrior(double p, double pmean, double psig)
{
  // calculate sigma
  double s = sqrt(p*(1.0-p)/psig);
  // calculate prior - excluding 2pi which is constant
  double prior = 0.5 * ( p - pmean ) * ( p - pmean ) / s / s + std::log(s);

  return kbt_ * prior;
}

// calculate priors on all phi
double SSRestraint::getPriors(map <string, double> ps,
                              map <string, double> pms, map <string, double> pss)
{
  double ene = 0.0;
  // cycle on map
  for(unsigned i=0; i<phi_label_.size(); ++i) {
    // map key string
    string s = phi_label_[i];
    // add contribution
    ene += getPrior(ps[s], pms[s], pss[s]);
  }
  return ene;
}

map<string, double> SSRestraint::get_p_coefficients(map <string, double> phi)
{

  // 3x3 coefficient matrix
  map<string, double>  coeff;
  // these are directly sampled
  coeff["HH"] = phi["HH"]; coeff["HE"] = phi["HE"]; coeff["HC"] = phi["HC"];
  coeff["EH"] = phi["EH"]; coeff["EE"] = phi["EE"]; coeff["EC"] = phi["EC"];
  // these are from normalization
  coeff["CH"] = 1.0 - phi["HH"] - phi["EH"];
  coeff["CE"] = 1.0 - phi["HE"] - phi["EE"];
  coeff["CC"] = 1.0 - phi["HC"] - phi["EC"];

  return coeff;
}


// used to update Bayesian nuisance phi parameters
double SSRestraint::getEnergy(map <string, double> phi,
                              const vector<double> &pHl, const vector<double> &pEl, const vector<double> &pCl)
{
  // get p(d|s) coefficients
  map<string, double> coeff = get_p_coefficients(phi);

  // calculate energy
  double ene = 0.0;
  // cycle on arguments
  for(unsigned i=rank_; i<ss_pred_.size(); i=i+nrep_) {
    // get ss type
    string ss = ss_pred_[i];
    // retrieve p(d|s) coefficients
    double a = coeff[ss+"H"];
    double b = coeff[ss+"E"];
    double c = coeff[ss+"C"];
    // calculate likelihood
    double like = ( a * pHl[i] + b * pEl[i] + c * pCl[i] );
    // add to energy
    ene += -kbt_ * std::log(like);
  }

  // sum energy
  comm.Sum(&ene, 1);

  // add priors on phi
  ene += getPriors(phi, phi_mean_, phi_sig_);

  return ene;
}

bool SSRestraint::doAccept(double oldE, double newE) {
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

double SSRestraint::proposeMove(double x, double xmin, double xmax, double dxmax)
{
  double r = static_cast<double>(rand()) / RAND_MAX;
  double dx = -dxmax + r * 2.0 * dxmax;
  double x_new = x + dx;
  // check boundaries
  if(x_new > xmax) {
    double delta = floor((x_new - xmax)/(xmax-xmin));
    double dx = x_new - xmax - delta * (xmax-xmin);
    x_new = xmax - dx;
  }
  if(x_new < xmin) {
    double delta = floor((xmin - x_new)/(xmax-xmin));
    double dx = xmin - x_new - delta * (xmax-xmin);
    x_new = xmin + dx;
  }
  return x_new;
}

void SSRestraint::proposeMoveCouple(double &x1, double &x2, double dx)
{
  if(rand()%2==0) {
    // propose move in x1
    x1 = proposeMove(x1, 0.0, 1.0, dx);
    // propose move in x2
    x2 = proposeMove(x2, 0.0, 1.0-x1, dx);
  } else {
    // propose move in x2
    x2 = proposeMove(x2, 0.0, 1.0, dx);
    // propose move in x1
    x1 = proposeMove(x1, 0.0, 1.0-x2, dx);
  }
}

void SSRestraint::doMonteCarlo(double oldE, long int step,
                               const vector<double> &pHl, const vector<double> &pEl, const vector<double> &pCl)
{
  // phi are moved in pairs to ensure that sum is lower or equal to one
  vector< pair< string, string > > pairs;
  pairs.push_back(make_pair("HH","EH"));
  pairs.push_back(make_pair("HE","EE"));
  pairs.push_back(make_pair("HC","EC"));
  // cycle on MC steps
  for(unsigned i=0; i<MCsteps_; ++i) {
    // cycle on phi couples
    for(unsigned j=0; j<pairs.size(); ++j) {
      // new map phi_ - copying old one
      map <string, double> phi_new(phi_);
      // change phi couple
      proposeMoveCouple(phi_new[pairs[j].first], phi_new[pairs[j].second], Dphi_);
      // calculate new energy
      double newE = getEnergy(phi_new, pHl, pEl, pCl);
      // accept or reject
      bool accept = doAccept(oldE, newE);
      if(accept) {
        // update value of phi for the couple
        phi_[pairs[j].first]  = phi_new[pairs[j].first];
        phi_[pairs[j].second] = phi_new[pairs[j].second];
        // increment acceptance
        MCaccphi_[pairs[j].first]++;
        MCaccphi_[pairs[j].second]++;
        // update energy
        oldE = newE;
      }
    }
  }
// this is needed when restarting simulations
  if(MCfirst_==-1) MCfirst_=step;
// calculate number of trials
  double MCtrials = std::floor(static_cast<double>(step-MCfirst_) / static_cast<double>(MCstride_))+1.0;
// phi acceptances
  for(unsigned i=0; i<phi_label_.size(); ++i) {
    double accphi = static_cast<double>(MCaccphi_[phi_label_[i]]) / static_cast<double>(MCsteps_) / MCtrials;
    getPntrToComponent("accphi_"+phi_label_[i])->set(accphi);
  }
}

void SSRestraint::calculate()
{
  // allocate force vector
  vector<double> force(getNumberOfArguments(), 0.0);

  // allocate auxiliary vectors
  vector<double> pHl(ss_pred_.size(), 0.0);
  vector<double> pEl(ss_pred_.size(), 0.0);
  vector<double> pCl(ss_pred_.size(), 0.0);

  // get p(d|s) coefficients
  map<string, double> coeff = get_p_coefficients(phi_);

  // calculate energy
  double ene = 0.0;
  // cycle on sequence
  for(unsigned i=rank_; i<ss_pred_.size(); i=i+nrep_) {
    // get dihedrals
    double phid = getArgument(i);
    double psid = getArgument(i+ss_pred_.size());
    // calculate auxiliary forward model stuff
    double tmpH = exp(cos(phid-phi_m_["H"])/phi_s_["H"]/phi_s_["H"]+cos(psid-psi_m_["H"])/psi_s_["H"]/psi_s_["H"]);
    double tmpE = exp(cos(phid-phi_m_["E"])/phi_s_["E"]/phi_s_["E"]+cos(psid-psi_m_["E"])/psi_s_["E"]/psi_s_["E"]);
    // calculate forward models
    double pH = (tmpH - fmod_a_["H"]) / fmod_b_["H"];
    double pE = (tmpE - fmod_a_["E"]) / fmod_b_["E"];
    double pC = (1.0 - pH) * (1.0 - pE);
    // and normalization
    double norm = pH + pE + pC;
    // store in lists (with normalization)
    pHl[i] = pH/norm; pEl[i] = pE/norm; pCl[i] = pC/norm;
    // get ss type
    string ss = ss_pred_[i];
    // retrieve p(d|s) coefficients
    double a = coeff[ss+"H"];
    double b = coeff[ss+"E"];
    double c = coeff[ss+"C"];
    // calculate likelihood
    double like = a * pHl[i] + b * pEl[i] + c * pCl[i];
    // add to energy
    ene += -kbt_ * std::log(like);
    // calculate partial derivatives
    double dene_dlike = -kbt_ / like;
    //
    double dlike_dpH = ( a  - ( a * pHl[i] + b * pEl[i] + c * pCl[i] ) ) / norm;
    double dlike_dpE = ( b  - ( a * pHl[i] + b * pEl[i] + c * pCl[i] ) ) / norm;
    double dlike_dpC = ( c  - ( a * pHl[i] + b * pEl[i] + c * pCl[i] ) ) / norm;
    //
    double dpH_dphid   = -tmpH * sin(phid-phi_m_["H"]) / phi_s_["H"] / phi_s_["H"] / fmod_b_["H"];
    double dpH_dpsid   = -tmpH * sin(psid-psi_m_["H"]) / psi_s_["H"] / psi_s_["H"] / fmod_b_["H"];
    double dpE_dphid   = -tmpE * sin(phid-phi_m_["E"]) / phi_s_["E"] / phi_s_["E"] / fmod_b_["E"];
    double dpE_dpsid   = -tmpE * sin(psid-psi_m_["E"]) / psi_s_["E"] / psi_s_["E"] / fmod_b_["E"];
    double dpC_dphid   = -dpH_dphid * (1.0 - pE) - (1.0 - pH) * dpE_dphid;
    double dpC_dpsid   = -dpH_dpsid * (1.0 - pE) - (1.0 - pH) * dpE_dpsid;
    // apply chain rule
    force[i]                 = -dene_dlike * ( dlike_dpH * dpH_dphid + dlike_dpE * dpE_dphid + dlike_dpC * dpC_dphid );
    force[i+ss_pred_.size()] = -dene_dlike * ( dlike_dpH * dpH_dpsid + dlike_dpE * dpE_dpsid + dlike_dpC * dpC_dpsid );
  }

  // sum energy and derivatives
  comm.Sum(&force[0], force.size());
  comm.Sum(&ene, 1);
  // sum auxiliary lists
  comm.Sum(&pHl[0], pHl.size());
  comm.Sum(&pEl[0], pEl.size());
  comm.Sum(&pCl[0], pCl.size());

  // apply forces
  for(unsigned i=0; i<force.size(); ++i) setOutputForce(i, force[i]);

  // add priors on phi
  ene += getPriors(phi_, phi_mean_, phi_sig_);

  // set value of the bias
  setBias(ene);

  // set values of phi
  for(unsigned i=0; i<phi_label_.size(); ++i)
    getPntrToComponent("phi_"+phi_label_[i])->set(phi_[phi_label_[i]]);

  // get time step
  long int step = getStep();
  // do MC stuff at the right time step
  if(step%MCstride_==0&&!getExchangeStep()) doMonteCarlo(ene, step, pHl, pEl, pCl);
}


}
}


