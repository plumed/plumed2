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
  // number of predictions for each secondary structure type
  map <string, double> n_;
  // phi prior parameters
  map <string, double> phi_mean_;
  map <string, double> phi_sig_;
  // and labels
  const vector<string> phi_label_ = {"HH","HE","EE","EC","CC","CH"};
  // reference dihedrals
  map <string, double> phi_ref_;
  map <string, double> psi_ref_;
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
  bool   doAccept(double oldE, double newE);
  map< pair<string,string>, double > get_p_coefficients(map <string, double> phi);

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
  keys.add("compulsory","PHI_EE0","initial value of the phi_EE parameter");
  keys.add("compulsory","PHI_EC0","initial value of the phi_EC parameter");
  keys.add("compulsory","PHI_CC0","initial value of the phi_CC parameter");
  keys.add("compulsory","PHI_CH0","initial value of the phi_CH parameter");
  keys.add("compulsory","DPHI","maximum MC move of the phi parameters");
  keys.add("compulsory","SS","secondary structure prediction");
  keys.add("optional","TEMP","temperature in energy units");
  keys.add("optional","MC_STEPS","number of MC steps");
  keys.add("optional","MC_STRIDE","MC stride");
  componentsAreNotOptional(keys);
  useCustomisableComponents(keys);
  keys.addOutputComponent("phi","default","phi parameters");
  keys.addOutputComponent("accphi","default","MC acceptance phi");
}

SSRestraint::SSRestraint(const ActionOptions&ao):
  PLUMED_BIAS_INIT(ao),
  MCsteps_(1), MCstride_(1), MCfirst_(-1)
{
  // parse phi dihedrals
  vector<Value*> phi_arg;
  parseArgumentList("ARG",1,phi_arg);
  // parse psi dihedrals
  vector<Value*> psi_arg;
  parseArgumentList("ARG",2,psi_arg);
  // check length
  if(phi_arg.size()!=psi_arg.size()) error("The number of arguments in ARG1 and ARG2 should be the same");

  // merge lists into global list
  vector<Value*> arg;
  for(unsigned i=0; i<phi_arg.size(); ++i) arg.push_back(phi_arg[i]);
  for(unsigned i=0; i<psi_arg.size(); ++i) arg.push_back(psi_arg[i]);
  // read initial values of the phi parameters
  for(unsigned i=0; i<phi_label_.size(); ++i) parse("PHI_"+phi_label_[i]+"0", phi_[phi_label_[i]]);

  // secondary structure prediction
  string ss;
  parse("SS", ss);
  // check length
  if(ss.size()!= phi_arg.size()) error("Length of prediction should be equal to the number of arguments in ARG1/ARG2");
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

  // ask for arguments
  requestArguments(arg);

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
  // set up reference dihedrals
  // ALPHA-HELIX            BETA SHEETS
  phi_ref_["H"] = -1.05;    phi_ref_["E"] = -2.36;
  psi_ref_["H"] = -0.79;    psi_ref_["E"] =  2.36;

  // total number of predictions for each secondary structure type
  n_["H"] = 0.0; n_["E"] = 0.0; n_["C"] = 0.0;
  for(unsigned i=0; i<ss_pred_.size(); ++i) n_[ss_pred_[i]] += 1.0;

  // set prior parameters - mean and sigma of truncated normal
  phi_mean_["HH"] = 0.  ; phi_sig_["HH"] = 0.   ;
  phi_mean_["HE"] = 0.  ; phi_sig_["HE"] = 0.   ;
  phi_mean_["EE"] = 0.  ; phi_sig_["EE"] = 0.   ;
  phi_mean_["EC"] = 0.  ; phi_sig_["EC"] = 0.   ;
  phi_mean_["CC"] = 0.  ; phi_sig_["CC"] = 0.   ;
  phi_mean_["CH"] = 0.  ; phi_sig_["CH"] = 0.   ;

  // reset acceptances
  for(unsigned i=0; i<phi_label_.size(); ++i) MCaccphi_[phi_label_[i]] = 0;
}

// get truncated Gaussian prior for a single phi
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

// calculate priors on all phi
double SSRestraint::getPriors(map <string, double> ps,
                              map <string, double> pms, map <string, double> pss)
{
  double ene = 0.0;
  // cycle on map
  for(map <string, double>::iterator it=ps.begin(); it!=ps.end(); ++it) {
    // map key string
    string s = it->first;
    // add contribution
    ene += getPrior(ps[s], pms[s], pss[s]);
  }
  return ene;
}

// used to update Bayesian nuisance phi parameters
double SSRestraint::getEnergy(map <string, double> phi,
                              const vector<double> &pHl, const vector<double> &pEl, const vector<double> &pCl)
{
  // get p(d|s) coefficients
  map< pair<string,string>, double > coeff = get_p_coefficients(phi);

  // calculate energy
  double ene = 0.0;
  // cycle on positive arguments
  for(unsigned i=rank_; i<ss_pred_.size(); i=i+nrep_) {
    // get ss type
    string ss = ss_pred_[i];
    // and normalization
    double norm = pHl[i] + pEl[i] + pCl[i];
    // retrieve p(d|s) coefficients
    double a = coeff[make_pair(ss,"H")];
    double b = coeff[make_pair(ss,"E")];
    double c = coeff[make_pair(ss,"C")];
    // calculate likelihood
    double like = ( a * pHl[i] + b * pEl[i] + c * pCl[i] ) / norm;
    // add to energy
    ene += -kbt_ * std::log(like);
  }

  // sum energy
  comm.Sum(&ene, 1);

  // add priors on phi
  ene += getPriors(phi, phi_mean_, phi_sig_);

  return ene;
}

double SSRestraint::proposeMove(double x, double xmin, double xmax, double dxmax)
{
  double r = static_cast<double>(rand()) / RAND_MAX;
  double dx = -dxmax + r * 2.0 * dxmax;
  double x_new = x + dx;
// check boundaries
  if(x_new > xmax) {x_new = 2.0 * xmax - x_new;}
  if(x_new < xmin) {x_new = 2.0 * xmin - x_new;}
  return x_new;
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

void SSRestraint::doMonteCarlo(double oldE, long int step,
                               const vector<double> &pHl, const vector<double> &pEl, const vector<double> &pCl)
{
// cycle on MC steps
  for(unsigned i=0; i<MCsteps_; ++i) {
    // cycle on phi
    for(unsigned j=0; j<phi_label_.size(); ++j) {
      // new map phi_
      map <string, double> phi_new = phi_;
      // propose move for phi
      phi_new[phi_label_[j]] = proposeMove(phi_[phi_label_[j]], 0.0, 1.0, Dphi_);
      // calculate new energy
      double newE = getEnergy(phi_new, pHl, pEl, pCl);
      // accept or reject
      bool accept = doAccept(oldE, newE);
      if(accept) {
        phi_ = phi_new;
        MCaccphi_[phi_label_[j]]++;
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

map< pair<string,string>, double > SSRestraint::get_p_coefficients(map <string, double> phi)
{
  map< pair<string,string>, double >  coeff;
  // normalization factors
  double norm_H = phi["HH"]*n_["H"]+(1.0-phi["EC"]-phi["EE"])*n_["E"]+phi["CH"]*n_["C"];
  double norm_E = phi["HE"]*n_["H"]+phi["EE"]*n_["E"]+(1.0-phi["CH"]-phi["CC"])*n_["C"];
  double norm_C = (1.0-phi["HH"]-phi["HE"])*n_["H"]+phi["EC"]*n_["E"]+phi["CC"]*n_["C"];

  // 3x3 coefficient matrix
  coeff[make_pair("H","H")] = phi["HH"] * n_["H"] / norm_H;
  coeff[make_pair("H","E")] = phi["HE"] * n_["H"] / norm_E;
  coeff[make_pair("H","C")] = (1.0 - phi["HH"] - phi["HE"]) * n_["H"] / norm_C;
  coeff[make_pair("E","H")] = (1.0 - phi["EE"] - phi["EC"]) * n_["E"] / norm_H;
  coeff[make_pair("E","E")] = phi["EE"] * n_["E"] / norm_E;
  coeff[make_pair("E","C")] = phi["EC"] * n_["E"] / norm_C;
  coeff[make_pair("C","H")] = phi["CH"] * n_["C"] / norm_H;
  coeff[make_pair("C","E")] = (1.0 - phi["CH"] - phi["CC"]) * n_["C"] / norm_E;
  coeff[make_pair("C","C")] = phi["CC"] * n_["C"] / norm_C;

  return coeff;
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
  map< pair<string,string>, double > coeff = get_p_coefficients(phi_);

  // calculate energy
  double ene = 0.0;
  // cycle on sequence
  for(unsigned i=rank_; i<ss_pred_.size(); i=i+nrep_) {
    // get dihedrals
    double phid = getArgument(i);
    double psid = getArgument(i+ss_pred_.size());
    // get ss type
    string ss = ss_pred_[i];
    // calculate forward models
    double pH = 0.25*(1.0+cos(phid-phi_ref_["H"]))*(1.0+cos(psid-psi_ref_["H"]));
    double pE = 0.25*(1.0+cos(phid-phi_ref_["E"]))*(1.0+cos(psid-psi_ref_["E"]));
    double pC = (1.0 - pH) * (1.0 - pE);
    // store in lists
    pHl[i] = pH; pEl[i] = pE; pCl[i] = pC;
    // and normalization
    double norm = pH + pE + pC;
    // retrieve p(d|s) coefficients
    double a = coeff[make_pair(ss,"H")];
    double b = coeff[make_pair(ss,"E")];
    double c = coeff[make_pair(ss,"C")];
    // calculate likelihood
    double like = ( a * pH + b * pE + c * pC ) / norm;
    // add to energy
    ene += -kbt_ * std::log(like);
    // calculate partial derivatives
    double dene_dlike = -kbt_ / like;
    //
    double dlike_dpH = ( a * norm - ( a * pH + b * pE + c * pC ) ) / norm / norm;
    double dlike_dpE = ( b * norm - ( a * pH + b * pE + c * pC ) ) / norm / norm;
    double dlike_dpC = ( c * norm - ( a * pH + b * pE + c * pC ) ) / norm / norm;
    //
    double dpH_dphid   = -0.25 * sin(phid-phi_ref_["H"]) * (1.0+cos(psid-psi_ref_["H"]));
    double dpH_dpsid   = -0.25 * (1.0+cos(phid-phi_ref_["H"])) * sin(psid-psi_ref_["H"]);
    double dpE_dphid   = -0.25 * sin(phid-phi_ref_["E"]) * (1.0+cos(psid-psi_ref_["E"]));
    double dpE_dpsid   = -0.25 * (1.0+cos(phid-phi_ref_["E"])) * sin(psid-psi_ref_["E"]);
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


