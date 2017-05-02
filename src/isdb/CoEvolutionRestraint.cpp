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
#include <fstream>
#include <sstream>
#include <map>


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
  vector<double> R0_;
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

  void   setup_restraint(double R0, string res_file);
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
  keys.add("compulsory","R0","forward model parameter R0");
  keys.add("compulsory","P0","forward model parameter P0");
  keys.add("compulsory","GAMMA","forward model parameter gamma");
  keys.add("compulsory","ALPHA_P_MEAN","alpha_p prior mean");
  keys.add("compulsory","ALPHA_P_SIG","alpha_p prior sigma constant");
  keys.add("compulsory","ALPHA_N_MEAN","alpha_n prior mean");
  keys.add("compulsory","ALPHA_N_SIG","alpha_n prior sigma constant");
  keys.add("compulsory","NPOS","number of positives");
  keys.add("compulsory","NNEG","number of negatives");
  keys.add("compulsory","RES_FILE","file with residue ids for each argument");
  keys.add("optional","TEMP","temperature in energy units");
  keys.add("optional","MC_STEPS","number of MC steps");
  keys.add("optional","MC_STRIDE","MC stride");
  componentsAreNotOptional(keys);
  useCustomisableComponents(keys);
  keys.addOutputComponent("alphap",   "default","alpha positive parameter");
  keys.addOutputComponent("alphan",   "default","alpha negative parameter");
  keys.addOutputComponent("accalphap","default","MC acceptance alpha positive");
  keys.addOutputComponent("accalphan","default","MC acceptance alpha negative");
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

  // forward model parameters
  double R0;
  parse("R0", R0);
  if(R0<=0.) error("R0 should be strictly positive");
  parse("P0", P0_);
  if(P0_<=0. or P0_>1.) error("P0 should be strictly positive and lower than 1");
  parse("GAMMA", gamma_);
  if(gamma_<=0.) error("GAMMA should be strictly positive");

  string res_file;
  parse("RES_FILE", res_file);

  // prior parameters
  parse("ALPHA_P_MEAN", alpha_p_mean_);
  parse("ALPHA_P_SIG",  alpha_p_sig_);
  parse("ALPHA_N_MEAN", alpha_n_mean_);
  parse("ALPHA_N_SIG",  alpha_n_sig_);

  // number of positives and negatives
  parse("NPOS", npos_);
  if(npos_<0) error("NPOS should be positive");
  parse("NNEG", nneg_);
  if(nneg_<0) error("NNEG should be positive");
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

  // setup restraint
  setup_restraint(R0, res_file);
  if(R0_.size()!=getNumberOfArguments()) error("Check residue pair file");

  // adjust for multiple-time steps
  MCstride_ *= getStride();

  log.printf("  initial value of alpha_p %f\n",alpha_p_);
  log.printf("  initial value of alpha_n %f\n",alpha_n_);
  log.printf("  maximum MC move of the alpha parameter %f\n",Dalpha_);
  log.printf("  forward model parameter R0 %f\n",R0);
  log.printf("  forward model parameter P0 %f\n",P0_);
  log.printf("  forward model parameter gamma %f\n",gamma_);
  log.printf("  alpha_p prior mean %f\n",alpha_p_mean_);
  log.printf("  alpha_p prior sigma constant %f\n",alpha_p_sig_);
  log.printf("  alpha_n prior mean %f\n",alpha_n_mean_);
  log.printf("  alpha_n prior sigma constant %f\n",alpha_n_sig_);
  log.printf("  number of positive data points %d\n",npos_);
  log.printf("  number of negative data points %d\n",nneg_);
  log.printf("  temperature of the system in energy unit %f\n",kbt_);
  log.printf("  number of MC steps %d\n",MCsteps_);
  log.printf("  do MC every %d steps\n", MCstride_);

  addComponent("alphap");    componentIsNotPeriodic("alphap");
  addComponent("accalphap"); componentIsNotPeriodic("accalphap");
  addComponent("alphan");    componentIsNotPeriodic("alphan");
  addComponent("accalphan"); componentIsNotPeriodic("accalphan");

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

void  CoEvolutionRestraint::setup_restraint(double R0, string res_file)
{
  // reference distance for each pairs of residues
  map< pair<string,string>, double > DREF_;
  // load from Baker's paper
  DREF_[make_pair("G","G")]=4.467;    DREF_[make_pair("G","A")]=5.201;    DREF_[make_pair("G","S")]=5.510;    DREF_[make_pair("G","V")]=5.671;    DREF_[make_pair("G","C")]=5.777;
  DREF_[make_pair("G","T")]=5.619;    DREF_[make_pair("G","P")]=6.140;    DREF_[make_pair("G","D")]=6.135;    DREF_[make_pair("G","N")]=6.321;    DREF_[make_pair("G","I")]=6.413;
  DREF_[make_pair("G","L")]=6.554;    DREF_[make_pair("G","E")]=7.036;    DREF_[make_pair("G","Q")]=7.297;    DREF_[make_pair("G","M")]=7.383;    DREF_[make_pair("G","H")]=7.472;
  DREF_[make_pair("G","K")]=8.216;    DREF_[make_pair("G","F")]=7.966;    DREF_[make_pair("G","Y")]=9.098;    DREF_[make_pair("G","R")]=9.166;    DREF_[make_pair("G","W")]=8.966;
  DREF_[make_pair("A","A")]=5.381;    DREF_[make_pair("A","S")]=5.829;    DREF_[make_pair("A","V")]=5.854;    DREF_[make_pair("A","C")]=6.057;    DREF_[make_pair("A","T")]=5.982;
  DREF_[make_pair("A","P")]=6.412;    DREF_[make_pair("A","D")]=6.388;    DREF_[make_pair("A","N")]=6.766;    DREF_[make_pair("A","I")]=6.587;    DREF_[make_pair("A","L")]=6.707;
  DREF_[make_pair("A","E")]=7.124;    DREF_[make_pair("A","Q")]=7.583;    DREF_[make_pair("A","M")]=7.605;    DREF_[make_pair("A","H")]=7.591;    DREF_[make_pair("A","K")]=8.327;
  DREF_[make_pair("A","F")]=8.162;    DREF_[make_pair("A","Y")]=9.121;    DREF_[make_pair("A","R")]=9.365;    DREF_[make_pair("A","W")]=9.252;    DREF_[make_pair("S","S")]=6.190;
  DREF_[make_pair("S","V")]=6.567;    DREF_[make_pair("S","C")]=6.590;    DREF_[make_pair("S","T")]=6.450;    DREF_[make_pair("S","P")]=6.937;    DREF_[make_pair("S","D")]=6.760;
  DREF_[make_pair("S","N")]=7.081;    DREF_[make_pair("S","I")]=7.142;    DREF_[make_pair("S","L")]=7.394;    DREF_[make_pair("S","E")]=7.483;    DREF_[make_pair("S","Q")]=7.807;
  DREF_[make_pair("S","M")]=8.010;    DREF_[make_pair("S","H")]=8.051;    DREF_[make_pair("S","K")]=8.792;    DREF_[make_pair("S","F")]=8.694;    DREF_[make_pair("S","Y")]=9.594;
  DREF_[make_pair("S","R")]=9.753;    DREF_[make_pair("S","W")]=9.770;    DREF_[make_pair("V","V")]=6.759;    DREF_[make_pair("V","C")]=6.941;    DREF_[make_pair("V","T")]=6.791;
  DREF_[make_pair("V","P")]=7.063;    DREF_[make_pair("V","D")]=6.972;    DREF_[make_pair("V","N")]=7.219;    DREF_[make_pair("V","I")]=7.441;    DREF_[make_pair("V","L")]=7.633;
  DREF_[make_pair("V","E")]=7.404;    DREF_[make_pair("V","Q")]=8.008;    DREF_[make_pair("V","M")]=8.335;    DREF_[make_pair("V","H")]=8.179;    DREF_[make_pair("V","K")]=8.077;
  DREF_[make_pair("V","F")]=9.057;    DREF_[make_pair("V","Y")]=9.442;    DREF_[make_pair("V","R")]=9.513;    DREF_[make_pair("V","W")]=10.021;   DREF_[make_pair("C","C")]=6.426;
  DREF_[make_pair("C","T")]=6.801;    DREF_[make_pair("C","P")]=7.157;    DREF_[make_pair("C","D")]=6.985;    DREF_[make_pair("C","N")]=7.205;    DREF_[make_pair("C","I")]=7.476;
  DREF_[make_pair("C","L")]=7.685;    DREF_[make_pair("C","E")]=7.449;    DREF_[make_pair("C","Q")]=7.962;    DREF_[make_pair("C","M")]=8.265;    DREF_[make_pair("C","H")]=8.422;
  DREF_[make_pair("C","K")]=8.494;    DREF_[make_pair("C","F")]=9.026;    DREF_[make_pair("C","Y")]=9.362;    DREF_[make_pair("C","R")]=9.460;    DREF_[make_pair("C","W")]=9.752;
  DREF_[make_pair("T","T")]=6.676;    DREF_[make_pair("T","P")]=7.062;    DREF_[make_pair("T","D")]=6.971;    DREF_[make_pair("T","N")]=7.159;    DREF_[make_pair("T","I")]=7.442;
  DREF_[make_pair("T","L")]=7.642;    DREF_[make_pair("T","E")]=7.628;    DREF_[make_pair("T","Q")]=8.055;    DREF_[make_pair("T","M")]=8.397;    DREF_[make_pair("T","H")]=8.221;
  DREF_[make_pair("T","K")]=8.715;    DREF_[make_pair("T","F")]=9.030;    DREF_[make_pair("T","Y")]=9.813;    DREF_[make_pair("T","R")]=9.764;    DREF_[make_pair("T","W")]=9.980;
  DREF_[make_pair("P","P")]=7.288;    DREF_[make_pair("P","D")]=7.321;    DREF_[make_pair("P","N")]=7.497;    DREF_[make_pair("P","I")]=7.554;    DREF_[make_pair("P","L")]=7.751;
  DREF_[make_pair("P","E")]=7.938;    DREF_[make_pair("P","Q")]=8.308;    DREF_[make_pair("P","M")]=8.247;    DREF_[make_pair("P","H")]=8.537;    DREF_[make_pair("P","K")]=9.198;
  DREF_[make_pair("P","F")]=8.895;    DREF_[make_pair("P","Y")]=9.965;    DREF_[make_pair("P","R")]=10.266;   DREF_[make_pair("P","W")]=9.719;    DREF_[make_pair("D","D")]=8.001;
  DREF_[make_pair("D","N")]=7.672;    DREF_[make_pair("D","I")]=7.472;    DREF_[make_pair("D","L")]=7.696;    DREF_[make_pair("D","E")]=8.945;    DREF_[make_pair("D","Q")]=8.601;
  DREF_[make_pair("D","M")]=8.401;    DREF_[make_pair("D","H")]=8.634;    DREF_[make_pair("D","K")]=9.306;    DREF_[make_pair("D","F")]=9.111;    DREF_[make_pair("D","Y")]=9.979;
  DREF_[make_pair("D","R")]=10.123;   DREF_[make_pair("D","W")]=9.867;    DREF_[make_pair("N","N")]=7.682;    DREF_[make_pair("N","I")]=7.631;    DREF_[make_pair("N","L")]=7.889;
  DREF_[make_pair("N","E")]=8.485;    DREF_[make_pair("N","Q")]=8.502;    DREF_[make_pair("N","M")]=8.550;    DREF_[make_pair("N","H")]=8.672;    DREF_[make_pair("N","K")]=9.319;
  DREF_[make_pair("N","F")]=9.168;    DREF_[make_pair("N","Y")]=10.039;   DREF_[make_pair("N","R")]=10.135;   DREF_[make_pair("N","W")]=9.976;    DREF_[make_pair("I","I")]=8.096;
  DREF_[make_pair("I","L")]=8.342;    DREF_[make_pair("I","E")]=7.949;    DREF_[make_pair("I","Q")]=8.302;    DREF_[make_pair("I","M")]=8.874;    DREF_[make_pair("I","H")]=8.523;
  DREF_[make_pair("I","K")]=8.329;    DREF_[make_pair("I","F")]=9.602;    DREF_[make_pair("I","Y")]=9.719;    DREF_[make_pair("I","R")]=9.746;    DREF_[make_pair("I","W")]=10.470;
  DREF_[make_pair("L","L")]=8.522;    DREF_[make_pair("L","E")]=8.077;    DREF_[make_pair("L","Q")]=8.480;    DREF_[make_pair("L","M")]=9.122;    DREF_[make_pair("L","H")]=8.676;
  DREF_[make_pair("L","K")]=8.479;    DREF_[make_pair("L","F")]=9.900;    DREF_[make_pair("L","Y")]=9.889;    DREF_[make_pair("L","R")]=9.852;    DREF_[make_pair("L","W")]=10.707;
  DREF_[make_pair("E","E")]=9.863;    DREF_[make_pair("E","Q")]=9.328;    DREF_[make_pair("E","M")]=8.870;    DREF_[make_pair("E","H")]=9.454;    DREF_[make_pair("E","K")]=9.842;
  DREF_[make_pair("E","F")]=9.403;    DREF_[make_pair("E","Y")]=10.544;   DREF_[make_pair("E","R")]=10.713;   DREF_[make_pair("E","W")]=10.303;   DREF_[make_pair("Q","Q")]=9.074;
  DREF_[make_pair("Q","M")]=9.102;    DREF_[make_pair("Q","H")]=9.391;    DREF_[make_pair("Q","K")]=9.667;    DREF_[make_pair("Q","F")]=9.506;    DREF_[make_pair("Q","Y")]=10.534;
  DREF_[make_pair("Q","R")]=10.610;   DREF_[make_pair("Q","W")]=10.429;   DREF_[make_pair("M","M")]=9.530;    DREF_[make_pair("M","H")]=9.396;    DREF_[make_pair("M","K")]=9.096;
  DREF_[make_pair("M","F")]=10.253;   DREF_[make_pair("M","Y")]=10.400;   DREF_[make_pair("M","R")]=10.250;   DREF_[make_pair("M","W")]=11.110;   DREF_[make_pair("H","H")]=10.606;
  DREF_[make_pair("H","K")]=9.582;    DREF_[make_pair("H","F")]=9.602;    DREF_[make_pair("H","Y")]=10.843;   DREF_[make_pair("H","R")]=10.879;   DREF_[make_pair("H","W")]=10.661;
  DREF_[make_pair("K","K")]=10.662;   DREF_[make_pair("K","F")]=9.344;    DREF_[make_pair("K","Y")]=10.627;   DREF_[make_pair("K","R")]=11.322;   DREF_[make_pair("K","W")]=10.136;
  DREF_[make_pair("F","F")]=10.903;   DREF_[make_pair("F","Y")]=10.999;   DREF_[make_pair("F","R")]=10.577;   DREF_[make_pair("F","W")]=11.758;   DREF_[make_pair("Y","Y")]=11.536;
  DREF_[make_pair("Y","R")]=11.615;   DREF_[make_pair("Y","W")]=11.807;   DREF_[make_pair("R","R")]=12.050;   DREF_[make_pair("R","W")]=11.355;   DREF_[make_pair("W","W")]=12.806;

// open residue pair file
  ifstream rfile;
  rfile.open(res_file);
  string line, res0, res1;
// iterator for residue map
  map< pair<string,string>, double >::iterator it;
// read file
  if (rfile.is_open()) {
    // read line by line
    while ( getline (rfile, line) )
    {
      // split line into strings separated by a space
      stringstream ss(line);
      string buf; // Have a buffer string
      // read residue pair
      ss >> res0; ss >> res1;
      // find entry in DREF_
      pair<string,string> key1 = make_pair(res0, res1);
      pair<string,string> key2 = make_pair(res1, res0);
      // look for key1 in DREF_
      it = DREF_.find(key1);
      // and add to reference distances - convert to nm, multiply by R0
      if (it != DREF_.end()) R0_.push_back( 0.1 * DREF_[key1] * R0 );
      else                   R0_.push_back( 0.1 * DREF_[key2] * R0 );
    }
    rfile.close();
  }
  else error("Unable to open residue file");
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
  getPntrToComponent("accalphap")->set(accalpha_p);
  // alpha_n acceptance
  double accalpha_n = static_cast<double>(MCaccalpha_n_) / static_cast<double>(MCsteps_) / MCtrials;
  getPntrToComponent("accalphan")->set(accalpha_n);
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
    double tmp = exp(-gamma_*(dist-R0_[i]));
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
    double tmp = exp(-gamma_*(dist-R0_[i+npos_]));
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
  getPntrToComponent("alphap")->set(alpha_p_);
  // set values of alpha_n
  getPntrToComponent("alphan")->set(alpha_n_);

  // get time step
  long int step = getStep();
  // do MC stuff at the right time step
  if(step%MCstride_==0&&!getExchangeStep()) doMonteCarlo(ene, fmod, step);

}


}
}


