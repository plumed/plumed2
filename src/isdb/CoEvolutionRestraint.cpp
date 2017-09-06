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

//+PLUMEDOC ISDB_BIAS COEVOLUTION
/*


*/
//+ENDPLUMEDOC


class CoEvolutionRestraint : public bias::Bias
{
  // slope
  double slope_;
  // alpha parameters
  double alpha_p_;
  double alpha_n_;
  double Dalpha_;
  // alpha prior parameters
  double alpha_p_mean_;
  double alpha_p_sig_;
  double alpha_n_mean_;
  double alpha_n_sig_;
  // sigma parameters
  double sigma_p_;
  double sigma_n_;
  double Dsigma_;
  // "cutoff", gamma, and P0 parameters;
  vector<double> R0_;
  vector<double> gamma_;
  double P0_;
  // temperature in kbt
  double kbt_;
  // Monte Carlo stuff
  int MCsteps_;
  int MCstride_;
  unsigned int MCaccalpha_p_;
  unsigned int MCaccalpha_n_;
  long int MCfirst_;
  unsigned int MCaccsigma_p_;
  unsigned int MCaccsigma_n_;
  // parallel stuff
  unsigned rank_;
  unsigned nrep_;
  // use jeffreys
  bool   doJeffreys_;

  void   setup_restraint(double R0, double gamma, string res_file);
  double getPrior(double a, double amean, double asig);
  double getJeffreysPrior(double a, double amean, double s);
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
  keys.add("compulsory","NPOS","number of positives");
  keys.add("compulsory","NNEG","number of negatives");
  keys.add("compulsory","RES_FILE","file with residue ids for each argument");
  keys.add("optional","SLOPE","add logarithmic slope");
  keys.add("optional","TEMP","temperature in energy units");
  keys.add("optional","MC_STEPS","number of MC steps");
  keys.add("optional","MC_STRIDE","MC stride");
  keys.addFlag("JEFFREYS",false,"use Jeffreys prior for sigma pos and neg");
  keys.add("optional","SIGMA_P0","initial value of the sigma positive parameter");
  keys.add("optional","SIGMA_N0","initial value of the sigma negative parameter");
  keys.add("optional","DSIGMA","maximum MC move of the sigma parameters");
  componentsAreNotOptional(keys);
  useCustomisableComponents(keys);
  keys.addOutputComponent("alphap",   "default","alpha positive parameter");
  keys.addOutputComponent("alphan",   "default","alpha negative parameter");
  keys.addOutputComponent("accalphap","default","MC acceptance alpha positive");
  keys.addOutputComponent("accalphan","default","MC acceptance alpha negative");
  keys.addOutputComponent("sigmap","JEFFREYS","sigma positive parameter");
  keys.addOutputComponent("sigman","JEFFREYS","sigma negative parameter");
  keys.addOutputComponent("accsigmap","JEFFREYS","MC acceptance sigma positive");
  keys.addOutputComponent("accsigman","JEFFREYS","MC acceptance sigma negative");
}

CoEvolutionRestraint::CoEvolutionRestraint(const ActionOptions&ao):
  PLUMED_BIAS_INIT(ao), slope_(0.0),
  sigma_p_(0.0), sigma_n_(0.0), Dsigma_(0.0),
  MCsteps_(1), MCstride_(1), MCaccalpha_p_(0),
  MCaccalpha_n_(0), MCfirst_(-1),
  MCaccsigma_p_(0), MCaccsigma_n_(0),
  doJeffreys_(false)
{
  // additional slope
  parse("SLOPE", slope_);

  // alpha stuff
  parse("ALPHA_P0", alpha_p_);
  parse("ALPHA_N0", alpha_n_);
  parse("DALPHA",   Dalpha_);

  // forward model parameters
  double R0;
  parse("R0", R0);
  if(R0<=0.) error("R0 should be strictly positive");
  double gamma;
  parse("GAMMA", gamma);
  if(gamma<=0.) error("GAMMA should be strictly positive");
  parse("P0", P0_);
  if(P0_<=0. or P0_>1.) error("P0 should be strictly positive and lower than 1");

  string res_file;
  parse("RES_FILE", res_file);

  // set constant prior parameters
  alpha_p_mean_ = 0.2399;
  alpha_p_sig_ = 4269552.0;
  alpha_n_mean_ = 0.0064;
  alpha_n_sig_ = 165805402.0;

  // check if using Jeffreys prior
  parseFlag("JEFFREYS", doJeffreys_);
  // parse initial values of sigma positive and negative
  parse("SIGMA_P0", sigma_p_);
  parse("SIGMA_N0", sigma_n_);
  parse("DSIGMA",   Dsigma_);
  // check for errors
  if(doJeffreys_) {
    if(sigma_p_ <=0.0 || sigma_p_ >=1.0) error("With JEFFREYS, SIGMA_P0 should be specified and strictly between 0 and 1");
    if(sigma_n_ <=0.0 || sigma_n_ >=1.0) error("With JEFFREYS, SIGMA_N0 should be specified and strictly between 0 and 1");
    if(Dsigma_  <=0.0 || Dsigma_  >=1.0) error("With JEFFREYS, DSIGMA should be specified and strictly between 0 and 1");
  }

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

  // setup restraint
  setup_restraint(R0, gamma, res_file);
  if(R0_.size()!=getNumberOfArguments()) error("Check residue pair file");

  // adjust for multiple-time steps
  MCstride_ *= getStride();

  log.printf("  initial value of alpha_p %f\n",alpha_p_);
  log.printf("  initial value of alpha_n %f\n",alpha_n_);
  log.printf("  maximum MC move of the alpha parameter %f\n",Dalpha_);
  if(doJeffreys_) {
    log.printf("  using Jeffreys prior for sigma\n");
    log.printf("  initial value of sigma_p %f\n",sigma_p_);
    log.printf("  initial value of sigma_n %f\n",sigma_n_);
    log.printf("  maximum MC move of the sigma parameters %f\n",Dsigma_);
  }
  log.printf("  forward model parameter R0 %f\n",R0);
  log.printf("  forward model parameter P0 %f\n",P0_);
  log.printf("  forward model parameter gamma %f\n",gamma);
  log.printf("  temperature of the system in energy unit %f\n",kbt_);
  log.printf("  number of MC steps %d\n",MCsteps_);
  log.printf("  do MC every %d steps\n", MCstride_);

  addComponent("alphap");    componentIsNotPeriodic("alphap");
  addComponent("accalphap"); componentIsNotPeriodic("accalphap");
  addComponent("alphan");    componentIsNotPeriodic("alphan");
  addComponent("accalphan"); componentIsNotPeriodic("accalphan");
  if(doJeffreys_) {
    addComponent("sigmap");    componentIsNotPeriodic("sigmap");
    addComponent("accsigmap"); componentIsNotPeriodic("accsigmap");
    addComponent("sigman");    componentIsNotPeriodic("sigman");
    addComponent("accsigman"); componentIsNotPeriodic("accsigman");
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

void  CoEvolutionRestraint::setup_restraint(double R0, double gamma, string res_file)
{
  // reference distance and slope for each pairs of residues
  map< pair<string,string>, pair<double,double> > DREF_;
  // load from Baker's paper - in Angstrom
  DREF_[make_pair("G","G")]=make_pair(4.467,0.017); DREF_[make_pair("G","A")]=make_pair(5.201,0.269); DREF_[make_pair("G","S")]=make_pair(5.510,0.153);
  DREF_[make_pair("G","V")]=make_pair(5.671,0.107); DREF_[make_pair("G","C")]=make_pair(5.777,0.129); DREF_[make_pair("G","T")]=make_pair(5.619,0.120);
  DREF_[make_pair("G","P")]=make_pair(6.140,0.245); DREF_[make_pair("G","D")]=make_pair(6.135,0.193); DREF_[make_pair("G","N")]=make_pair(6.321,0.169);
  DREF_[make_pair("G","I")]=make_pair(6.413,0.179); DREF_[make_pair("G","L")]=make_pair(6.554,0.125); DREF_[make_pair("G","E")]=make_pair(7.036,0.249);
  DREF_[make_pair("G","Q")]=make_pair(7.297,0.216); DREF_[make_pair("G","M")]=make_pair(7.383,0.255); DREF_[make_pair("G","H")]=make_pair(7.472,0.206);
  DREF_[make_pair("G","K")]=make_pair(8.216,0.358); DREF_[make_pair("G","F")]=make_pair(7.966,0.219); DREF_[make_pair("G","Y")]=make_pair(9.098,0.267);
  DREF_[make_pair("G","R")]=make_pair(9.166,0.334); DREF_[make_pair("G","W")]=make_pair(8.966,0.239); DREF_[make_pair("A","A")]=make_pair(5.381,0.262);
  DREF_[make_pair("A","S")]=make_pair(5.829,0.291); DREF_[make_pair("A","V")]=make_pair(5.854,0.312); DREF_[make_pair("A","C")]=make_pair(6.057,0.394);
  DREF_[make_pair("A","T")]=make_pair(5.982,0.378); DREF_[make_pair("A","P")]=make_pair(6.412,0.399); DREF_[make_pair("A","D")]=make_pair(6.388,0.289);
  DREF_[make_pair("A","N")]=make_pair(6.766,0.349); DREF_[make_pair("A","I")]=make_pair(6.587,0.214); DREF_[make_pair("A","L")]=make_pair(6.707,0.250);
  DREF_[make_pair("A","E")]=make_pair(7.124,0.340); DREF_[make_pair("A","Q")]=make_pair(7.583,0.356); DREF_[make_pair("A","M")]=make_pair(7.605,0.394);
  DREF_[make_pair("A","H")]=make_pair(7.591,0.380); DREF_[make_pair("A","K")]=make_pair(8.327,0.550); DREF_[make_pair("A","F")]=make_pair(8.162,0.260);
  DREF_[make_pair("A","Y")]=make_pair(9.121,0.443); DREF_[make_pair("A","R")]=make_pair(9.365,0.485); DREF_[make_pair("A","W")]=make_pair(9.252,0.290);
  DREF_[make_pair("S","S")]=make_pair(6.190,0.292); DREF_[make_pair("S","V")]=make_pair(6.567,0.205); DREF_[make_pair("S","C")]=make_pair(6.590,0.240);
  DREF_[make_pair("S","T")]=make_pair(6.450,0.214); DREF_[make_pair("S","P")]=make_pair(6.937,0.321); DREF_[make_pair("S","D")]=make_pair(6.760,0.323);
  DREF_[make_pair("S","N")]=make_pair(7.081,0.305); DREF_[make_pair("S","I")]=make_pair(7.142,0.342); DREF_[make_pair("S","L")]=make_pair(7.394,0.287);
  DREF_[make_pair("S","E")]=make_pair(7.483,0.446); DREF_[make_pair("S","Q")]=make_pair(7.807,0.408); DREF_[make_pair("S","M")]=make_pair(8.010,0.369);
  DREF_[make_pair("S","H")]=make_pair(8.051,0.435); DREF_[make_pair("S","K")]=make_pair(8.792,0.445); DREF_[make_pair("S","F")]=make_pair(8.694,0.394);
  DREF_[make_pair("S","Y")]=make_pair(9.594,0.467); DREF_[make_pair("S","R")]=make_pair(9.753,0.483); DREF_[make_pair("S","W")]=make_pair(9.770,0.497);
  DREF_[make_pair("V","V")]=make_pair(6.759,0.145); DREF_[make_pair("V","C")]=make_pair(6.941,0.173); DREF_[make_pair("V","T")]=make_pair(6.791,0.138);
  DREF_[make_pair("V","P")]=make_pair(7.063,0.298); DREF_[make_pair("V","D")]=make_pair(6.972,0.287); DREF_[make_pair("V","N")]=make_pair(7.219,0.232);
  DREF_[make_pair("V","I")]=make_pair(7.441,0.242); DREF_[make_pair("V","L")]=make_pair(7.633,0.179); DREF_[make_pair("V","E")]=make_pair(7.404,0.510);
  DREF_[make_pair("V","Q")]=make_pair(8.008,0.359); DREF_[make_pair("V","M")]=make_pair(8.335,0.295); DREF_[make_pair("V","H")]=make_pair(8.179,0.383);
  DREF_[make_pair("V","K")]=make_pair(8.077,0.634); DREF_[make_pair("V","F")]=make_pair(9.057,0.246); DREF_[make_pair("V","Y")]=make_pair(9.442,0.535);
  DREF_[make_pair("V","R")]=make_pair(9.513,0.514); DREF_[make_pair("V","W")]=make_pair(10.021,0.271); DREF_[make_pair("C","C")]=make_pair(6.426,0.178);
  DREF_[make_pair("C","T")]=make_pair(6.801,0.181); DREF_[make_pair("C","P")]=make_pair(7.157,0.259); DREF_[make_pair("C","D")]=make_pair(6.985,0.299);
  DREF_[make_pair("C","N")]=make_pair(7.205,0.240); DREF_[make_pair("C","I")]=make_pair(7.476,0.295); DREF_[make_pair("C","L")]=make_pair(7.685,0.206);
  DREF_[make_pair("C","E")]=make_pair(7.449,0.538); DREF_[make_pair("C","Q")]=make_pair(7.962,0.347); DREF_[make_pair("C","M")]=make_pair(8.265,0.439);
  DREF_[make_pair("C","H")]=make_pair(8.422,0.203); DREF_[make_pair("C","K")]=make_pair(8.494,0.521); DREF_[make_pair("C","F")]=make_pair(9.026,0.286);
  DREF_[make_pair("C","Y")]=make_pair(9.362,0.585); DREF_[make_pair("C","R")]=make_pair(9.460,0.491); DREF_[make_pair("C","W")]=make_pair(9.752,0.417);
  DREF_[make_pair("T","T")]=make_pair(6.676,0.188); DREF_[make_pair("T","P")]=make_pair(7.062,0.320); DREF_[make_pair("T","D")]=make_pair(6.971,0.307);
  DREF_[make_pair("T","N")]=make_pair(7.159,0.262); DREF_[make_pair("T","I")]=make_pair(7.442,0.259); DREF_[make_pair("T","L")]=make_pair(7.642,0.190);
  DREF_[make_pair("T","E")]=make_pair(7.628,0.409); DREF_[make_pair("T","Q")]=make_pair(8.055,0.378); DREF_[make_pair("T","M")]=make_pair(8.397,0.292);
  DREF_[make_pair("T","H")]=make_pair(8.221,0.417); DREF_[make_pair("T","K")]=make_pair(8.715,0.464); DREF_[make_pair("T","F")]=make_pair(9.030,0.264);
  DREF_[make_pair("T","Y")]=make_pair(9.813,0.430); DREF_[make_pair("T","R")]=make_pair(9.764,0.477); DREF_[make_pair("T","W")]=make_pair(9.980,0.315);
  DREF_[make_pair("P","P")]=make_pair(7.288,0.339); DREF_[make_pair("P","D")]=make_pair(7.321,0.416); DREF_[make_pair("P","N")]=make_pair(7.497,0.334);
  DREF_[make_pair("P","I")]=make_pair(7.554,0.336); DREF_[make_pair("P","L")]=make_pair(7.751,0.317); DREF_[make_pair("P","E")]=make_pair(7.938,0.475);
  DREF_[make_pair("P","Q")]=make_pair(8.308,0.410); DREF_[make_pair("P","M")]=make_pair(8.247,0.388); DREF_[make_pair("P","H")]=make_pair(8.537,0.457);
  DREF_[make_pair("P","K")]=make_pair(9.198,0.550); DREF_[make_pair("P","F")]=make_pair(8.895,0.425); DREF_[make_pair("P","Y")]=make_pair(9.965,0.506);
  DREF_[make_pair("P","R")]=make_pair(10.266,0.506); DREF_[make_pair("P","W")]=make_pair(9.719,0.462); DREF_[make_pair("D","D")]=make_pair(8.001,0.392);
  DREF_[make_pair("D","N")]=make_pair(7.672,0.337); DREF_[make_pair("D","I")]=make_pair(7.472,0.341); DREF_[make_pair("D","L")]=make_pair(7.696,0.348);
  DREF_[make_pair("D","E")]=make_pair(8.945,0.354); DREF_[make_pair("D","Q")]=make_pair(8.601,0.357); DREF_[make_pair("D","M")]=make_pair(8.401,0.361);
  DREF_[make_pair("D","H")]=make_pair(8.634,0.325); DREF_[make_pair("D","K")]=make_pair(9.306,0.343); DREF_[make_pair("D","F")]=make_pair(9.111,0.351);
  DREF_[make_pair("D","Y")]=make_pair(9.979,0.676); DREF_[make_pair("D","R")]=make_pair(10.123,0.327); DREF_[make_pair("D","W")]=make_pair(9.867,0.475);
  DREF_[make_pair("N","N")]=make_pair(7.682,0.249); DREF_[make_pair("N","I")]=make_pair(7.631,0.341); DREF_[make_pair("N","L")]=make_pair(7.889,0.279);
  DREF_[make_pair("N","E")]=make_pair(8.485,0.423); DREF_[make_pair("N","Q")]=make_pair(8.502,0.373); DREF_[make_pair("N","M")]=make_pair(8.550,0.310);
  DREF_[make_pair("N","H")]=make_pair(8.672,0.289); DREF_[make_pair("N","K")]=make_pair(9.319,0.398); DREF_[make_pair("N","F")]=make_pair(9.168,0.393);
  DREF_[make_pair("N","Y")]=make_pair(10.039,0.586); DREF_[make_pair("N","R")]=make_pair(10.135,0.372); DREF_[make_pair("N","W")]=make_pair(9.976,0.458);
  DREF_[make_pair("I","I")]=make_pair(8.096,0.321); DREF_[make_pair("I","L")]=make_pair(8.342,0.261); DREF_[make_pair("I","E")]=make_pair(7.949,0.453);
  DREF_[make_pair("I","Q")]=make_pair(8.302,0.406); DREF_[make_pair("I","M")]=make_pair(8.874,0.327); DREF_[make_pair("I","H")]=make_pair(8.523,0.379);
  DREF_[make_pair("I","K")]=make_pair(8.329,0.582); DREF_[make_pair("I","F")]=make_pair(9.602,0.347); DREF_[make_pair("I","Y")]=make_pair(9.719,0.589);
  DREF_[make_pair("I","R")]=make_pair(9.746,0.557); DREF_[make_pair("I","W")]=make_pair(10.470,0.397); DREF_[make_pair("L","L")]=make_pair(8.522,0.198);
  DREF_[make_pair("L","E")]=make_pair(8.077,0.475); DREF_[make_pair("L","Q")]=make_pair(8.480,0.411); DREF_[make_pair("L","M")]=make_pair(9.122,0.318);
  DREF_[make_pair("L","H")]=make_pair(8.676,0.401); DREF_[make_pair("L","K")]=make_pair(8.479,0.591); DREF_[make_pair("L","F")]=make_pair(9.900,0.260);
  DREF_[make_pair("L","Y")]=make_pair(9.889,0.611); DREF_[make_pair("L","R")]=make_pair(9.852,0.578); DREF_[make_pair("L","W")]=make_pair(10.707,0.331);
  DREF_[make_pair("E","E")]=make_pair(9.863,0.389); DREF_[make_pair("E","Q")]=make_pair(9.328,0.450); DREF_[make_pair("E","M")]=make_pair(8.870,0.511);
  DREF_[make_pair("E","H")]=make_pair(9.454,0.443); DREF_[make_pair("E","K")]=make_pair(9.842,0.434); DREF_[make_pair("E","F")]=make_pair(9.403,0.512);
  DREF_[make_pair("E","Y")]=make_pair(10.544,0.469); DREF_[make_pair("E","R")]=make_pair(10.713,0.363); DREF_[make_pair("E","W")]=make_pair(10.303,0.493);
  DREF_[make_pair("Q","Q")]=make_pair(9.074,0.436); DREF_[make_pair("Q","M")]=make_pair(9.102,0.498); DREF_[make_pair("Q","H")]=make_pair(9.391,0.401);
  DREF_[make_pair("Q","K")]=make_pair(9.667,0.521); DREF_[make_pair("Q","F")]=make_pair(9.506,0.451); DREF_[make_pair("Q","Y")]=make_pair(10.534,0.547);
  DREF_[make_pair("Q","R")]=make_pair(10.610,0.535); DREF_[make_pair("Q","W")]=make_pair(10.429,0.490); DREF_[make_pair("M","M")]=make_pair(9.530,0.457);
  DREF_[make_pair("M","H")]=make_pair(9.396,0.342); DREF_[make_pair("M","K")]=make_pair(9.096,0.611); DREF_[make_pair("M","F")]=make_pair(10.253,0.377);
  DREF_[make_pair("M","Y")]=make_pair(10.400,0.661); DREF_[make_pair("M","R")]=make_pair(10.250,0.641); DREF_[make_pair("M","W")]=make_pair(11.110,0.397);
  DREF_[make_pair("H","H")]=make_pair(10.606,0.333); DREF_[make_pair("H","K")]=make_pair(9.582,0.714); DREF_[make_pair("H","F")]=make_pair(9.602,0.542);
  DREF_[make_pair("H","Y")]=make_pair(10.843,0.554); DREF_[make_pair("H","R")]=make_pair(10.879,0.595); DREF_[make_pair("H","W")]=make_pair(10.661,0.458);
  DREF_[make_pair("K","K")]=make_pair(10.662,0.738); DREF_[make_pair("K","F")]=make_pair(9.344,0.441); DREF_[make_pair("K","Y")]=make_pair(10.627,0.704);
  DREF_[make_pair("K","R")]=make_pair(11.322,0.648); DREF_[make_pair("K","W")]=make_pair(10.136,0.470); DREF_[make_pair("F","F")]=make_pair(10.903,0.460);
  DREF_[make_pair("F","Y")]=make_pair(10.999,0.767); DREF_[make_pair("F","R")]=make_pair(10.577,0.738); DREF_[make_pair("F","W")]=make_pair(11.758,0.447);
  DREF_[make_pair("Y","Y")]=make_pair(11.536,0.855); DREF_[make_pair("Y","R")]=make_pair(11.615,0.822); DREF_[make_pair("Y","W")]=make_pair(11.807,0.684);
  DREF_[make_pair("R","R")]=make_pair(12.050,0.704); DREF_[make_pair("R","W")]=make_pair(11.355,0.889); DREF_[make_pair("W","W")]=make_pair(12.806,0.473);

// open residue pair file
  ifstream rfile;
  rfile.open(res_file);
  string line, res0, res1;
// iterator for residue map
  map< pair<string,string>, pair<double,double> >::iterator it;
// read file
  if (rfile.is_open()) {
    // read line by line
    while ( getline (rfile, line) )
    {
      // split line into strings separated by a space
      stringstream ss(line);
      // read residue pair
      ss >> res0; ss >> res1;
      // find entry in DREF_
      pair<string,string> key1 = make_pair(res0, res1);
      pair<string,string> key2 = make_pair(res1, res0);
      // look for key1 in DREF_
      it = DREF_.find(key1);
      // and add to reference distances and slopes - convert to nm, multiply by prefactor
      if (it != DREF_.end()) {
        R0_.push_back( 0.1 * DREF_[key1].first * R0 );
        gamma_.push_back( 0.1 * DREF_[key1].second * gamma );
      } else {
        R0_.push_back( 0.1 * DREF_[key2].first * R0 );
        gamma_.push_back( 0.1 * DREF_[key2].second * gamma );
      }
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
  // calculate Gaussian prior - excluding 2pi which is constant
  double prior = 0.5 * ( a - amean ) * ( a - amean ) / s / s + std::log(s);

  return kbt_ * prior;
}

// calculate Gaussian + Jeffreys prior for a single alpha
double CoEvolutionRestraint::getJeffreysPrior(double a, double amean, double s)
{
  // calculate Gaussian prior + Jeffreys - excluding 2pi which is constant
  double prior = 0.5 * ( a - amean ) * ( a - amean ) / s / s + 2.0 * std::log(s);

  return kbt_ * prior;
}



// used to update Bayesian nuisance parameters
double CoEvolutionRestraint::getEnergy
(double alpha_p, double alpha_n, const vector<double> &fmod)
{
  // calculate energy
  double ene = 0.0;
  // cycle on arguments (= positive data)
  for(unsigned i=rank_; i<getNumberOfArguments(); i=i+nrep_) {
    // calculate data likelihood
    double like = alpha_p * fmod[i] + alpha_n * ( 1.0 - fmod[i] );
    // add to energy
    ene += -kbt_ * std::log(like);
  }

  // sum energy
  comm.Sum(&ene, 1);

  if(doJeffreys_ == false) {
    // add Gaussian prior on alpha_p
    ene += getPrior(alpha_p, alpha_p_mean_, alpha_p_sig_);
    // add Gaussian prior on alpha_n
    ene += getPrior(alpha_n, alpha_n_mean_, alpha_n_sig_);
  } else {
    // add Gaussian + Jeffreys prior on alpha_p
    ene += getJeffreysPrior(alpha_p, alpha_p_mean_, sigma_p_);
    // add Gaussian + Jeffreys prior on alpha_n
    ene += getJeffreysPrior(alpha_n, alpha_n_mean_, sigma_n_);
  }

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
// if Jeffreys prior, update sigma_p_ and sigma_n_
  if(doJeffreys_) {
    // cycle on MC steps
    for(unsigned i=0; i<MCsteps_; ++i) {
      // first store old energy, only the sigma_p_ part
      oldE = getJeffreysPrior(alpha_p_, alpha_p_mean_, sigma_p_);
      // 3) propose move in sigma_p_
      double new_sigma_p = proposeMove(sigma_p_, 0.0, 1.0, Dsigma_);
      // calculate new energy
      double newE = getJeffreysPrior(alpha_p_, alpha_p_mean_, new_sigma_p);
      // accept or reject
      bool accept = doAccept(oldE, newE);
      if(accept) {
        sigma_p_ = new_sigma_p;
        MCaccsigma_p_++;
      }
      // first store old energy, only the sigma_n_ part
      oldE = getJeffreysPrior(alpha_n_, alpha_n_mean_, sigma_n_);
      // 4) propose move in sigma_n_
      double new_sigma_n = proposeMove(sigma_n_, 0.0, 1.0, Dsigma_);
      // calculate new energy
      newE = getJeffreysPrior(alpha_n_, alpha_n_mean_, new_sigma_n);
      // accept or reject
      accept = doAccept(oldE, newE);
      if(accept) {
        sigma_n_ = new_sigma_n;
        MCaccsigma_n_++;
      }
    } // end on MC cycle
  } // endif Jeffreys

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
  // if Jeffreys
  if(doJeffreys_) {
    // sigma_p acceptance
    double accsigma_p = static_cast<double>(MCaccsigma_p_) / static_cast<double>(MCsteps_) / MCtrials;
    getPntrToComponent("accsigmap")->set(accsigma_p);
    // sigma_n acceptance
    double accsigma_n = static_cast<double>(MCaccsigma_n_) / static_cast<double>(MCsteps_) / MCtrials;
    getPntrToComponent("accsigman")->set(accsigma_n);
  }
}

void CoEvolutionRestraint::calculate()
{
  // allocate force vector
  vector<double> force(getNumberOfArguments(), 0.0);
  // and forward model vector
  vector<double> fmod(getNumberOfArguments(), 0.0);

  // calculate energy
  double ene = 0.0;
  // cycle on arguments (= positive data)
  for(unsigned i=rank_; i<getNumberOfArguments(); i=i+nrep_) {
    // get distance
    double dist = getArgument(i);
    // calculate forward model
    double p = P0_;
    if(dist >= R0_[i]) p = P0_ * exp(-0.5*(dist-R0_[i])*(dist-R0_[i])/gamma_[i]/gamma_[i]);
    // store forward model
    fmod[i] = p;
    // calculate data likelihood
    double like = alpha_p_ * p + alpha_n_ * ( 1.0 - p );
    // add to energy
    ene += -kbt_ * std::log(like) + kbt_ * slope_ * std::log(dist);
    // calculate force
    force[i] = - kbt_ * slope_ / dist;
    // second contribution
    if(dist >= R0_[i]) {
      double dene_dlike = -kbt_ / like;
      double dlike_dp   = alpha_p_ - alpha_n_;
      double dp_ddist   = -p * (dist-R0_[i]) / gamma_[i] / gamma_[i];
      // apply chain rule
      force[i] += -dene_dlike * dlike_dp * dp_ddist;
    }
  }

  // sum energy, fmod, and derivatives
  comm.Sum(&force[0], force.size());
  comm.Sum(&fmod[0], fmod.size());
  comm.Sum(&ene, 1);

  // apply forces
  for(unsigned i=0; i<force.size(); ++i) setOutputForce(i, force[i]);

  if(doJeffreys_ == false) {
    // add Gaussian prior on alpha_p
    ene += getPrior(alpha_p_, alpha_p_mean_, alpha_p_sig_);
    // add Gaussian prior on alpha_n
    ene += getPrior(alpha_n_, alpha_n_mean_, alpha_n_sig_);
  } else {
    // add Gaussian + Jeffreys prior on alpha_p
    ene += getJeffreysPrior(alpha_p_, alpha_p_mean_, sigma_p_);
    // add Gaussian + Jeffreys prior on alpha_n
    ene += getJeffreysPrior(alpha_n_, alpha_n_mean_, sigma_n_);
  }

  // set value of the bias
  setBias(ene);
  // set values of alpha_p
  getPntrToComponent("alphap")->set(alpha_p_);
  // set values of alpha_n
  getPntrToComponent("alphan")->set(alpha_n_);
  // if Jeffreys
  if(doJeffreys_) {
    // set values of sigma_p
    getPntrToComponent("sigmap")->set(sigma_p_);
    // set values of sigma_n
    getPntrToComponent("sigman")->set(sigma_n_);
  }

  // get time step
  long int step = getStep();
  // do MC stuff at the right time step
  if(step%MCstride_==0&&!getExchangeStep()) doMonteCarlo(ene, fmod, step);

}


}
}


