/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

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
#include "colvar/Colvar.h"
#include "colvar/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Matrix.h"
#include "core/SetupMolInfo.h"
#include "core/ActionSet.h"
#include "tools/File.h"

#include <string>
#include <cmath>
#include <map>
#include <numeric>
#include <ctime>

using namespace std;

namespace PLMD {
namespace isdb {

//+PLUMEDOC ISDB_COLVAR EMMI
/*
Calculate the fit of a structure or ensemble of structures with a cryo-EM density map.

This action implements the multi-scale Bayesian approach to cryo-EM data fitting introduced in  Ref. \cite Hanot113951 .
This method allows efficient and accurate structural modeling of cryo-electron microscopy density maps at multiple scales, from coarse-grained to atomistic resolution, by addressing the presence of random and systematic errors in the data, sample heterogeneity, data correlation, and noise correlation.

The experimental density map is fit by a Gaussian Mixture Model (GMM), which is provided as an external file specified by the keyword
GMM_FILE. We are currently working on a web server to perform
this operation. In the meantime, the user can request a stand-alone version of the GMM code at massimiliano.bonomi_AT_gmail.com.

When run in single-replica mode, this action allows atomistic, flexible refinement of an individual structure into a density map.
Combined with a multi-replica framework (such as the -multi option in GROMACS), the user can model an esemble of structures using
the Metainference approach \cite Bonomi:2016ip .

\warning
To use \ref EMMI, the user should always add a \ref MOLINFO line and specify a pdb file of the system.

\note
To enhance sampling in single-structure refinement, one can use a Replica Exchange Method, such as Parallel Tempering.
In this case, the user should add the NO_AVER flag to the input line.

\note
\ref EMMI can be used in combination with periodic and non-periodic systems. In the latter case, one should
add the NOPBC flag to the input line

\par Examples

In this example, we perform a single-structure refinement based on an experimental cryo-EM map. The map is fit with a GMM, whose
parameters are listed in the file GMM_fit.dat. This file contains one line per GMM component in the following format:

\plumedfile
#! FIELDS Id Weight Mean_0 Mean_1 Mean_2 Cov_00 Cov_01 Cov_02 Cov_11 Cov_12 Cov_22 Beta
     0  2.9993805e+01   6.54628 10.37820 -0.92988  2.078920e-02 1.216254e-03 5.990827e-04 2.556246e-02 8.411835e-03 2.486254e-02  1
     1  2.3468312e+01   6.56095 10.34790 -0.87808  1.879859e-02 6.636049e-03 3.682865e-04 3.194490e-02 1.750524e-03 3.017100e-02  1
     ...
\endplumedfile

To accelerate the computation of the Bayesian score, one can:
- use neighbor lists, specified by the keywords NL_CUTOFF and NL_STRIDE;
- calculate the restraint every other step (or more).

All the heavy atoms of the system are used to calculate the density map. This list can conveniently be provided
using a GROMACS index file.

The input file looks as follows:

\plumedfile
# include pdb info
MOLINFO STRUCTURE=prot.pdb

#  all heavy atoms
protein-h: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein-H

# create EMMI score
gmm: EMMI NOPBC SIGMA_MEAN=0.01 TEMP=300.0 NL_STRIDE=100 NL_CUTOFF=0.01 GMM_FILE=GMM_fit.dat ATOMS=protein-h

# translate into bias - apply every 2 steps
emr: BIASVALUE ARG=gmm.scoreb STRIDE=2

PRINT ARG=emr.* FILE=COLVAR STRIDE=500 FMT=%20.10f
\endplumedfile


*/
//+ENDPLUMEDOC

class EMMI : public Colvar {

private:

// temperature in kbt
  double kbt_;
// model GMM - atom types
  vector<unsigned> GMM_m_type_;
// model GMM - list of atom sigmas - one per atom type
  vector<double> GMM_m_s_;
// model GMM - list of atom weights - one per atom type
  vector<double> GMM_m_w_;
// data GMM - means, weights, and covariances + beta option
  vector<Vector>             GMM_d_m_;
  vector<double>             GMM_d_w_;
  vector< VectorGeneric<6> > GMM_d_cov_;
  vector<int>                GMM_d_beta_;
// overlaps
  vector<double> ovmd_;
  vector<double> ovdd_;
  vector<double> ovmd_ave_;
  vector<double> ov_cut_;
  vector<double> ovdd_cut_;
// and derivatives
  vector<Vector> ovmd_der_;
  vector<Vector> atom_der_;
  vector<double> err_f_;
  vector<double> exp_f_;
// constant quantities;
  double cfact_;
  double inv_sqrt2_, sqrt2_pi_;
// metainference
  unsigned nrep_;
  unsigned replica_;
  vector<double> sigma_;
  vector<double> sigma_mean_;
  vector<double> sigma_min_, sigma_max_;
  vector<double> dsigma_;
  vector<double> sigma0_;
// list of prefactors for overlap between two components of model and data GMM
// pre_fact = 1.0 / (2pi)**1.5 / sqrt(det_md) * Wm * Wd
  vector< double > pre_fact_;
// inverse of the sum of model and data covariances matrices
  vector< VectorGeneric<6> > inv_cov_md_;
// neighbor list
  double   nl_cutoff_;
  unsigned nl_stride_;
  bool first_time_;
  bool no_aver_;
  vector < unsigned > nl_;
// parallel stuff
  unsigned size_;
  unsigned rank_;
// analysis mode
  bool analysis_;
  OFile Devfile_;
  double nframe_;
// pbc
  bool pbc_;
// Monte Carlo stuff
  int MCstride_;
  long int MCfirst_;
  unsigned int MCaccept_;
  // status stuff
  unsigned int statusstride_;
  string       statusfilename_;
  OFile        statusfile_;
  bool         first_status_;
  // sampling or marginal?
  bool do_sampling_;
  // prior exponent
  double prior_;

// read and write status
  void read_status();
  void print_status(long int step);
// propose move
  double proposeMove(double x, double xmin, double xmax, double dxmax);
// accept or reject
  bool doAccept(double oldE, double newE);
// do MonteCarlo
  void doMonteCarlo(long int step);
// read error file
  vector<double> read_exp_errors(string errfile);
// calculate model GMM weights and covariances
  vector<double> get_GMM_m(vector<AtomNumber> &atoms);
// read data GMM file
  void get_GMM_d(string gmm_file);
// check GMM data
  void check_GMM_d(VectorGeneric<6> &cov, double w);
// auxiliary method
  void calculate_useful_stuff();
// get cutoff in overlap
  void get_cutoff_ov();
// get pref_fact and inv_cov_md
  double get_prefactor_inverse (const VectorGeneric<6> &GMM_cov_0, const VectorGeneric<6> &GMM_cov_1,
                                double &GMM_w_0, double &GMM_w_1,
                                VectorGeneric<6> &sum, VectorGeneric<6> &inv_sum);
// calculate self overlaps between data GMM components - ovdd_
  double get_self_overlap(unsigned id);
// calculate overlap between two components
  double get_overlap(const Vector &m_m, const Vector &d_m, double &pre_fact,
                     const VectorGeneric<6> &inv_cov_md, Vector &ov_der);
// calculate exp of overlap for neighbor list update
  double get_exp_overlap(const Vector &m_m, const Vector &d_m,
                         const VectorGeneric<6> &inv_cov_md);
// update the neighbor list
  void update_neighbor_list();
// calculate overlap
  void calculate_overlap();
// non-marginal version
  void calculate_sigma();
// marginal version
  void calculate_marginal();

public:
  static void registerKeywords( Keywords& keys );
  explicit EMMI(const ActionOptions&);
// active methods:
  void prepare();
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(EMMI,"EMMI")

void EMMI::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOMS","atoms for which we calculate the density map, typically all heavy atoms");
  keys.add("compulsory","GMM_FILE","file with the parameters of the GMM components");
  keys.add("compulsory","TEMP","temperature");
  keys.addFlag("NO_AVER",false,"don't do ensemble averaging in multi-replica mode");
  keys.addFlag("ANALYSIS",false,"run in analysis mode");
  keys.addFlag("SAMPLING",false,"do explicit sampling in uncertainty");
  keys.add("compulsory","NL_CUTOFF","The cutoff in overlap for the neighbor list");
  keys.add("compulsory","NL_STRIDE","The frequency with which we are updating the neighbor list");
  keys.add("compulsory","SIGMA_MEAN_H","the (hot) uncertainty in the mean estimate");
  keys.add("compulsory","SIGMA_MEAN_C","the (cold) uncertainty in the mean estimate");
  keys.add("optional","SIGMA0","initial value of the uncertainty");
  keys.add("optional","DSIGMA","MC step for uncertainties");
  keys.add("optional","MC_STRIDE","Monte Carlo stride");
  keys.add("optional","ERR_FILE","file with experimental overlaps");
  keys.add("optional","STATUS_FILE","write a file with all the data usefull for restart");
  keys.add("optional","WRITE_STRIDE","write the status to a file every N steps, this can be used for restart");
  keys.add("optional","PRIOR", "p(sigma)=1/sigma^n, where n = 2*prior-1");
  componentsAreNotOptional(keys);
  keys.addOutputComponent("scoreb", "default","Bayesian score");
  keys.addOutputComponent("acc",    "SAMPLING","MC acceptance");
}

EMMI::EMMI(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  inv_sqrt2_(0.707106781186548),
  sqrt2_pi_(0.797884560802865),
  nl_cutoff_(-1.0), nl_stride_(0),
  first_time_(true), no_aver_(false),
  analysis_(false), nframe_(0.0), pbc_(true),
  MCstride_(1), MCfirst_(-1), MCaccept_(0),
  first_status_(true), do_sampling_(false), prior_(1.0)
{
  // marginal or non-marginal version of the score
  parseFlag("SAMPLING",do_sampling_);

  bool nopbc=!pbc_;
  parseFlag("NOPBC",nopbc);
  pbc_=!nopbc;

  vector<AtomNumber> atoms;
  parseAtomList("ATOMS", atoms);

  string GMM_file;
  parse("GMM_FILE", GMM_file);

  // uncertainty in the mean estimate
  // hot (GMM_beta=1) and cold (GMM_beta=0)
  double sigma_mean_h;
  parse("SIGMA_MEAN_H", sigma_mean_h);
  double sigma_mean_c;
  parse("SIGMA_MEAN_C", sigma_mean_c);

  // initial value of the uncertainty
  double sigma_ini;
  parse("SIGMA0", sigma_ini);
  if(do_sampling_ && sigma_ini<=0) error("with SAMPLING you must specify a positive SIGMA0");

  // MC stuff
  double dsigma;
  parse("DSIGMA", dsigma);
  if(do_sampling_ && dsigma<=0) error("with SAMPLING you must specify a positive DSIGMA");
  parse("MC_STRIDE", MCstride_);

  // error file
  string errfile;
  parse("ERR_FILE", errfile);

  // temperature
  double temp=0.0;
  parse("TEMP",temp);
  // convert temp to kbt
  if(temp>0.0) kbt_=plumed.getAtoms().getKBoltzmann()*temp;
  else kbt_=plumed.getAtoms().getKbT();

  // prior exponent
  parse("PRIOR",prior_);

  // neighbor list stuff
  parse("NL_CUTOFF",nl_cutoff_);
  if(nl_cutoff_<=0.0) error("NL_CUTOFF should be explicitly specified and positive");
  parse("NL_STRIDE",nl_stride_);
  if(nl_stride_<=0) error("NL_STRIDE should be explicitly specified and positive");

  // various flags
  parseFlag("NO_AVER",no_aver_);
  parseFlag("ANALYSIS",analysis_);

  // writing status file
  parse("WRITE_STRIDE", statusstride_);
  if(do_sampling_ && statusstride_==0) error("with SAMPLING you must specify a positive WRITE_STRIDE");

  parse("STATUS_FILE",  statusfilename_);
  if(statusfilename_=="") statusfilename_ = "MISTATUS"+getLabel();
  else                    statusfilename_ = statusfilename_+getLabel();

  checkRead();

  // set parallel stuff
  size_=comm.Get_size();
  rank_=comm.Get_rank();

  // get number of replicas
  if(rank_==0) {
    if(no_aver_) {
      nrep_ = 1;
    } else {
      nrep_ = multi_sim_comm.Get_size();
    }
    replica_ = multi_sim_comm.Get_rank();
  } else {
    nrep_ = 0;
    replica_ = 0;
  }
  comm.Sum(&nrep_,1);
  comm.Sum(&replica_,1);

  log.printf("  atoms involved : ");
  for(unsigned i=0; i<atoms.size(); ++i) log.printf("%d ",atoms[i].serial());
  log.printf("\n");
  log.printf("  GMM data file : %s\n", GMM_file.c_str());
  if(no_aver_) log.printf("  without ensemble averaging\n");
  log.printf("  neighbor list overlap cutoff : %lf\n", nl_cutoff_);
  log.printf("  neighbor list stride : %u\n",  nl_stride_);
  log.printf("  (hot) uncertainty in the mean estimate : %f\n",sigma_mean_h);
  log.printf("  (cold) uncertainty in the mean estimate : %f\n",sigma_mean_c);
  if(do_sampling_) {
    log.printf("  initial value of the uncertainty : %f\n",sigma_ini);
    log.printf("  max MC move in uncertainty : %f\n",dsigma);
    log.printf("  MC stride : %u\n", MCstride_);
    log.printf("  reading/writing to status file : %s\n",statusfilename_.c_str());
    log.printf("  with stride : %u\n",statusstride_);
    log.printf("  prior exponent : %f\n",prior_);
  }
  if(errfile.size()>0) log.printf("  reading experimental overlaps from file : %s\n", errfile.c_str());
  log.printf("  temperature of the system in energy unit : %f\n",kbt_);
  log.printf("  number of replicas for averaging: %u\n",nrep_);
  log.printf("  id of the replica : %u\n",replica_);

  // set constant quantity before calculating stuff
  cfact_ = 1.0/pow( 2.0*pi, 1.5 );

  // calculate model GMM constant parameters
  vector<double> GMM_m_w = get_GMM_m(atoms);

  // read data GMM parameters
  get_GMM_d(GMM_file);
  log.printf("  number of GMM components : %u\n", static_cast<unsigned>(GMM_d_m_.size()));

  // normalize atom weight map
  double norm_d = accumulate(GMM_d_w_.begin(), GMM_d_w_.end(), 0.0);
  double norm_m = accumulate(GMM_m_w.begin(),  GMM_m_w.end(),  0.0);
  for(unsigned i=0; i<GMM_m_w_.size(); ++i) GMM_m_w_[i] *= norm_d / norm_m;

  // read experimental errors
  vector<double> exp_err;
  if(errfile.size()>0) exp_err = read_exp_errors(errfile);

  // get self overlaps between data GMM components
  // retrieve error and set sampling parameters
  double s0_ave = 0.0;
  vector<double> s0_median;
  for(unsigned i=0; i<GMM_d_m_.size(); ++i) {
    double ov = get_self_overlap(i);
    ovdd_.push_back(ov);
    // retrieve experimental error, if present
    double s0_exp = 0.0;
    if(errfile.size()>0) s0_exp = exp_err[i];
    // calculate average and median relative s0_exp
    s0_ave += s0_exp / ov;
    s0_median.push_back(s0_exp/ov);
    // add sigma_mean contribution
    if(GMM_d_beta_[i]==1) sigma_mean_.push_back(sigma_mean_h*ov);
    if(GMM_d_beta_[i]==0) sigma_mean_.push_back(sigma_mean_c*ov);
    // for non marginal version
    if(do_sampling_) {
      // add minimum value of sigma
      sigma_min_.push_back(s0_exp);
      // add MC move in sigma
      dsigma_.push_back(dsigma*ov);
      // set sigma max
      sigma_max_.push_back(1000.*ov);
      // initialize sigma
      sigma_.push_back(std::max(sigma_min_[i], std::min(sigma_max_[i], sigma_ini*ov)));
    } else {
      // for marginal version
      sigma0_.push_back(sqrt(s0_exp*s0_exp+sigma_mean_[i]*sigma_mean_[i]));
    }
  }
  // final calculation of average and median
  s0_ave /= static_cast<double>(GMM_d_m_.size());
  std::sort(s0_median.begin(), s0_median.end());
  if(errfile.size()>0) {
    log.printf("  average relative error : %f\n", s0_ave);
    log.printf("  median relative error  : %f\n", s0_median[s0_median.size()/2]);
  }
  // read status file if restarting
  if(do_sampling_ && getRestart()) read_status();

  // calculate auxiliary stuff
  calculate_useful_stuff();

  // prepare data and derivative vectors
  ovmd_.resize(GMM_d_m_.size());
  atom_der_.resize(GMM_m_type_.size());
  if(!do_sampling_) {
    err_f_.resize(GMM_d_w_.size());
    exp_f_.resize(GMM_d_w_.size());
  }

  // clear things that are no longer needed
  GMM_d_cov_.clear();
  GMM_d_w_.clear();
  ovdd_cut_.clear();

  // add components
  addComponentWithDerivatives("scoreb"); componentIsNotPeriodic("scoreb");
  if(do_sampling_) {
    addComponent("acc"); componentIsNotPeriodic("acc");
  }

  // request the atoms
  requestAtoms(atoms);

  log<<"  Bibliography "<<plumed.cite("Bonomi, Camilloni, Cavalli, Vendruscolo, Sci. Adv. 2, e150117 (2016)");
  log<<plumed.cite("Hanot, Bonomi, Greenberg, Sali, Nilges, Vendruscolo, Pellarin, bioRxiv doi: 10.1101/113951 (2017)");
  log<<"\n";
}

void EMMI::read_status()
{
  double MDtime;
// open file
  IFile *ifile = new IFile();
  if(ifile->FileExist(statusfilename_)) {
    ifile->open(statusfilename_);
    while(ifile->scanField("MD_time", MDtime)) {
      for(unsigned i=0; i<sigma_.size(); ++i) {
        // convert i to string
        std::string num; Tools::convert(i,num);
        // read entries
        ifile->scanField("s"+num, sigma_[i]);
      }
      // new line
      ifile->scanField();
    }
    ifile->close();
  } else {
    error("Cannot find status file "+statusfilename_+"\n");
  }
  delete ifile;
}

void EMMI::print_status(long int step)
{
// if first time open the file
  if(first_status_) {
    first_status_ = false;
    statusfile_.link(*this);
    statusfile_.open(statusfilename_);
    statusfile_.setHeavyFlush();
    statusfile_.fmtField("%6.3e ");
  }
// write fields
  double MDtime = static_cast<double>(step)*getTimeStep();
  statusfile_.printField("MD_time", MDtime);
  for(unsigned i=0; i<sigma_.size(); ++i) {
    // convert i to string
    std::string num; Tools::convert(i,num);
    // print entry
    statusfile_.printField("s"+num, sigma_[i]);
  }
  statusfile_.printField();
}

double EMMI::proposeMove(double x, double xmin, double xmax, double dxmax)
{
  double r = static_cast<double>(rand()) / RAND_MAX;
  double dx = -dxmax + r * 2.0 * dxmax;
  double x_new = x + dx;
// check boundaries
  if(x_new > xmax) {x_new = 2.0 * xmax - x_new;}
  if(x_new < xmin) {x_new = 2.0 * xmin - x_new;}
  return x_new;
}

bool EMMI::doAccept(double oldE, double newE) {
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

void EMMI::doMonteCarlo(long int step)
{
// prepare new sigma vector
  vector<double> new_sigma(sigma_.size(), 0.0);
// and acceptance
  double MCaccept = 0.0;

// cycle on sigmas - in parallel
  for(unsigned i=rank_; i<sigma_.size(); i+=size_) {
    // store prefactor
    double pre_fact = 0.5*kbt_*( ovmd_[i]-ovdd_[i] )*( ovmd_[i]-ovdd_[i] );
    // old stuff
    double old_s2 = sigma_mean_[i]*sigma_mean_[i]+sigma_[i]*sigma_[i];
    // deviations from data, normalization and prior
    double old_ene = pre_fact/old_s2 + kbt_*prior_*std::log(old_s2);
    // propose move
    double new_s = proposeMove(sigma_[i], sigma_min_[i], sigma_max_[i], dsigma_[i]);
    // new stuff
    double new_s2 = sigma_mean_[i]*sigma_mean_[i]+new_s*new_s;
    // deviations from data, normalization and Jeffret's prior
    double new_ene = pre_fact/new_s2 + kbt_*prior_*std::log(new_s2);
    // accept or reject
    new_sigma[i] = sigma_[i];
    bool accept = doAccept(old_ene, new_ene);
    if(accept) {
      new_sigma[i] = new_s;
      MCaccept++;
    }
  }
// collect all new sigmas
  comm.Sum(&new_sigma[0], new_sigma.size());
// collect acceptances
  comm.Sum(&MCaccept, 1);
// overwrite old sigmas
  for(unsigned i=0; i<sigma_.size(); ++i) sigma_[i] = new_sigma[i];
// increment total acceptances
  MCaccept_ += MCaccept;
}

vector<double> EMMI::read_exp_errors(string errfile)
{
  int nexp, idcomp;
  double err;
  vector<double> exp_err;
// open file
  IFile *ifile = new IFile();
  if(ifile->FileExist(errfile)) {
    ifile->open(errfile);
    // scan for number of experimental errors
    ifile->scanField("Nexp", nexp);
    // cycle on GMM components
    while(ifile->scanField("Id",idcomp)) {
      // cycle on number of experimental errors
      double err_tot = 0.0;
      for(unsigned i=0; i<nexp; ++i) {
        string ss; Tools::convert(i,ss);
        ifile->scanField("Err"+ss, err);
        err_tot += err*err;
      }
      // new line
      ifile->scanField();
      // calculate root mean squared error
      err_tot = sqrt(err_tot/static_cast<double>(nexp));
      // add to global vector
      exp_err.push_back(err_tot);
    }
    ifile->close();
  } else {
    error("Cannot find ERR_FILE "+errfile+"\n");
  }
  return exp_err;
}

vector<double> EMMI::get_GMM_m(vector<AtomNumber> &atoms)
{
  // list of weights - one per atom
  vector<double> GMM_m_w;

  vector<SetupMolInfo*> moldat=plumed.getActionSet().select<SetupMolInfo*>();
  // map of atom types to A and B coefficients of scattering factor
  // f(s) = A * exp(-B*s**2)
  // B is in Angstrom squared
  // map between an atom type and an index
  map<string, unsigned> type_map;
  type_map["C"]=0;
  type_map["O"]=1;
  type_map["N"]=2;
  type_map["S"]=3;
  // fill in sigma vector
  GMM_m_s_.push_back(15.146);  // type 0
  GMM_m_s_.push_back(8.59722); // type 1
  GMM_m_s_.push_back(11.1116); // type 2
  GMM_m_s_.push_back(15.8952); // type 3
  // fill in weight vector
  GMM_m_w_.push_back(2.49982); // type 0
  GMM_m_w_.push_back(1.97692); // type 1
  GMM_m_w_.push_back(2.20402); // type 2
  GMM_m_w_.push_back(5.14099); // type 3

  // check if MOLINFO line is present
  if( moldat.size()==1 ) {
    log<<"  MOLINFO DATA found, using proper atom names\n";
    for(unsigned i=0; i<atoms.size(); ++i) {
      // get atom name
      string name = moldat[0]->getAtomName(atoms[i]);
      char type;
      // get atom type
      char first = name.at(0);
      // GOLDEN RULE: type is first letter, if not a number
      if (!isdigit(first)) {
        type = first;
        // otherwise is the second
      } else {
        type = name.at(1);
      }
      // check if key in map
      std::string type_s = std::string(1,type);
      if(type_map.find(type_s) != type_map.end()) {
        // save atom type
        GMM_m_type_.push_back(type_map[type_s]);
        // this will be normalized in the final density
        GMM_m_w.push_back(GMM_m_w_[type_map[type_s]]);
      } else {
        error("Wrong atom type "+type_s+" from atom name "+name+"\n");
      }
    }
  } else {
    error("MOLINFO DATA not found\n");
  }
  return GMM_m_w;
}

void EMMI::check_GMM_d(VectorGeneric<6> &cov, double w)
{

// check if positive defined, by calculating the 3 leading principal minors
  double pm1 = cov[0];
  double pm2 = cov[0]*cov[3]-cov[1]*cov[1];
  double pm3 = cov[0]*(cov[3]*cov[5]-cov[4]*cov[4])-cov[1]*(cov[1]*cov[5]-cov[4]*cov[2])+cov[2]*(cov[1]*cov[4]-cov[3]*cov[2]);
// apply Sylvester’s criterion
  if(pm1<=0.0 || pm2<=0.0 || pm3<=0.0)
    error("check data GMM: covariance matrix is not positive defined");

// check if weight is positive
  if(w<0.0) error("check data GMM: weight must be positive");
}

// read GMM data file in PLUMED format:
void EMMI::get_GMM_d(string GMM_file)
{
  int idcomp, beta;
  double w, m0, m1, m2;
  VectorGeneric<6> cov;

// open file
  IFile *ifile = new IFile();
  if(ifile->FileExist(GMM_file)) {
    ifile->open(GMM_file);
    while(ifile->scanField("Id",idcomp)) {
      ifile->scanField("Weight",w);
      ifile->scanField("Mean_0",m0);
      ifile->scanField("Mean_1",m1);
      ifile->scanField("Mean_2",m2);
      ifile->scanField("Cov_00",cov[0]);
      ifile->scanField("Cov_01",cov[1]);
      ifile->scanField("Cov_02",cov[2]);
      ifile->scanField("Cov_11",cov[3]);
      ifile->scanField("Cov_12",cov[4]);
      ifile->scanField("Cov_22",cov[5]);
      ifile->scanField("Beta",beta);
      // check input
      check_GMM_d(cov, w);
      // check beta
      if(beta!=0 && beta!=1) error("Beta must be either 0 or 1");
      // center of the Gaussian
      GMM_d_m_.push_back(Vector(m0,m1,m2));
      // covariance matrix
      GMM_d_cov_.push_back(cov);
      // weight
      GMM_d_w_.push_back(w);
      // beta
      GMM_d_beta_.push_back(beta);
      // new line
      ifile->scanField();
    }
    ifile->close();
  } else {
    error("Cannot find GMM_FILE "+GMM_file+"\n");
  }
  delete ifile;
}

void EMMI::calculate_useful_stuff()
{
  VectorGeneric<6> cov, sum, inv_sum;
  // cycle on all atoms types
  for(unsigned i=0; i<GMM_m_s_.size(); ++i) {
    // the Gaussian in density (real) space is the FT of scattering factor
    // f(r) = A * (pi/B)**1.5 * exp(-pi**2/B*r**2)
    double s = sqrt ( 0.5 * GMM_m_s_[i] ) / pi * 0.1;
    // covariance matrix for spherical Gaussian
    cov[0]=s*s; cov[1]=0.0; cov[2]=0.0;
    cov[3]=s*s; cov[4]=0.0;
    cov[5]=s*s;
    // cycle on all data GMM
    for(unsigned j=0; j<GMM_d_m_.size(); ++j) {
      // we need the sum of the covariance matrices
      for(unsigned k=0; k<6; ++k) sum[k] = cov[k] + GMM_d_cov_[j][k];
      // and to calculate its determinant
      double det = sum[0]*(sum[3]*sum[5]-sum[4]*sum[4]);
      det -= sum[1]*(sum[1]*sum[5]-sum[4]*sum[2]);
      det += sum[2]*(sum[1]*sum[4]-sum[3]*sum[2]);
      // calculate prefactor
      double pre_fact =  cfact_ / sqrt(det) * GMM_d_w_[j] * GMM_m_w_[i];
      // and its inverse
      inv_sum[0] = (sum[3]*sum[5] - sum[4]*sum[4])/det;
      inv_sum[1] = (sum[2]*sum[4] - sum[1]*sum[5])/det;
      inv_sum[2] = (sum[1]*sum[4] - sum[2]*sum[3])/det;
      inv_sum[3] = (sum[0]*sum[5] - sum[2]*sum[2])/det;
      inv_sum[4] = (sum[2]*sum[1] - sum[0]*sum[4])/det;
      inv_sum[5] = (sum[0]*sum[3] - sum[1]*sum[1])/det;
      // now we store the prefactor
      pre_fact_.push_back(pre_fact);
      // and the inverse of the sum
      inv_cov_md_.push_back(inv_sum);
    }
  }
  // get cutoff for overlap calculation - avoid millions of exp calculations
  get_cutoff_ov();
}

// get prefactors
double EMMI::get_prefactor_inverse
(const VectorGeneric<6> &GMM_cov_0, const VectorGeneric<6> &GMM_cov_1,
 double &GMM_w_0, double &GMM_w_1,
 VectorGeneric<6> &sum, VectorGeneric<6> &inv_sum)
{
// we need the sum of the covariance matrices
  for(unsigned k=0; k<6; ++k) sum[k] = GMM_cov_0[k] + GMM_cov_1[k];

// and to calculate its determinant
  double det = sum[0]*(sum[3]*sum[5]-sum[4]*sum[4]);
  det -= sum[1]*(sum[1]*sum[5]-sum[4]*sum[2]);
  det += sum[2]*(sum[1]*sum[4]-sum[3]*sum[2]);

// the prefactor is
  double pre_fact =  cfact_ / sqrt(det) * GMM_w_0 * GMM_w_1;

// and its inverse
  inv_sum[0] = (sum[3]*sum[5] - sum[4]*sum[4])/det;
  inv_sum[1] = (sum[2]*sum[4] - sum[1]*sum[5])/det;
  inv_sum[2] = (sum[1]*sum[4] - sum[2]*sum[3])/det;
  inv_sum[3] = (sum[0]*sum[5] - sum[2]*sum[2])/det;
  inv_sum[4] = (sum[2]*sum[1] - sum[0]*sum[4])/det;
  inv_sum[5] = (sum[0]*sum[3] - sum[1]*sum[1])/det;

// return pre-factor
  return pre_fact;
}

double EMMI::get_self_overlap(unsigned id)
{
  double ov_tot = 0.0;
  VectorGeneric<6> sum, inv_sum;
  Vector ov_der;
// start loop
  for(unsigned i=0; i<GMM_d_m_.size(); ++i) {
    // call auxiliary method
    double pre_fact = get_prefactor_inverse(GMM_d_cov_[id], GMM_d_cov_[i],
                                            GMM_d_w_[id],   GMM_d_w_[i], sum, inv_sum);
    // add overlap to ov_tot
    ov_tot += get_overlap(GMM_d_m_[id], GMM_d_m_[i], pre_fact, inv_sum, ov_der);
  }
// store cutoff
  ovdd_cut_.push_back(ov_tot * nl_cutoff_);
// and return it
  return ov_tot;
}

// this is to avoid the calculation of millions of exp function
// when updating the neighbor list using calculate_overlap
void EMMI::get_cutoff_ov()
{
  // prepare ov_cut_
  for(unsigned i=0; i<pre_fact_.size(); ++i) ov_cut_.push_back(0.0);
  // temp stuff
  unsigned GMM_d_size = GMM_d_m_.size();
  unsigned GMM_m_size = GMM_m_type_.size();
  // cycle on all overlaps
  unsigned nover = GMM_d_size * GMM_m_size;
  for(unsigned k=0; k<nover; ++k) {
    // get data (id) and atom (im) indexes
    unsigned id = k / GMM_m_size;
    unsigned im = k % GMM_m_size;
    // get index in auxiliary lists
    unsigned kaux = GMM_m_type_[im] * GMM_d_size + id;
    // store cutoff for exponent of the overlap
    ov_cut_[kaux] = -2.0 * std::log( ovdd_cut_[id] / pre_fact_[kaux] );
  }
}

// get overlap and derivatives
double EMMI::get_overlap(const Vector &m_m, const Vector &d_m, double &pre_fact,
                         const VectorGeneric<6> &inv_cov_md, Vector &ov_der)
{
  Vector md;
  // calculate vector difference m_m-d_m with/without pbc
  if(pbc_) md = pbcDistance(d_m, m_m);
  else     md = delta(d_m, m_m);
  // calculate product of transpose of md and inv_cov_md
  double p_x = md[0]*inv_cov_md[0]+md[1]*inv_cov_md[1]+md[2]*inv_cov_md[2];
  double p_y = md[0]*inv_cov_md[1]+md[1]*inv_cov_md[3]+md[2]*inv_cov_md[4];
  double p_z = md[0]*inv_cov_md[2]+md[1]*inv_cov_md[4]+md[2]*inv_cov_md[5];
  // calculate product of prod and md
  double ov = md[0]*p_x+md[1]*p_y+md[2]*p_z;
  // final calculation
  ov = pre_fact * exp(-0.5*ov);
  // derivatives
  ov_der = ov * Vector(p_x, p_y, p_z);
  return ov;
}

// get the exponent of the overlap
double EMMI::get_exp_overlap(const Vector &m_m, const Vector &d_m,
                             const VectorGeneric<6> &inv_cov_md)
{
  Vector md;
  // calculate vector difference m_m-d_m with/without pbc
  if(pbc_) md = pbcDistance(d_m, m_m);
  else     md = delta(d_m, m_m);
  // calculate product of transpose of md and inv_cov_md
  double p_x = md[0]*inv_cov_md[0]+md[1]*inv_cov_md[1]+md[2]*inv_cov_md[2];
  double p_y = md[0]*inv_cov_md[1]+md[1]*inv_cov_md[3]+md[2]*inv_cov_md[4];
  double p_z = md[0]*inv_cov_md[2]+md[1]*inv_cov_md[4]+md[2]*inv_cov_md[5];
  // calculate product of prod and md
  double ov = md[0]*p_x+md[1]*p_y+md[2]*p_z;
  return ov;
}

void EMMI::update_neighbor_list()
{
  // temp stuff
  unsigned GMM_d_size = GMM_d_m_.size();
  unsigned GMM_m_size = GMM_m_type_.size();
  // local neighbor list
  vector < unsigned > nl_l;
  // clear old neighbor list
  nl_.clear();
  // cycle on all overlaps (in parallel)
  unsigned nover = GMM_d_size * GMM_m_size;
  for(unsigned k=rank_; k<nover; k=k+size_) {
    // get data (id) and atom (im) indexes
    unsigned id = k / GMM_m_size;
    unsigned im = k % GMM_m_size;
    // get index in auxiliary lists
    unsigned kaux = GMM_m_type_[im] * GMM_d_size + id;
    // calculate exponent of overlap
    double expov = get_exp_overlap(GMM_d_m_[id], getPosition(im), inv_cov_md_[kaux]);
    // fill the neighbor list
    if(expov <= ov_cut_[kaux]) nl_l.push_back(k);
  }
  // find total dimension of neighborlist
  vector <int> recvcounts(size_, 0);
  recvcounts[rank_] = nl_l.size();
  comm.Sum(&recvcounts[0], size_);
  int tot_size = accumulate(recvcounts.begin(), recvcounts.end(), 0);
  // resize neighbor stuff
  nl_.resize(tot_size);
  // calculate vector of displacement
  vector<int> disp(size_);
  disp[0] = 0;
  int rank_size = 0;
  for(unsigned i=0; i<size_-1; ++i) {
    rank_size += recvcounts[i];
    disp[i+1] = rank_size;
  }
  // Allgather neighbor list
  comm.Allgatherv(&nl_l[0], recvcounts[rank_], &nl_[0], &recvcounts[0], &disp[0]);
  // now resize derivatives
  ovmd_der_.resize(tot_size);
}

void EMMI::prepare()
{
  if(getExchangeStep()) first_time_=true;
}

// overlap calculator
void EMMI::calculate_overlap() {

  if(first_time_ || getExchangeStep() || getStep()%nl_stride_==0) {
    update_neighbor_list();
    first_time_=false;
  }

  // clean temporary vectors
  for(unsigned i=0; i<ovmd_.size(); ++i)     ovmd_[i] = 0.0;
  for(unsigned i=0; i<ovmd_der_.size(); ++i) ovmd_der_[i] = Vector(0,0,0);

  // we have to cycle over all model and data GMM components in the neighbor list
  unsigned GMM_d_size = GMM_d_m_.size();
  unsigned GMM_m_size = GMM_m_type_.size();
  for(unsigned i=rank_; i<nl_.size(); i=i+size_) {
    // get data (id) and atom (im) indexes
    unsigned id = nl_[i] / GMM_m_size;
    unsigned im = nl_[i] % GMM_m_size;
    // get index in auxiliary lists
    unsigned kaux = GMM_m_type_[im] * GMM_d_size + id;
    // add overlap with im component of model GMM
    ovmd_[id] += get_overlap(GMM_d_m_[id], getPosition(im), pre_fact_[kaux],
                             inv_cov_md_[kaux], ovmd_der_[i]);
  }
  // communicate stuff
  comm.Sum(&ovmd_[0], ovmd_.size());
  comm.Sum(&ovmd_der_[0][0], 3*ovmd_der_.size());
}


void EMMI::calculate() {

// calculate CV
  calculate_overlap();

  if(!analysis_) {

    // BIASING MODE
    // sampling sigma or marginal version?
    if(do_sampling_) calculate_sigma();
    else             calculate_marginal();

  } else {

    // ANALYSIS MODE
    // prepare stuff for the first time
    if(nframe_ <= 0.0) {
      Devfile_.link(*this);
      Devfile_.open("ovmd_deviations.dat");
      Devfile_.setHeavyFlush();
      Devfile_.fmtField("%12.6f");
      ovmd_ave_.resize(ovmd_.size());
      for(unsigned i=0; i<ovmd_ave_.size(); ++i) ovmd_ave_[i] = 0.0;
    }

    // increment number of frames
    nframe_ += 1.0;

    // add average ovmd_
    for(unsigned i=0; i<ovmd_.size(); ++i) ovmd_ave_[i] += ovmd_[i];

    // print stuff
    for(unsigned i=0; i<ovmd_.size(); ++i) {
      // convert i to string
      string ss; Tools::convert(i,ss);
      // print entry
      double ave = ovmd_ave_[i] / nframe_;
      double dev2 = (ave-ovdd_[i])*(ave-ovdd_[i])/ovdd_[i]/ovdd_[i];
      double dev = sqrt(dev2);
      Devfile_.printField("ovmd_" + ss, dev);
    }
    Devfile_.printField();
  }

}

void EMMI::calculate_sigma()
{
  // non-marginal version

  // rescale factor for ensemble average
  double escale = 1.0 / static_cast<double>(nrep_);

  // prepare vector of inv_s2
  vector<double> inv_s2;
  for(unsigned i=0; i<ovmd_.size(); ++i)
    inv_s2.push_back( 1.0 / ( sigma_mean_[i]*sigma_mean_[i]+sigma_[i]*sigma_[i] ) );

  // calculate average of ovmd_ across replicas
  // and collect sigma contributions from replicas
  if(!no_aver_ && nrep_>1) {
    if(rank_==0) {
      multi_sim_comm.Sum(&ovmd_[0],  ovmd_.size());
      multi_sim_comm.Sum(&inv_s2[0], inv_s2.size());
      for(unsigned i=0; i<ovmd_.size(); ++i) ovmd_[i] *= escale;
    } else {
      for(unsigned i=0; i<ovmd_.size(); ++i) {
        ovmd_[i]  = 0.0;
        inv_s2[i] = 0.0;
      }
    }
    comm.Sum(&ovmd_[0],  ovmd_.size());
    comm.Sum(&inv_s2[0], ovmd_.size());
  }

  // calculate score and reweighting score
  double ene = 0.0;
  for(unsigned i=0; i<ovmd_.size(); ++i) {
    // increment energy
    ene += ( ovmd_[i]-ovdd_[i] ) * ( ovmd_[i]-ovdd_[i] ) * inv_s2[i];
  }

  // multiply by constant factors
  ene *= 0.5 * kbt_;

  // clean temporary vector
  for(unsigned i=0; i<atom_der_.size(); ++i) atom_der_[i] = Vector(0,0,0);

  // virial
  Tensor virial;
  // get derivatives of bias with respect to atoms
  for(unsigned i=rank_; i<nl_.size(); i=i+size_) {
    // get indexes of data and model component
    unsigned id = nl_[i] / GMM_m_type_.size();
    unsigned im = nl_[i] % GMM_m_type_.size();
    // derivative
    double der = kbt_ * ( ovmd_[id]-ovdd_[id] ) * inv_s2[id];
    // chain rule + replica normalization
    Vector tot_der = der * ovmd_der_[i] * escale;
    // atom's position in GMM cell
    Vector pos;
    if(pbc_) pos = pbcDistance(GMM_d_m_[id], getPosition(im)) + GMM_d_m_[id];
    else     pos = getPosition(im);
    // increment derivative and virial
    atom_der_[im] += tot_der;
    virial += Tensor(pos, -tot_der);
  }

  // communicate stuff
  comm.Sum(&atom_der_[0][0], 3*atom_der_.size());
  comm.Sum(virial);

  // set derivatives, virial, and score
  for(unsigned i=0; i<atom_der_.size(); ++i) setAtomsDerivatives(getPntrToComponent("scoreb"), i, atom_der_[i]);
  setBoxDerivatives(getPntrToComponent("scoreb"), virial);
  getPntrToComponent("scoreb")->set(ene);

  // get time step
  long int step = getStep();

  // move sigmas
  if(step%MCstride_==0) doMonteCarlo(step);

  // print status
  if(step%statusstride_==0) print_status(step);

  // calculate acceptance ratio
  // this is needed when restarting simulations
  if(MCfirst_==-1) MCfirst_=step;
  double MCtrials = std::floor(static_cast<double>(step-MCfirst_) / static_cast<double>(MCstride_))+1.0;
  // average acceptance across all sigmas
  double acc = static_cast<double>(MCaccept_) / MCtrials / static_cast<double>(sigma_.size());
  // set value
  getPntrToComponent("acc")->set(acc);

}

void EMMI::calculate_marginal()
{
  // marginal version
  // rescale factor for ensemble average
  double escale = 1.0 / static_cast<double>(nrep_);

  // calculate average of ovmd_ across replicas
  if(!no_aver_ && nrep_>1) {
    if(comm.Get_rank()==0) {
      multi_sim_comm.Sum(&ovmd_[0], ovmd_.size());
      for(unsigned i=0; i<ovmd_.size(); ++i) ovmd_[i] *= escale;
    } else {
      for(unsigned i=0; i<ovmd_.size(); ++i) ovmd_[i]  = 0.0;
    }
    comm.Sum(&ovmd_[0], ovmd_.size());
  }

  // calculate score
  double ene = 0.0;
  for(unsigned i=0; i<ovmd_.size(); ++i) {
    // useful quantity
    double dev = ( ovmd_[i]-ovdd_[i] ) / sigma0_[i];
    // calculate and store err and exp function
    err_f_[i] = erf ( dev * inv_sqrt2_ );
    exp_f_[i] = exp( - 0.5 * dev * dev );
    // increment energy
    ene += -std::log ( 0.5 / dev * err_f_[i]);
  }

  // multiply by constant factors
  ene *= kbt_ / escale;

  // clean temporary vector
  for(unsigned i=0; i<atom_der_.size(); ++i) atom_der_[i] = Vector(0,0,0);

  // virial
  Tensor virial;
  // get derivatives of bias with respect to atoms
  for(unsigned i=rank_; i<nl_.size(); i=i+size_) {
    // get indexes of data and model component
    unsigned id = nl_[i] / GMM_m_type_.size();
    unsigned im = nl_[i] % GMM_m_type_.size();
    // first part of derivative
    double der = -kbt_ / err_f_[id] * sqrt2_pi_ * exp_f_[id] / sigma0_[id];
    // second part
    der += kbt_ / (ovmd_[id]-ovdd_[id]);
    // chain rule
    Vector tot_der = der * ovmd_der_[i];
    // atom's position in GMM cell
    Vector pos;
    if(pbc_) pos = pbcDistance(GMM_d_m_[id], getPosition(im)) + GMM_d_m_[id];
    else     pos = getPosition(im);
    // increment derivative and virial
    atom_der_[im] += tot_der;
    virial += Tensor(pos, -tot_der);
  }

  // communicate stuff
  comm.Sum(&atom_der_[0][0], 3*atom_der_.size());
  comm.Sum(virial);

  // set derivatives, virial, and score
  for(unsigned i=0; i<atom_der_.size(); ++i) setAtomsDerivatives(getPntrToComponent("scoreb"), i, atom_der_[i]);
  setBoxDerivatives(getPntrToComponent("scoreb"), virial);
  getPntrToComponent("scoreb")->set(ene);
}

}
}
