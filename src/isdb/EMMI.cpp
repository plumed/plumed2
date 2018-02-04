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
#include "tools/Random.h"

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
  vector < vector<int> >     GMM_d_grps_;
// overlaps
  vector<double> ovmd_;
  vector<double> ovdd_;
  vector<double> ovmd_ave_;
  vector<double> ov_cut_;
// and derivatives
  vector<Vector> ovmd_der_;
  vector<Vector> atom_der_;
// constant quantities;
  double cfact_;
// metainference
  unsigned nrep_;
  unsigned replica_;
  vector<double> sigma_;
  vector<double> sigma_mean_;
  vector<double> sigma_min_;
  vector<double> sigma_max_;
  vector<double> dsigma_;
// list of prefactors for overlap between two components of model and data GMM
// pre_fact = 1.0 / (2pi)**1.5 / sqrt(det_md) * Wm * Wd
  vector<double> pre_fact_;
// inverse of the sum of model and data covariances matrices
  vector< VectorGeneric<6> > inv_cov_md_;
// neighbor list
  double   nl_cutoff_;
  unsigned nl_stride_;
  bool first_time_;
  bool no_aver_;
  vector<unsigned> nl_;
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
  int      MCstride_;
  long int MCfirst_;
  unsigned MCaccept_;
  Random   random_;
  // status stuff
  unsigned int statusstride_;
  string       statusfilename_;
  OFile        statusfile_;
  bool         first_status_;
  // regression
  unsigned nregres_;
  double scale_;
  double off_;
  double corr_;
  // tabulated exponential
  double dpcutoff_;
  double dexp_;
  unsigned nexp_;
  vector<double> tab_exp_;
  // simulated annealing
  unsigned nanneal_;
  double   kanneal_;
  double   anneal_;

// annealing
  double get_annealing(long int step);
// do regression
  void doRegression();
// read and write status
  void read_status();
  void print_status(long int step);
// accept or reject
  bool doAccept(double oldE, double newE);
// do MonteCarlo
  void doMonteCarlo();
// read error file
  vector<double> read_exp_errors(string errfile);
// calculate model GMM weights and covariances
  vector<double> get_GMM_m(vector<AtomNumber> &atoms);
// read data GMM file
  void get_GMM_d(string gmm_file);
// check GMM data
  void check_GMM_d(VectorGeneric<6> &cov, double w);
// auxiliary method
  void calculate_useful_stuff(double reso);
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
  keys.add("compulsory","NL_CUTOFF","The cutoff in overlap for the neighbor list");
  keys.add("compulsory","NL_STRIDE","The frequency with which we are updating the neighbor list");
  keys.add("compulsory","SIGMA_MEAN","the uncertainty in the mean estimate");
  keys.add("compulsory","SIGMA0","initial value of the uncertainty");
  keys.add("compulsory","DSIGMA","MC step for uncertainties");
  keys.add("compulsory","RESOLUTION", "Cryo-EM map resolution");
  keys.add("compulsory","WRITE_STRIDE","write the status to a file every N steps, this can be used for restart");
  keys.add("optional","MC_STRIDE", "Monte Carlo stride");
  keys.add("optional","ERR_FILE","file with experimental overlaps");
  keys.add("optional","STATUS_FILE","write a file with all the data useful for restart");
  keys.add("optional","REGRESSION","regression stride");
  keys.add("optional","SCALE","scale factor");
  keys.add("optional","OFFSET","offset");
  keys.add("optional","ANNEAL", "Length of annealing cycle");
  keys.add("optional","ANNEAL_FACT", "Annealing temperature factor");
  keys.add("optional","TEMP","temperature");
  keys.addFlag("NO_AVER",false,"don't do ensemble averaging in multi-replica mode");
  keys.addFlag("ANALYSIS",false,"run in analysis mode");
  componentsAreNotOptional(keys);
  keys.addOutputComponent("scoreb","default","Bayesian score");
  keys.addOutputComponent("acc",   "default","MC acceptance");
  keys.addOutputComponent("scale", "REGRESSION","scaling factor");
  keys.addOutputComponent("offset","REGRESSION","offset");
  keys.addOutputComponent("corr",  "REGRESSION","Pearson's correlation coefficient");
  keys.addOutputComponent("anneal","ANNEAL","annealing factor");
}

EMMI::EMMI(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  first_time_(true), no_aver_(false),
  analysis_(false), nframe_(0.0), pbc_(true),
  MCstride_(1), MCfirst_(-1), MCaccept_(0),
  statusstride_(0), first_status_(true),
  nregres_(0), scale_(1.0), off_(0.0),
  dpcutoff_(15.0), nexp_(1000000),
  nanneal_(0), kanneal_(0.), anneal_(1.)
{
  bool nopbc=!pbc_;
  parseFlag("NOPBC",nopbc);
  pbc_=!nopbc;

  vector<AtomNumber> atoms;
  parseAtomList("ATOMS", atoms);

  string GMM_file;
  parse("GMM_FILE", GMM_file);

  // uncertainty in the mean estimate
  double sigma_mean;
  parse("SIGMA_MEAN", sigma_mean);
  if(sigma_mean<0) error("SIGMA_MEAN should be greater or equal to zero");

  // initial value of the uncertainty
  double sigma_ini;
  parse("SIGMA0", sigma_ini);
  if(sigma_ini<=0) error("you must specify a positive SIGMA0");

  // MC stuff
  double dsigma;
  parse("DSIGMA", dsigma);
  if(dsigma<0) error("you must specify a positive DSIGMA");
  parse("MC_STRIDE", MCstride_);
  if(dsigma>0 && MCstride_<=0) error("you must specify a positive MC_STRIDE");

  // error file
  string errfile;
  parse("ERR_FILE", errfile);

  // temperature
  double temp=0.0;
  parse("TEMP",temp);
  // convert temp to kbt
  if(temp>0.0) kbt_=plumed.getAtoms().getKBoltzmann()*temp;
  else kbt_=plumed.getAtoms().getKbT();

  // simulated annealing stuff
  parse("ANNEAL", nanneal_);
  parse("ANNEAL_FACT", kanneal_);
  if(nanneal_>0 && kanneal_<=1.0) error("with ANNEAL, ANNEAL_FACT must be greater than 1");

  // regression stride
  parse("REGRESSION",nregres_);
  // scale and offset
  parse("SCALE", scale_);
  parse("OFFSET", off_);

  // read map resolution
  double reso;
  parse("RESOLUTION", reso);
  if(reso<=0.) error("RESOLUTION should be strictly positive");

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
  if(statusstride_<=0) error("you must specify a positive WRITE_STRIDE");

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
  log.printf("  neighbor list cutoff : %lf\n", nl_cutoff_);
  log.printf("  neighbor list stride : %u\n",  nl_stride_);
  log.printf("  uncertainty in the mean estimate : %f\n",sigma_mean);
  if(nregres_>0) log.printf("  regression stride : %u\n", nregres_);
  log.printf("  scale factor : %lf\n",scale_);
  log.printf("  constant offset : %lf\n",off_);
  log.printf("  initial value of the uncertainty : %f\n",sigma_ini);
  log.printf("  max MC move in uncertainty : %f\n",dsigma);
  log.printf("  MC stride : %u\n", MCstride_);
  log.printf("  reading/writing to status file : %s\n",statusfilename_.c_str());
  log.printf("  with stride : %u\n",statusstride_);
  if(errfile.size()>0) log.printf("  reading experimental errors from file : %s\n", errfile.c_str());
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

  // normalize atom weight map - not really needed with REGRESSION
  double norm_d = accumulate(GMM_d_w_.begin(), GMM_d_w_.end(), 0.0);
  double norm_m = accumulate(GMM_m_w.begin(),  GMM_m_w.end(),  0.0);
  for(unsigned i=0; i<GMM_m_w_.size(); ++i) GMM_m_w_[i] *= norm_d / norm_m;

  // read experimental overlaps
  vector<double> exp_err;
  if(errfile.size()>0) exp_err = read_exp_errors(errfile);

  // get self overlaps between data GMM components
  for(unsigned i=0; i<GMM_d_m_.size(); ++i) {
    double ov = get_self_overlap(i);
    ovdd_.push_back(ov);
  }

  log.printf("  number of GMM groups : %u\n", static_cast<unsigned>(GMM_d_grps_.size()));
  // cycle on GMM groups
  for(unsigned Gid=0; Gid<GMM_d_grps_.size(); ++Gid) {
    log.printf("    group %d\n", Gid);
    // calculate median overlap and experimental error
    vector<double> ovdd;
    vector<double> err;
    // cycle on the group members
    for(unsigned i=0; i<GMM_d_grps_[Gid].size(); ++i) {
      // GMM id
      int GMMid = GMM_d_grps_[Gid][i];
      // add to experimental error
      if(errfile.size()>0) err.push_back(exp_err[GMMid]);
      else                 err.push_back(0.);
      // add to GMM overlap
      ovdd.push_back(ovdd_[GMMid]);
    }
    // calculate median quantities
    std::sort(ovdd.begin(), ovdd.end());
    double ovdd_m = ovdd[ovdd.size()/2];
    std::sort(err.begin(), err.end());
    double err_m = err[err.size()/2];
    // print out statistics
    log.printf("     # of members : %u\n", GMM_d_grps_[Gid].size());
    log.printf("     median overlap : %lf\n", ovdd_m);
    log.printf("     median error : %lf\n", err_m);
    // set sigma mean
    sigma_mean_.push_back(sigma_mean * ovdd_m);
    // set dsigma
    dsigma_.push_back(dsigma * ovdd_m);
    // add minimum value of sigma for this group of GMMs
    sigma_min_.push_back(err_m);
    // set sigma max
    sigma_max_.push_back(10.0*ovdd_m + sigma_min_[Gid] + dsigma_[Gid]);
    // initialize sigma
    sigma_.push_back(std::max(sigma_min_[Gid], sigma_ini*ovdd_m));
  }

  // read status file if restarting
  if(getRestart()) read_status();

  // calculate auxiliary stuff
  calculate_useful_stuff(reso);

  // prepare data and derivative vectors
  ovmd_.resize(GMM_d_m_.size());
  atom_der_.resize(GMM_m_type_.size());

  // clear things that are no longer needed
  GMM_d_cov_.clear();

  // add components
  addComponentWithDerivatives("scoreb"); componentIsNotPeriodic("scoreb");
  addComponent("acc"); componentIsNotPeriodic("acc");

  if(nregres_>0) {
    addComponent("scale");  componentIsNotPeriodic("scale");
    addComponent("offset"); componentIsNotPeriodic("offset");
    addComponent("corr");   componentIsNotPeriodic("corr");
  }
  if(nanneal_>0) {addComponent("anneal"); componentIsNotPeriodic("anneal");}

  // initialize random seed
  unsigned iseed;
  if(rank_==0) iseed = time(NULL)+replica_;
  else iseed = 0;
  comm.Sum(&iseed, 1);
  random_.setSeed(-iseed);

  // request the atoms
  requestAtoms(atoms);

  log<<"  Bibliography "<<plumed.cite("Bonomi, Camilloni, Cavalli, Vendruscolo, Sci. Adv. 2, e150117 (2016)");
  log<<plumed.cite("Hanot, Bonomi, Greenberg, Sali, Nilges, Vendruscolo, Pellarin, bioRxiv doi: 10.1101/113951 (2017)");
  log<<plumed.cite("Bonomi, Camilloni, Bioinformatics, 33, 3999 (2017)");
  log<<"\n";
}

void EMMI::read_status()
{
  double MDtime;
// open file
  IFile *ifile = new IFile();
  ifile->link(*this);
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

bool EMMI::doAccept(double oldE, double newE) {
  bool accept = false;
  // calculate delta energy
  double delta = ( newE - oldE ) / kbt_;
  // if delta is negative always accept move
  if( delta < 0.0 ) {
    accept = true;
  } else {
    // otherwise extract random number
    double s = random_.RandU01();
    if( s < exp(-delta) ) { accept = true; }
  }
  return accept;
}

void EMMI::doMonteCarlo()
{
  // extract random GMM group
  unsigned nGMM = static_cast<unsigned>(floor(random_.RandU01()*static_cast<double>(GMM_d_grps_.size())));
  if(nGMM==GMM_d_grps_.size()) nGMM -= 1;

  // generate random move
  double shift = dsigma_[nGMM] * ( 2.0 * random_.RandU01() - 1.0 );
  // new sigma
  double new_s = sigma_[nGMM] + shift;
  // check boundaries
  if(new_s > sigma_max_[nGMM]) {new_s = 2.0 * sigma_max_[nGMM] - new_s;}
  if(new_s < sigma_min_[nGMM]) {new_s = 2.0 * sigma_min_[nGMM] - new_s;}
  // old s2
  double inv_old_s2 = 1.0 / ( sigma_mean_[nGMM]*sigma_mean_[nGMM]+sigma_[nGMM]*sigma_[nGMM] );
  // new s2
  double inv_new_s2 = 1.0 / ( sigma_mean_[nGMM]*sigma_mean_[nGMM]+new_s*new_s );

  // cycle on GMM group and calculate chi2
  double chi2 = 0.0;
  for(unsigned i=0; i<GMM_d_grps_[nGMM].size(); ++i) {
    // id GMM component
    int GMMid = GMM_d_grps_[nGMM][i];
    // deviation
    double dev = ( scale_*ovmd_[GMMid]+off_-ovdd_[GMMid] );
    // add to chi2
    chi2 += dev * dev;
  }
  // final energy calculation: add normalization and prior
  double ng = static_cast<double>(GMM_d_grps_[nGMM].size());
  double old_ene = 0.5 * kbt_ * ( chi2 * inv_old_s2 - (ng+1.0) * std::log(inv_old_s2) );
  double new_ene = 0.5 * kbt_ * ( chi2 * inv_new_s2 - (ng+1.0) * std::log(inv_new_s2) );

  // accept or reject
  bool accept = doAccept(old_ene/anneal_, new_ene/anneal_);
  if(accept) {
    sigma_[nGMM] = new_s;
    MCaccept_++;
  }
  // local communication
  if(rank_!=0) {
    for(unsigned i=0; i<sigma_.size(); ++i) sigma_[i] = 0.0;
    MCaccept_ = 0;
  }
  if(size_>1) {
    comm.Sum(&sigma_[0], sigma_.size());
    comm.Sum(&MCaccept_, 1);
  }
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
      // total experimental error
      double err_tot = 0.0;
      // cycle on number of experimental overlaps
      for(unsigned i=0; i<nexp; ++i) {
        string ss; Tools::convert(i,ss);
        ifile->scanField("Err"+ss, err);
        // add to total error
        err_tot += err*err;
      }
      // new line
      ifile->scanField();
      // calculate RMSE
      err_tot = sqrt(err_tot/static_cast<double>(nexp));
      // add to global
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
  if(w<=0) error("check data GMM: weight must be positive");
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
      if(beta<0) error("Beta must be positive!");
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
  // now create a set from beta (unique set of values)
  set<int> bu(GMM_d_beta_.begin(), GMM_d_beta_.end());
  // now prepare the group vector
  GMM_d_grps_.resize(bu.size());
  // and fill it in
  for(unsigned i=0; i<GMM_d_beta_.size(); ++i) {
    if(GMM_d_beta_[i]>=GMM_d_grps_.size()) error("Check Beta values");
    GMM_d_grps_[GMM_d_beta_[i]].push_back(i);
  }
}

void EMMI::calculate_useful_stuff(double reso)
{
  // calculate effective resolution
  double ave_s2 = 0.0;
  for(unsigned i=0; i<GMM_m_type_.size(); ++i) {
    // Gaussian sigma in real space (and in nm)
    double s = sqrt ( 0.5 * GMM_m_s_[GMM_m_type_[i]] ) / pi * 0.1;
    // add to average
    ave_s2 += s*s;
  }
  ave_s2 /= static_cast<double>(GMM_m_type_.size());
  // target effective sigma squared
  double st = reso*reso/4.0;
  // calculate blur factor
  double blur = 0.0;
  if(st>ave_s2) blur = sqrt(st-ave_s2);
  log.printf("  map resolution : %f\n", reso);
  log.printf("  Gaussian blur factor : %f\n", blur);

  // now calculate useful stuff
  VectorGeneric<6> cov, sum, inv_sum;
  // cycle on all atoms types (4 for the moment)
  for(unsigned i=0; i<GMM_m_s_.size(); ++i) {
    // the Gaussian in density (real) space is the FT of scattering factor
    // f(r) = A * (pi/B)**1.5 * exp(-pi**2/B*r**2)
    double s = sqrt ( 0.5 * GMM_m_s_[i] ) / pi * 0.1;
    // calculate s2 and add Gaussian blur
    double s2 = s*s + blur*blur;
    // covariance matrix for spherical Gaussian
    cov[0]=s2; cov[1]=0.0; cov[2]=0.0;
    cov[3]=s2; cov[4]=0.0;
    cov[5]=s2;
    // cycle on all data GMM
    for(unsigned j=0; j<GMM_d_m_.size(); ++j) {
      // we need the sum of the covariance matrices
      for(unsigned k=0; k<6; ++k) sum[k] = cov[k] + GMM_d_cov_[j][k];
      // and to calculate its determinant
      double det = sum[0]*(sum[3]*sum[5]-sum[4]*sum[4]);
      det -= sum[1]*(sum[1]*sum[5]-sum[4]*sum[2]);
      det += sum[2]*(sum[1]*sum[4]-sum[3]*sum[2]);
      // calculate prefactor - model weights are already normalized
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
  // tabulate exponential
  dexp_ = dpcutoff_ / static_cast<double> (nexp_-1);
  for(unsigned i=0; i<nexp_; ++i) {
    tab_exp_.push_back(exp(-static_cast<double>(i) * dexp_));
  }
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
// and return it
  return ov_tot;
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
  // dimension of GMM and atom vectors
  unsigned GMM_d_size = GMM_d_m_.size();
  unsigned GMM_m_size = GMM_m_type_.size();
  // local neighbor list
  vector < unsigned > nl_l;
  // clear old neighbor list
  nl_.clear();

  // cycle on GMM components - in parallel
  for(unsigned id=rank_; id<GMM_d_size; id+=size_) {
    // overlap lists and map
    vector<double> ov_l;
    map<double, unsigned> ov_m;
    // total overlap with id
    double ov_tot = 0.0;
    // cycle on all atoms
    for(unsigned im=0; im<GMM_m_size; ++im) {
      // get index in auxiliary lists
      unsigned kaux = GMM_m_type_[im] * GMM_d_size + id;
      // calculate exponent of overlap
      double expov = get_exp_overlap(GMM_d_m_[id], getPosition(im), inv_cov_md_[kaux]);
      // get index of 0.5*expov in tabulated exponential
      unsigned itab = static_cast<unsigned> (round( 0.5*expov/dexp_ ));
      // check boundaries and skip atom in case
      if(itab >= tab_exp_.size()) continue;
      // in case calculate overlap overlap
      double ov = pre_fact_[kaux] * tab_exp_[itab];
      // add to list
      ov_l.push_back(ov);
      // and map to retrieve atom index
      ov_m[ov] = im;
      // increase ov_tot
      ov_tot += ov;
    }
    // check if zero size -> ov_tot = 0
    if(ov_l.size()==0) continue;
    // define cutoff
    double ov_cut = ov_tot * nl_cutoff_;
    // sort ov_l in ascending order
    std::sort(ov_l.begin(), ov_l.end());
    // integrate ov_l
    double res = 0.0;
    for(unsigned i=0; i<ov_l.size(); ++i) {
      res += ov_l[i];
      // if exceeding the cutoff for overlap, stop
      if(res >= ov_cut) break;
      else ov_m.erase(ov_l[i]);
    }
    // now add atoms to neighborlist
    for(map<double, unsigned>::iterator it=ov_m.begin(); it!=ov_m.end(); ++it)
      nl_l.push_back(id*GMM_m_size+it->second);
    // end cycle on GMM components in parallel
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
  if(size_>1) {
    comm.Sum(&ovmd_[0], ovmd_.size());
    comm.Sum(&ovmd_der_[0][0], 3*ovmd_der_.size());
  }
}

void EMMI::doRegression()
{
// calculate averages
  double ovdd_ave = 0.0;
  double ovmd_ave = 0.0;
  for(unsigned i=0; i<ovdd_.size(); ++i) {
    ovdd_ave += ovdd_[i];
    ovmd_ave += ovmd_[i];
  }
  ovdd_ave /= static_cast<double>(ovdd_.size());
  ovmd_ave /= static_cast<double>(ovmd_.size());
// calculate scaling and offset
  double Bn = 0.0; double Bd = 0.0; double C = 0.0;
  for(unsigned i=0; i<ovmd_.size(); ++i) {
    // add to sum
    Bn += ( ovmd_[i] - ovmd_ave ) * ( ovdd_[i] - ovdd_ave );
    Bd += ( ovmd_[i] - ovmd_ave ) * ( ovmd_[i] - ovmd_ave );
    C  += ( ovdd_[i] - ovdd_ave ) * ( ovdd_[i] - ovdd_ave );
  }
// reset scale_
  if(Bd<=0 || Bn<=0) {
    scale_ = 1.;
    off_ = 0.;
    corr_ = 1.0;
  } else {
    scale_ = Bn / Bd;
    off_ = ovdd_ave - scale_*ovmd_ave;
    corr_ = Bn / sqrt(Bd) / sqrt(C);
  }
}

void EMMI::calculate()
{

// calculate CV
  calculate_overlap();

  if(!analysis_) {

    // BIASING MODE
    calculate_sigma();

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

double EMMI::get_annealing(long int step)
{
// default no annealing
  double fact = 1.0;
// position in annealing cycle
  unsigned nc = step%(4*nanneal_);
// useful doubles
  double ncd = static_cast<double>(nc);
  double nn  = static_cast<double>(nanneal_);
// set fact
  if(nc>=nanneal_   && nc<2*nanneal_) fact = (kanneal_-1.0) / nn * ( ncd - nn ) + 1.0;
  if(nc>=2*nanneal_ && nc<3*nanneal_) fact = kanneal_;
  if(nc>=3*nanneal_)                  fact = (1.0-kanneal_) / nn * ( ncd - 3.0*nn) + kanneal_;
  return fact;
}

void EMMI::calculate_sigma()
{
  // NOTE: all ranks and replicas know local ovmd and sigma at this point

  // rescale factor for ensemble average
  double escale = 1.0 / static_cast<double>(nrep_);

  // prepare vector of inv_s2
  vector<double> inv_s2(sigma_.size(), 0.0);
  for(unsigned i=0; i<sigma_.size(); ++i)
    inv_s2[i] = 1.0 / ( sigma_mean_[i]*sigma_mean_[i]+sigma_[i]*sigma_[i] );

  // in case of ensemble averaging
  if(!no_aver_ && nrep_>1) {
    // if master node, calculate average across replicas and sum inv_s2
    if(rank_==0) {
      multi_sim_comm.Sum(&ovmd_[0],  ovmd_.size());
      multi_sim_comm.Sum(&inv_s2[0], inv_s2.size());
      // divide model overlap by number of replicas
      for(unsigned i=0; i<ovmd_.size(); ++i) ovmd_[i] *= escale;
    } else {
      // set model overlaps and uncertainties to zero
      for(unsigned i=0; i<ovmd_.size(); ++i)  ovmd_[i] = 0.0;
      for(unsigned i=0; i<inv_s2.size(); ++i) inv_s2[i] = 0.0;
    }
    // local communication
    if(size_>1) {
      comm.Sum(&ovmd_[0],  ovmd_.size());
      comm.Sum(&inv_s2[0], inv_s2.size());
    }
    // NOTE: all ranks and replicas know average ovmd and global inv_s2 at this point
  }

  // get time step
  long int step = getStep();

  // do regression
  if(nregres_>0 && step%nregres_==0 && !getExchangeStep()) doRegression();

  // get annealing rescale factor
  if(nanneal_>0) {
    anneal_ = get_annealing(step);
    getPntrToComponent("anneal")->set(anneal_);
  }

  // calculate score and reweighting score
  double ene = 0.0;
  // cycle on GMM groups
  for(unsigned i=0; i<GMM_d_grps_.size(); ++i) {
    double chi2 = 0.0;
    // cycle on GMM members of the group
    for(unsigned j=0; j<GMM_d_grps_[i].size(); ++j) {
      // id of the GMM component
      int GMMid = GMM_d_grps_[i][j];
      // deviation
      double dev = ( scale_*ovmd_[GMMid]+off_-ovdd_[GMMid] );
      // add to chi2
      chi2 += dev * dev;
    }
    // add to energy
    ene += 0.5 * kbt_ * ( chi2 * inv_s2[i] - (static_cast<double>(GMM_d_grps_[i].size())+1.0) * std::log(inv_s2[i]) );
  }

  // annealing rescale
  ene /= anneal_;

  // clean temporary vector
  for(unsigned i=0; i<atom_der_.size(); ++i) atom_der_[i] = Vector(0,0,0);

  // virial
  Tensor virial;
  // get derivatives of bias with respect to atoms
  for(unsigned i=rank_; i<nl_.size(); i=i+size_) {
    // get indexes of data and model component
    unsigned id = nl_[i] / GMM_m_type_.size();
    unsigned im = nl_[i] % GMM_m_type_.size();
    // GMM group id
    int Gid = GMM_d_beta_[id];
    // derivative
    double der = kbt_ * ( scale_*ovmd_[id]+off_-ovdd_[id] ) * inv_s2[Gid];
    // chain rule + replica normalization
    Vector tot_der = der * ovmd_der_[i] * escale * scale_ / anneal_;
    Vector pos;
    if(pbc_) pos = pbcDistance(GMM_d_m_[id], getPosition(im)) + GMM_d_m_[id];
    else     pos = getPosition(im);
    // increment derivatives and virial
    atom_der_[im] += tot_der;
    virial += Tensor(pos, -tot_der);
  }

  // communicate local derivatives and virial
  if(size_>1) {
    comm.Sum(&atom_der_[0][0], 3*atom_der_.size());
    comm.Sum(virial);
  }

  // set derivatives, virial, and score
  for(unsigned i=0; i<atom_der_.size(); ++i) setAtomsDerivatives(getPntrToComponent("scoreb"), i, atom_der_[i]);
  setBoxDerivatives(getPntrToComponent("scoreb"), virial);
  getPntrToComponent("scoreb")->set(ene);

  // do Montecarlo
  if(dsigma_[0]>0 && step%MCstride_==0 && !getExchangeStep()) doMonteCarlo();

  // print status
  if(step%statusstride_==0) print_status(step);

  // calculate acceptance ratio
  // this is needed when restarting simulations
  if(MCfirst_==-1) MCfirst_=step;
  // acceptance
  double MCtrials = std::floor(static_cast<double>(step-MCfirst_) / static_cast<double>(MCstride_))+1.0;
  double acc = static_cast<double>(MCaccept_) / MCtrials;
  // set value
  getPntrToComponent("acc")->set(acc);

  // print scale
  if(nregres_>0) {
    getPntrToComponent("scale")->set(scale_);
    getPntrToComponent("offset")->set(off_);
    getPntrToComponent("corr")->set(corr_);
  }
}

}
}
