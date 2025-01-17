/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2023 The plumed team
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
#include "core/ActionWithValue.h"
#include "core/ActionAtomistic.h"
#include "core/ActionWithArguments.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Matrix.h"
#include "core/GenericMolInfo.h"
#include "core/ActionSet.h"
#include "tools/File.h"
#include "tools/Random.h"

#include <string>
#include <map>
#include <numeric>
#include <ctime>

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
Combined with a multi-replica framework (such as the -multi option in GROMACS), the user can model an ensemble of structures using
the Metainference approach \cite Bonomi:2016ip .

\warning
To use \ref EMMI, the user should always add a \ref MOLINFO line and specify a pdb file of the system.

\note
To enhance sampling in single-structure refinement, one can use a Replica Exchange Method, such as Parallel Tempering.
In this case, the user should add the NO_AVER flag to the input line. To use a replica-based enhanced sampling scheme such as
Parallel-Bias Metadynamics (\ref PBMETAD), one should use the REWEIGHT flag and pass the Metadynamics bias using the ARG keyword.

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
gmm: EMMI NOPBC SIGMA_MIN=0.01 TEMP=300.0 NL_STRIDE=100 NL_CUTOFF=0.01 GMM_FILE=GMM_fit.dat ATOMS=protein-h

# translate into bias - apply every 2 steps
emr: BIASVALUE ARG=gmm.scoreb STRIDE=2

PRINT ARG=emr.* FILE=COLVAR STRIDE=500 FMT=%20.10f
\endplumedfile


*/
//+ENDPLUMEDOC

class EMMI :
  public ActionAtomistic,
  public ActionWithArguments,
  public ActionWithValue
{
private:

// temperature in kbt
  double kbt_;
// model GMM - atom types
  std::vector<unsigned> GMM_m_type_;
// model GMM - list of atom sigmas - one per atom type
  std::vector<double> GMM_m_s_;
// model GMM - list of atom weights - one per atom type
  std::vector<double> GMM_m_w_;
// data GMM - means, weights, and covariances + beta option
  std::vector<Vector>             GMM_d_m_;
  std::vector<double>             GMM_d_w_;
  std::vector< VectorGeneric<6> > GMM_d_cov_;
  std::vector<int>                GMM_d_beta_;
  std::vector < std::vector<int> >     GMM_d_grps_;
// overlaps
  std::vector<double> ovmd_;
  std::vector<double> ovdd_;
  std::vector<double> ovmd_ave_;
// and derivatives
  std::vector<Vector> ovmd_der_;
  std::vector<Vector> atom_der_;
  std::vector<double> GMMid_der_;
// constants
  double cfact_;
  double inv_sqrt2_, sqrt2_pi_;
// metainference
  unsigned nrep_;
  unsigned replica_;
  std::vector<double> sigma_;
  std::vector<double> sigma_min_;
  std::vector<double> sigma_max_;
  std::vector<double> dsigma_;
// list of prefactors for overlap between two Gaussians
// pre_fact = 1.0 / (2pi)**1.5 / std::sqrt(det_md) * Wm * Wd
  std::vector<double> pre_fact_;
// inverse of the sum of model and data covariances matrices
  std::vector< VectorGeneric<6> > inv_cov_md_;
// neighbor list
  double   nl_cutoff_;
  unsigned nl_stride_;
  bool first_time_;
  bool no_aver_;
  std::vector<unsigned> nl_;
// parallel stuff
  unsigned size_;
  unsigned rank_;
// pbc
  bool pbc_;
// Monte Carlo stuff
  int      MCstride_;
  double   MCaccept_;
  double   MCtrials_;
  Random   random_;
  // status stuff
  unsigned int statusstride_;
  std::string       statusfilename_;
  OFile        statusfile_;
  bool         first_status_;
  // regression
  unsigned nregres_;
  double scale_;
  double scale_min_;
  double scale_max_;
  double dscale_;
  // tabulated exponential
  double dpcutoff_;
  double dexp_;
  unsigned nexp_;
  std::vector<double> tab_exp_;
  // simulated annealing
  unsigned nanneal_;
  double   kanneal_;
  double   anneal_;
  // prior exponent
  double prior_;
  // noise type
  unsigned noise_;
  // total score and virial;
  double ene_;
  Tensor virial_;
  // model overlap file
  unsigned int ovstride_;
  std::string       ovfilename_;

  // Reweighting additions
  bool do_reweight_;
  bool first_time_w_;
  std::vector<double> forces;
  std::vector<double> forcesToApply;

  // average weights
  double decay_w_;
  std::vector<double> average_weights_;

// write file with model overlap
  void write_model_overlap(long long int step);
// get median of std::vector
  double get_median(std::vector<double> &v);
// annealing
  double get_annealing(long long int step);
// do regression
  double scaleEnergy(double s);
  double doRegression();
// read and write status
  void read_status();
  void print_status(long long int step);
// accept or reject
  bool doAccept(double oldE, double newE, double kbt);
// do MonteCarlo
  void doMonteCarlo();
// read error file
  std::vector<double> read_exp_errors(const std::string & errfile);
// read experimental overlaps
  std::vector<double> read_exp_overlaps(const std::string & ovfile);
// calculate model GMM weights and covariances
  std::vector<double> get_GMM_m(std::vector<AtomNumber> &atoms);
// read data GMM file
  void get_GMM_d(std::string gmm_file);
// check GMM data
  void check_GMM_d(const VectorGeneric<6> &cov, double w);
// auxiliary method
  void calculate_useful_stuff(double reso);
// get pref_fact and inv_cov_md
  double get_prefactor_inverse (const VectorGeneric<6> &GMM_cov_0, const VectorGeneric<6> &GMM_cov_1,
                                double GMM_w_0, double GMM_w_1,
                                VectorGeneric<6> &sum, VectorGeneric<6> &inv_sum);
// calculate self overlaps between data GMM components - ovdd_
  double get_self_overlap(unsigned id);
// calculate overlap between two Gaussians
  double get_overlap(const Vector &m_m, const Vector &d_m, double pre_fact,
                     const VectorGeneric<6> &inv_cov_md, Vector &ov_der);
// calculate exponent of overlap for neighbor list update
  double get_exp_overlap(const Vector &m_m, const Vector &d_m,
                         const VectorGeneric<6> &inv_cov_md);
// update the neighbor list
  void update_neighbor_list();
// calculate overlap
  void calculate_overlap();
// Gaussian noise
  void calculate_Gauss();
// Outliers noise
  void calculate_Outliers();
// Marginal noise
  void calculate_Marginal();

  // See MetainferenceBase
  void get_weights(double &weight, double &norm, double &neff);

public:
  static void registerKeywords( Keywords& keys );
  explicit EMMI(const ActionOptions&);
  // needed for reweighting
  void setDerivatives();
  void turnOnDerivatives() override;
  unsigned getNumberOfDerivatives() override;
  void lockRequests() override;
  void unlockRequests() override;
  void calculateNumericalDerivatives( ActionWithValue* a ) override;
  void apply() override;
  void setArgDerivatives(Value *v, const double &d);
  void setAtomsDerivatives(Value*v, const unsigned i, const Vector&d);
  void setBoxDerivatives(Value*v, const Tensor&d);
// active methods:
  void prepare() override;
  void calculate() override;
};

inline
void EMMI::setDerivatives() {
  // Get appropriate number of derivatives
  // Derivatives are first for arguments and then for atoms
  unsigned nder;
  if( getNumberOfAtoms()>0 ) {
    nder = 3*getNumberOfAtoms() + 9 + getNumberOfArguments();
  } else {
    nder = getNumberOfArguments();
  }

  // Resize all derivative arrays
  forces.resize( nder ); forcesToApply.resize( nder );
  for(int i=0; i<getNumberOfComponents(); ++i) getPntrToComponent(i)->resizeDerivatives(nder);
}

inline
void EMMI::turnOnDerivatives() {
  ActionWithValue::turnOnDerivatives();
}

inline
unsigned EMMI::getNumberOfDerivatives() {
  if( getNumberOfAtoms()>0 ) {
    return 3*getNumberOfAtoms() + 9 + getNumberOfArguments();
  }
  return getNumberOfArguments();
}

inline
void EMMI::lockRequests() {
  ActionAtomistic::lockRequests();
  ActionWithArguments::lockRequests();
}

inline
void EMMI::unlockRequests() {
  ActionAtomistic::unlockRequests();
  ActionWithArguments::unlockRequests();
}

inline
void EMMI::calculateNumericalDerivatives( ActionWithValue* a=NULL ) {
  if( getNumberOfArguments()>0 ) {
    ActionWithArguments::calculateNumericalDerivatives( a );
  }
  if( getNumberOfAtoms()>0 ) {
    Matrix<double> save_derivatives( getNumberOfComponents(), getNumberOfArguments() );
    for(int j=0; j<getNumberOfComponents(); ++j) {
      for(unsigned i=0; i<getNumberOfArguments(); ++i) if(getPntrToComponent(j)->hasDerivatives()) save_derivatives(j,i)=getPntrToComponent(j)->getDerivative(i);
    }
    calculateAtomicNumericalDerivatives( a, getNumberOfArguments() );
    for(int j=0; j<getNumberOfComponents(); ++j) {
      for(unsigned i=0; i<getNumberOfArguments(); ++i) if(getPntrToComponent(j)->hasDerivatives()) getPntrToComponent(j)->addDerivative( i, save_derivatives(j,i) );
    }
  }
}

inline
void EMMI::apply() {
  bool wasforced=false; forcesToApply.assign(forcesToApply.size(),0.0);
  for(int i=0; i<getNumberOfComponents(); ++i) {
    if( getPntrToComponent(i)->applyForce( forces ) ) {
      wasforced=true;
      for(unsigned i=0; i<forces.size(); ++i) forcesToApply[i]+=forces[i];
    }
  }
  if( wasforced ) {
    addForcesOnArguments( forcesToApply );
    if( getNumberOfAtoms()>0 ) setForcesOnAtoms( forcesToApply, getNumberOfArguments() );
  }
}

inline
void EMMI::setArgDerivatives(Value *v, const double &d) {
  v->addDerivative(0,d);
}

inline
void EMMI::setAtomsDerivatives(Value*v, const unsigned i, const Vector&d) {
  const unsigned noa=getNumberOfArguments();
  v->addDerivative(noa+3*i+0,d[0]);
  v->addDerivative(noa+3*i+1,d[1]);
  v->addDerivative(noa+3*i+2,d[2]);
}

inline
void EMMI::setBoxDerivatives(Value* v,const Tensor&d) {
  const unsigned noa=getNumberOfArguments();
  const unsigned nat=getNumberOfAtoms();
  v->addDerivative(noa+3*nat+0,d(0,0));
  v->addDerivative(noa+3*nat+1,d(0,1));
  v->addDerivative(noa+3*nat+2,d(0,2));
  v->addDerivative(noa+3*nat+3,d(1,0));
  v->addDerivative(noa+3*nat+4,d(1,1));
  v->addDerivative(noa+3*nat+5,d(1,2));
  v->addDerivative(noa+3*nat+6,d(2,0));
  v->addDerivative(noa+3*nat+7,d(2,1));
  v->addDerivative(noa+3*nat+8,d(2,2));
}

PLUMED_REGISTER_ACTION(EMMI,"EMMI")

void EMMI::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  keys.use("ARG");
  keys.add("atoms","ATOMS","atoms for which we calculate the density map, typically all heavy atoms");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.add("compulsory","GMM_FILE","file with the parameters of the GMM components");
  keys.add("compulsory","NL_CUTOFF","The cutoff in overlap for the neighbor list");
  keys.add("compulsory","NL_STRIDE","The frequency with which we are updating the neighbor list");
  keys.add("compulsory","SIGMA_MIN","minimum uncertainty");
  keys.add("compulsory","RESOLUTION", "Cryo-EM map resolution");
  keys.add("compulsory","NOISETYPE","functional form of the noise (GAUSS, OUTLIERS, MARGINAL)");
  keys.add("optional","SIGMA0","initial value of the uncertainty");
  keys.add("optional","DSIGMA","MC step for uncertainties");
  keys.add("optional","MC_STRIDE", "Monte Carlo stride");
  keys.add("optional","ERR_FILE","file with experimental or GMM fit errors");
  keys.add("optional","OV_FILE","file with experimental overlaps");
  keys.add("optional","NORM_DENSITY","integral of the experimental density");
  keys.add("optional","STATUS_FILE","write a file with all the data useful for restart");
  keys.add("optional","WRITE_STRIDE","write the status to a file every N steps, this can be used for restart");
  keys.add("optional","REGRESSION","regression stride");
  keys.add("optional","REG_SCALE_MIN","regression minimum scale");
  keys.add("optional","REG_SCALE_MAX","regression maximum scale");
  keys.add("optional","REG_DSCALE","regression maximum scale MC move");
  keys.add("optional","SCALE","scale factor");
  keys.add("optional","ANNEAL", "Length of annealing cycle");
  keys.add("optional","ANNEAL_FACT", "Annealing temperature factor");
  keys.add("optional","TEMP","temperature");
  keys.add("optional","PRIOR", "exponent of uncertainty prior");
  keys.add("optional","WRITE_OV_STRIDE","write model overlaps every N steps");
  keys.add("optional","WRITE_OV","write a file with model overlaps");
  keys.add("optional","AVERAGING", "Averaging window for weights");
  keys.addFlag("NO_AVER",false,"don't do ensemble averaging in multi-replica mode");
  keys.addFlag("REWEIGHT",false,"simple REWEIGHT using the ARG as energy");
  componentsAreNotOptional(keys);
  keys.addOutputComponent("scoreb","default","Bayesian score");
  keys.addOutputComponent("acc",   "NOISETYPE","MC acceptance for uncertainty");
  keys.addOutputComponent("scale", "REGRESSION","scale factor");
  keys.addOutputComponent("accscale", "REGRESSION","MC acceptance for scale regression");
  keys.addOutputComponent("enescale", "REGRESSION","MC energy for scale regression");
  keys.addOutputComponent("anneal","ANNEAL","annealing factor");
  keys.addOutputComponent("weight",       "REWEIGHT",     "weights of the weighted average");
  keys.addOutputComponent("biasDer",      "REWEIGHT",     "derivatives with respect to the bias");
  keys.addOutputComponent("sigma",      "NOISETYPE",     "uncertainty in the forward models and experiment");
  keys.addOutputComponent("neff",         "default",      "effective number of replicas");
}

EMMI::EMMI(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionWithArguments(ao),
  ActionWithValue(ao),
  inv_sqrt2_(0.707106781186548),
  sqrt2_pi_(0.797884560802865),
  first_time_(true), no_aver_(false), pbc_(true),
  MCstride_(1), MCaccept_(0.), MCtrials_(0.),
  statusstride_(0), first_status_(true),
  nregres_(0), scale_(1.),
  dpcutoff_(15.0), nexp_(1000000), nanneal_(0),
  kanneal_(0.), anneal_(1.), prior_(1.), ovstride_(0),
  do_reweight_(false), first_time_w_(true), decay_w_(1.)
{
  // periodic boundary conditions
  bool nopbc=!pbc_;
  parseFlag("NOPBC",nopbc);
  pbc_=!nopbc;

  // list of atoms
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS", atoms);

  // file with data GMM
  std::string GMM_file;
  parse("GMM_FILE", GMM_file);

  // type of data noise
  std::string noise;
  parse("NOISETYPE",noise);
  if      (noise=="GAUSS")   noise_ = 0;
  else if(noise=="OUTLIERS") noise_ = 1;
  else if(noise=="MARGINAL") noise_ = 2;
  else error("Unknown noise type!");

  // minimum value for error
  double sigma_min;
  parse("SIGMA_MIN", sigma_min);
  if(sigma_min<0) error("SIGMA_MIN should be greater or equal to zero");

  // the following parameters must be specified with noise type 0 and 1
  double sigma_ini, dsigma;
  if(noise_!=2) {
    // initial value of the uncertainty
    parse("SIGMA0", sigma_ini);
    if(sigma_ini<=0) error("you must specify a positive SIGMA0");
    // MC parameters
    parse("DSIGMA", dsigma);
    if(dsigma<0) error("you must specify a positive DSIGMA");
    parse("MC_STRIDE", MCstride_);
    if(dsigma>0 && MCstride_<=0) error("you must specify a positive MC_STRIDE");
    // status file parameters
    parse("WRITE_STRIDE", statusstride_);
    if(statusstride_==0) error("you must specify a positive WRITE_STRIDE");
    parse("STATUS_FILE",  statusfilename_);
    if(statusfilename_=="") statusfilename_ = "MISTATUS"+getLabel();
    else                    statusfilename_ = statusfilename_+getLabel();
  }

  // error file
  std::string errfile;
  parse("ERR_FILE", errfile);

  // file with experimental overlaps
  std::string ovfile;
  parse("OV_FILE", ovfile);

  // integral of the experimetal density
  double norm_d = 0.0;
  parse("NORM_DENSITY", norm_d);

  // temperature
  double temp=0.0;
  parse("TEMP",temp);
  // convert temp to kbt
  if(temp>0.0) kbt_=plumed.getAtoms().getKBoltzmann()*temp;
  else kbt_=plumed.getAtoms().getKbT();

  // exponent of uncertainty prior
  parse("PRIOR",prior_);

  // simulated annealing stuff
  parse("ANNEAL", nanneal_);
  parse("ANNEAL_FACT", kanneal_);
  if(nanneal_>0 && kanneal_<=1.0) error("with ANNEAL, ANNEAL_FACT must be greater than 1");

  // regression stride
  parse("REGRESSION",nregres_);
  // other regression parameters
  if(nregres_>0) {
    parse("REG_SCALE_MIN",scale_min_);
    parse("REG_SCALE_MAX",scale_max_);
    parse("REG_DSCALE",dscale_);
    // checks
    if(scale_max_<=scale_min_) error("with REGRESSION, REG_SCALE_MAX must be greater than REG_SCALE_MIN");
    if(dscale_<=0.) error("with REGRESSION, REG_DSCALE must be positive");
  }

  // scale factor
  parse("SCALE", scale_);

  // read map resolution
  double reso;
  parse("RESOLUTION", reso);
  if(reso<=0.) error("RESOLUTION should be strictly positive");

  // neighbor list stuff
  parse("NL_CUTOFF",nl_cutoff_);
  if(nl_cutoff_<=0.0) error("NL_CUTOFF should be explicitly specified and positive");
  parse("NL_STRIDE",nl_stride_);
  if(nl_stride_==0) error("NL_STRIDE should be explicitly specified and positive");

  // averaging or not
  parseFlag("NO_AVER",no_aver_);

  // write overlap file
  parse("WRITE_OV_STRIDE", ovstride_);
  parse("WRITE_OV", ovfilename_);
  if(ovstride_>0 && ovfilename_=="") error("With WRITE_OV_STRIDE you must specify WRITE_OV");

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

  // Reweighting flag
  parseFlag("REWEIGHT", do_reweight_);
  if(do_reweight_&&getNumberOfArguments()!=1) error("To REWEIGHT one must provide one single bias as an argument");
  if(do_reweight_&&no_aver_) error("REWEIGHT cannot be used with NO_AVER");
  if(do_reweight_&&nrep_<2) error("REWEIGHT can only be used in parallel with 2 or more replicas");
  if(!getRestart()) average_weights_.resize(nrep_, 1./static_cast<double>(nrep_));
  else average_weights_.resize(nrep_, 0.);

  unsigned averaging=0;
  parse("AVERAGING", averaging);
  if(averaging>0) {
    decay_w_ = 1./static_cast<double> (averaging);
  }

  checkRead();

  log.printf("  atoms involved : ");
  for(unsigned i=0; i<atoms.size(); ++i) log.printf("%d ",atoms[i].serial());
  log.printf("\n");
  log.printf("  GMM data file : %s\n", GMM_file.c_str());
  if(no_aver_) log.printf("  without ensemble averaging\n");
  log.printf("  type of data noise : %s\n", noise.c_str());
  log.printf("  neighbor list cutoff : %lf\n", nl_cutoff_);
  log.printf("  neighbor list stride : %u\n",  nl_stride_);
  log.printf("  minimum uncertainty : %f\n",sigma_min);
  log.printf("  scale factor : %lf\n",scale_);
  if(nregres_>0) {
    log.printf("  regression stride : %u\n", nregres_);
    log.printf("  regression minimum scale : %lf\n", scale_min_);
    log.printf("  regression maximum scale : %lf\n", scale_max_);
    log.printf("  regression maximum scale MC move : %lf\n", dscale_);
  }
  if(noise_!=2) {
    log.printf("  initial value of the uncertainty : %f\n",sigma_ini);
    log.printf("  max MC move in uncertainty : %f\n",dsigma);
    log.printf("  MC stride : %u\n", MCstride_);
    log.printf("  reading/writing to status file : %s\n",statusfilename_.c_str());
    log.printf("  with stride : %u\n",statusstride_);
  }
  if(errfile.size()>0) log.printf("  reading experimental errors from file : %s\n", errfile.c_str());
  if(ovfile.size()>0)  log.printf("  reading experimental overlaps from file : %s\n", ovfile.c_str());
  log.printf("  temperature of the system in energy unit : %f\n",kbt_);
  log.printf("  prior exponent : %f\n",prior_);
  log.printf("  number of replicas for averaging: %u\n",nrep_);
  log.printf("  id of the replica : %u\n",replica_);
  if(nanneal_>0) {
    log.printf("  length of annealing cycle : %u\n",nanneal_);
    log.printf("  annealing factor : %f\n",kanneal_);
  }
  if(ovstride_>0) {
    log.printf("  stride for writing model overlaps : %u\n",ovstride_);
    log.printf("  file for writing model overlaps : %s\n", ovfilename_.c_str());
  }

  // set constant quantity before calculating stuff
  cfact_ = 1.0/pow( 2.0*pi, 1.5 );

  // calculate model GMM constant parameters
  std::vector<double> GMM_m_w = get_GMM_m(atoms);

  // read data GMM parameters
  get_GMM_d(GMM_file);
  log.printf("  number of GMM components : %u\n", static_cast<unsigned>(GMM_d_m_.size()));

  // normalize atom weight map
  if(norm_d <= 0.0) norm_d = accumulate(GMM_d_w_.begin(), GMM_d_w_.end(), 0.0);
  double norm_m = accumulate(GMM_m_w.begin(),  GMM_m_w.end(),  0.0);
  // renormalization
  for(unsigned i=0; i<GMM_m_w_.size(); ++i) GMM_m_w_[i] *= norm_d / norm_m;

  // read experimental errors
  std::vector<double> exp_err;
  if(errfile.size()>0) exp_err = read_exp_errors(errfile);

  // get self overlaps between data GMM components
  if(ovfile.size()>0) {
    ovdd_ = read_exp_overlaps(ovfile);
  } else {
    for(unsigned i=0; i<GMM_d_m_.size(); ++i) {
      double ov = get_self_overlap(i);
      ovdd_.push_back(ov);
    }
  }

  log.printf("  number of GMM groups : %u\n", static_cast<unsigned>(GMM_d_grps_.size()));
  // cycle on GMM groups
  for(unsigned Gid=0; Gid<GMM_d_grps_.size(); ++Gid) {
    log.printf("    group %d\n", Gid);
    // calculate median overlap and experimental error
    std::vector<double> ovdd;
    std::vector<double> err;
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
    double ovdd_m = get_median(ovdd);
    double err_m  = get_median(err);
    // print out statistics
    log.printf("     # of members : %zu\n", GMM_d_grps_[Gid].size());
    log.printf("     median overlap : %lf\n", ovdd_m);
    log.printf("     median error : %lf\n", err_m);
    // add minimum value of sigma for this group of GMMs
    sigma_min_.push_back(std::sqrt(err_m*err_m+sigma_min*ovdd_m*sigma_min*ovdd_m));
    // these are only needed with Gaussian and Outliers noise models
    if(noise_!=2) {
      // set dsigma
      dsigma_.push_back(dsigma * ovdd_m);
      // set sigma max
      sigma_max_.push_back(10.0*ovdd_m + sigma_min_[Gid] + dsigma_[Gid]);
      // initialize sigma
      sigma_.push_back(std::max(sigma_min_[Gid],std::min(sigma_ini*ovdd_m,sigma_max_[Gid])));
    }
  }

  // read status file if restarting
  if(getRestart() && noise_!=2) read_status();

  // calculate auxiliary stuff
  calculate_useful_stuff(reso);

  // prepare data and derivative std::vectors
  ovmd_.resize(ovdd_.size());
  ovmd_ave_.resize(ovdd_.size());
  atom_der_.resize(GMM_m_type_.size());
  GMMid_der_.resize(ovdd_.size());

  // clear things that are no longer needed
  GMM_d_cov_.clear();

  // add components
  addComponentWithDerivatives("scoreb"); componentIsNotPeriodic("scoreb");

  if(noise_!=2) {addComponent("acc"); componentIsNotPeriodic("acc");}

  if(nregres_>0) {
    addComponent("scale");     componentIsNotPeriodic("scale");
    addComponent("accscale");  componentIsNotPeriodic("accscale");
    addComponent("enescale");  componentIsNotPeriodic("enescale");
  }

  if(nanneal_>0) {addComponent("anneal"); componentIsNotPeriodic("anneal");}

  if(do_reweight_) {
    addComponent("biasDer");
    componentIsNotPeriodic("biasDer");
    addComponent("weight");
    componentIsNotPeriodic("weight");
  }

  addComponent("neff");
  componentIsNotPeriodic("neff");

  for(unsigned i=0; i<sigma_.size(); ++i) {
    std::string num; Tools::convert(i, num);
    addComponent("sigma-"+num); componentIsNotPeriodic("sigma-"+num);
    getPntrToComponent("sigma-"+num)->set(sigma_[i]);
  }

  // initialize random seed
  unsigned iseed;
  if(rank_==0) iseed = time(NULL)+replica_;
  else iseed = 0;
  comm.Sum(&iseed, 1);
  random_.setSeed(-iseed);

  // request the atoms
  requestAtoms(atoms);
  setDerivatives();

  // print bibliography
  log<<"  Bibliography "<<plumed.cite("Bonomi, Camilloni, Bioinformatics, 33, 3999 (2017)");
  log<<plumed.cite("Hanot, Bonomi, Greenberg, Sali, Nilges, Vendruscolo, Pellarin, bioRxiv doi: 10.1101/113951 (2017)");
  log<<plumed.cite("Bonomi, Pellarin, Vendruscolo, Biophys. J. 114, 1604 (2018)");
  if(do_reweight_) log<<plumed.cite("Bonomi, Camilloni, Vendruscolo, Sci. Rep. 6, 31232 (2016)");
  if(!no_aver_ && nrep_>1)log<<plumed.cite("Bonomi, Camilloni, Cavalli, Vendruscolo, Sci. Adv. 2, e150117 (2016)");
  log<<"\n";
}

void EMMI::write_model_overlap(long long int step)
{
  OFile ovfile;
  ovfile.link(*this);
  std::string num; Tools::convert(step,num);
  std::string name = ovfilename_+"-"+num;
  ovfile.open(name);
  ovfile.setHeavyFlush();
  ovfile.fmtField("%10.7e ");
// write overlaps
  for(int i=0; i<ovmd_.size(); ++i) {
    ovfile.printField("Model", ovmd_[i]);
    ovfile.printField("ModelScaled", scale_ * ovmd_[i]);
    ovfile.printField("Data", ovdd_[i]);
    ovfile.printField();
  }
  ovfile.close();
}

double EMMI::get_median(std::vector<double> &v)
{
// dimension of std::vector
  unsigned size = v.size();
// in case of only one entry
  if (size==1) {
    return v[0];
  } else {
    // reorder std::vector
    sort(v.begin(), v.end());
    // odd or even?
    if (size%2==0) {
      return (v[size/2-1]+v[size/2])/2.0;
    } else {
      return v[size/2];
    }
  }
}

void EMMI::read_status()
{
  double MDtime;
// open file
  auto ifile = Tools::make_unique<IFile>();
  ifile->link(*this);
  if(ifile->FileExist(statusfilename_)) {
    ifile->open(statusfilename_);
    while(ifile->scanField("MD_time", MDtime)) {
      for(unsigned i=0; i<sigma_.size(); ++i) {
        // convert i to std::string
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
}

void EMMI::print_status(long long int step)
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
    // convert i to std::string
    std::string num; Tools::convert(i,num);
    // print entry
    statusfile_.printField("s"+num, sigma_[i]);
  }
  statusfile_.printField();
}

bool EMMI::doAccept(double oldE, double newE, double kbt) {
  bool accept = false;
  // calculate delta energy
  double delta = ( newE - oldE ) / kbt;
  // if delta is negative always accept move
  if( delta < 0.0 ) {
    accept = true;
  } else {
    // otherwise extract random number
    double s = random_.RandU01();
    if( s < std::exp(-delta) ) { accept = true; }
  }
  return accept;
}

void EMMI::doMonteCarlo()
{
  // extract random GMM group
  unsigned nGMM = static_cast<unsigned>(std::floor(random_.RandU01()*static_cast<double>(GMM_d_grps_.size())));
  if(nGMM==GMM_d_grps_.size()) nGMM -= 1;

  // generate random move
  double shift = dsigma_[nGMM] * ( 2.0 * random_.RandU01() - 1.0 );
  // new sigma
  double new_s = sigma_[nGMM] + shift;
  // check boundaries
  if(new_s > sigma_max_[nGMM]) {new_s = 2.0 * sigma_max_[nGMM] - new_s;}
  if(new_s < sigma_min_[nGMM]) {new_s = 2.0 * sigma_min_[nGMM] - new_s;}
  // old s2
  double old_inv_s2 = 1.0 / sigma_[nGMM] / sigma_[nGMM];
  // new s2
  double new_inv_s2 = 1.0 / new_s / new_s;

  // cycle on GMM group and calculate old and new energy
  double old_ene = 0.0;
  double new_ene = 0.0;
  double ng = static_cast<double>(GMM_d_grps_[nGMM].size());

  // in case of Gaussian noise
  if(noise_==0) {
    double chi2 = 0.0;
    for(unsigned i=0; i<GMM_d_grps_[nGMM].size(); ++i) {
      // id GMM component
      int GMMid = GMM_d_grps_[nGMM][i];
      // deviation
      double dev = ( scale_*ovmd_ave_[GMMid]-ovdd_[GMMid] );
      // add to chi2
      chi2 += dev * dev;
    }
    // final energy calculation: add normalization and prior
    old_ene = 0.5 * kbt_ * ( chi2 * old_inv_s2 - (ng+prior_) * std::log(old_inv_s2) );
    new_ene = 0.5 * kbt_ * ( chi2 * new_inv_s2 - (ng+prior_) * std::log(new_inv_s2) );
  }

  // in case of Outliers noise
  if(noise_==1) {
    for(unsigned i=0; i<GMM_d_grps_[nGMM].size(); ++i) {
      // id GMM component
      int GMMid = GMM_d_grps_[nGMM][i];
      // calculate deviation
      double dev = ( scale_*ovmd_ave_[GMMid]-ovdd_[GMMid] );
      // add to energies
      old_ene += std::log1p( 0.5 * dev * dev * old_inv_s2);
      new_ene += std::log1p( 0.5 * dev * dev * new_inv_s2);
    }
    // final energy calculation: add normalization and prior
    old_ene = kbt_ * ( old_ene + (ng+prior_) * std::log(sigma_[nGMM]) );
    new_ene = kbt_ * ( new_ene + (ng+prior_) * std::log(new_s) );
  }

  // increment number of trials
  MCtrials_ += 1.0;

  // accept or reject
  bool accept = doAccept(old_ene/anneal_, new_ene/anneal_, kbt_);
  if(accept) {
    sigma_[nGMM] = new_s;
    MCaccept_ += 1.0;
  }
  // local communication
  if(rank_!=0) {
    for(unsigned i=0; i<sigma_.size(); ++i) sigma_[i] = 0.0;
    MCaccept_ = 0.0;
  }
  if(size_>1) {
    comm.Sum(&sigma_[0], sigma_.size());
    comm.Sum(&MCaccept_, 1);
  }

  // update sigma output
  std::string num; Tools::convert(nGMM, num);
  getPntrToComponent("sigma-"+num)->set(sigma_[nGMM]);
}

std::vector<double> EMMI::read_exp_errors(const std::string & errfile)
{
  int nexp, idcomp;
  double err;
  std::vector<double> exp_err;
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
        std::string ss; Tools::convert(i,ss);
        ifile->scanField("Err"+ss, err);
        // add to total error
        err_tot += err*err;
      }
      // new line
      ifile->scanField();
      // calculate RMSE
      err_tot = std::sqrt(err_tot/static_cast<double>(nexp));
      // add to global
      exp_err.push_back(err_tot);
    }
    ifile->close();
  } else {
    error("Cannot find ERR_FILE "+errfile+"\n");
  }
  return exp_err;
}

std::vector<double> EMMI::read_exp_overlaps(const std::string & ovfile)
{
  int idcomp;
  double ov;
  std::vector<double> ovdd;
// open file
  IFile *ifile = new IFile();
  if(ifile->FileExist(ovfile)) {
    ifile->open(ovfile);
    // cycle on GMM components
    while(ifile->scanField("Id",idcomp)) {
      // read experimental overlap
      ifile->scanField("Overlap", ov);
      // add to ovdd
      ovdd.push_back(ov);
      // new line
      ifile->scanField();
    }
    ifile->close();
  } else {
    error("Cannot find OV_FILE "+ovfile+"\n");
  }
  return ovdd;
}

std::vector<double> EMMI::get_GMM_m(std::vector<AtomNumber> &atoms)
{
  // list of weights - one per atom
  std::vector<double> GMM_m_w;

  auto* moldat=plumed.getActionSet().selectLatest<GenericMolInfo*>(this);
  // map of atom types to A and B coefficients of scattering factor
  // f(s) = A * exp(-B*s**2)
  // B is in Angstrom squared
  // map between an atom type and an index
  std::map<std::string, unsigned> type_map;
  type_map["C"]=0;
  type_map["O"]=1;
  type_map["N"]=2;
  type_map["S"]=3;
  // fill in sigma std::vector
  GMM_m_s_.push_back(0.01*15.146);  // type 0
  GMM_m_s_.push_back(0.01*8.59722); // type 1
  GMM_m_s_.push_back(0.01*11.1116); // type 2
  GMM_m_s_.push_back(0.01*15.8952); // type 3
  // fill in weight std::vector
  GMM_m_w_.push_back(2.49982); // type 0
  GMM_m_w_.push_back(1.97692); // type 1
  GMM_m_w_.push_back(2.20402); // type 2
  GMM_m_w_.push_back(5.14099); // type 3

  // check if MOLINFO line is present
  if( moldat ) {
    log<<"  MOLINFO DATA found with label " <<moldat->getLabel()<<", using proper atom names\n";
    for(unsigned i=0; i<atoms.size(); ++i) {
      // get atom name
      std::string name = moldat->getAtomName(atoms[i]);
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

void EMMI::check_GMM_d(const VectorGeneric<6> &cov, double w)
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
void EMMI::get_GMM_d(std::string GMM_file)
{
  VectorGeneric<6> cov;

// open file
  auto ifile=Tools::make_unique<IFile>();
  if(ifile->FileExist(GMM_file)) {
    ifile->open(GMM_file);
    int idcomp;
    while(ifile->scanField("Id",idcomp)) {
      int beta;
      double w, m0, m1, m2;
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
  } else {
    error("Cannot find GMM_FILE "+GMM_file+"\n");
  }
  // now create a set from beta (unique set of values)
  std::set<int> bu(GMM_d_beta_.begin(), GMM_d_beta_.end());
  // now prepare the group std::vector
  GMM_d_grps_.resize(bu.size());
  // and fill it in
  for(unsigned i=0; i<GMM_d_beta_.size(); ++i) {
    if(GMM_d_beta_[i]>=GMM_d_grps_.size()) error("Check Beta values");
    GMM_d_grps_[GMM_d_beta_[i]].push_back(i);
  }
}

void EMMI::calculate_useful_stuff(double reso)
{
  // We use the following definition for resolution:
  // the Fourier transform of the density distribution in real space
  // f(s) falls to 1/e of its maximum value at wavenumber 1/resolution
  // i.e. from f(s) = A * exp(-B*s**2) -> Res = std::sqrt(B).
  // average value of B
  double Bave = 0.0;
  for(unsigned i=0; i<GMM_m_type_.size(); ++i) {
    Bave += GMM_m_s_[GMM_m_type_[i]];
  }
  Bave /= static_cast<double>(GMM_m_type_.size());
  // calculate blur factor
  double blur = 0.0;
  if(reso*reso>Bave) blur = reso*reso-Bave;
  else warning("PLUMED should not be used with maps at resolution better than 0.3 nm");
  // add blur to B
  for(unsigned i=0; i<GMM_m_s_.size(); ++i) GMM_m_s_[i] += blur;
  // calculate average resolution
  double ave_res = 0.0;
  for(unsigned i=0; i<GMM_m_type_.size(); ++i) {
    ave_res += std::sqrt(GMM_m_s_[GMM_m_type_[i]]);
  }
  ave_res = ave_res / static_cast<double>(GMM_m_type_.size());
  log.printf("  experimental map resolution : %3.2f\n", reso);
  log.printf("  predicted map resolution : %3.2f\n", ave_res);
  log.printf("  blur factor : %f\n", blur);
  // now calculate useful stuff
  VectorGeneric<6> cov, sum, inv_sum;
  // cycle on all atoms types (4 for the moment)
  for(unsigned i=0; i<GMM_m_s_.size(); ++i) {
    // the Gaussian in density (real) space is the FT of scattering factor
    // f(r) = A * (pi/B)**1.5 * exp(-pi**2/B*r**2)
    double s = std::sqrt ( 0.5 * GMM_m_s_[i] ) / pi;
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
      // calculate prefactor - model weights are already normalized
      double pre_fact =  cfact_ / std::sqrt(det) * GMM_d_w_[j] * GMM_m_w_[i];
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
    tab_exp_.push_back(std::exp(-static_cast<double>(i) * dexp_));
  }
}

// get prefactors
double EMMI::get_prefactor_inverse
(const VectorGeneric<6> &GMM_cov_0, const VectorGeneric<6> &GMM_cov_1,
 double GMM_w_0, double GMM_w_1,
 VectorGeneric<6> &sum, VectorGeneric<6> &inv_sum)
{
// we need the sum of the covariance matrices
  for(unsigned k=0; k<6; ++k) sum[k] = GMM_cov_0[k] + GMM_cov_1[k];

// and to calculate its determinant
  double det = sum[0]*(sum[3]*sum[5]-sum[4]*sum[4]);
  det -= sum[1]*(sum[1]*sum[5]-sum[4]*sum[2]);
  det += sum[2]*(sum[1]*sum[4]-sum[3]*sum[2]);

// the prefactor is
  double pre_fact =  cfact_ / std::sqrt(det) * GMM_w_0 * GMM_w_1;

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
double EMMI::get_overlap(const Vector &m_m, const Vector &d_m, double pre_fact,
                         const VectorGeneric<6> &inv_cov_md, Vector &ov_der)
{
  Vector md;
  // calculate std::vector difference m_m-d_m with/without pbc
  if(pbc_) md = pbcDistance(d_m, m_m);
  else     md = delta(d_m, m_m);
  // calculate product of transpose of md and inv_cov_md
  double p_x = md[0]*inv_cov_md[0]+md[1]*inv_cov_md[1]+md[2]*inv_cov_md[2];
  double p_y = md[0]*inv_cov_md[1]+md[1]*inv_cov_md[3]+md[2]*inv_cov_md[4];
  double p_z = md[0]*inv_cov_md[2]+md[1]*inv_cov_md[4]+md[2]*inv_cov_md[5];
  // calculate product of prod and md
  double ov = md[0]*p_x+md[1]*p_y+md[2]*p_z;
  // final calculation
  ov = pre_fact * std::exp(-0.5*ov);
  // derivatives
  ov_der = ov * Vector(p_x, p_y, p_z);
  return ov;
}

// get the exponent of the overlap
double EMMI::get_exp_overlap(const Vector &m_m, const Vector &d_m,
                             const VectorGeneric<6> &inv_cov_md)
{
  Vector md;
  // calculate std::vector difference m_m-d_m with/without pbc
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
  // dimension of GMM and atom std::vectors
  unsigned GMM_d_size = GMM_d_m_.size();
  unsigned GMM_m_size = GMM_m_type_.size();
  // local neighbor list
  std::vector < unsigned > nl_l;
  // clear old neighbor list
  nl_.clear();

  // cycle on GMM components - in parallel
  for(unsigned id=rank_; id<GMM_d_size; id+=size_) {
    // overlap lists and map
    std::vector<double> ov_l;
    std::map<double, unsigned> ov_m;
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
      // in case calculate overlap
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
    for(std::map<double, unsigned>::iterator it=ov_m.begin(); it!=ov_m.end(); ++it)
      nl_l.push_back(id*GMM_m_size+it->second);
    // end cycle on GMM components in parallel
  }
  // find total dimension of neighborlist
  std::vector <int> recvcounts(size_, 0);
  recvcounts[rank_] = nl_l.size();
  comm.Sum(&recvcounts[0], size_);
  int tot_size = accumulate(recvcounts.begin(), recvcounts.end(), 0);
  // resize neighbor stuff
  nl_.resize(tot_size);
  // calculate std::vector of displacement
  std::vector<int> disp(size_);
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

  // clean temporary std::vectors
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

double EMMI::scaleEnergy(double s)
{
  double ene = 0.0;
  for(unsigned i=0; i<ovdd_.size(); ++i) {
    ene += std::log( abs ( s * ovmd_ave_[i] - ovdd_[i] ) );
  }
  return ene;
}

double EMMI::doRegression()
{
// standard MC parameters
  unsigned MCsteps = 100000;
  double kbtmin = 1.0;
  double kbtmax = 10.0;
  unsigned ncold = 5000;
  unsigned nhot = 2000;
  double MCacc = 0.0;
  double kbt, ebest, scale_best;

// initial value of scale factor and energy
  double scale = random_.RandU01() * ( scale_max_ - scale_min_ ) + scale_min_;
  double ene = scaleEnergy(scale);
// set best energy
  ebest = ene;

// MC loop
  for(unsigned istep=0; istep<MCsteps; ++istep) {
    // get temperature
    if(istep%(ncold+nhot)<ncold) kbt = kbtmin;
    else kbt = kbtmax;
    // propose move in scale
    double ds = dscale_ * ( 2.0 * random_.RandU01() - 1.0 );
    double new_scale = scale + ds;
    // check boundaries
    if(new_scale > scale_max_) {new_scale = 2.0 * scale_max_ - new_scale;}
    if(new_scale < scale_min_) {new_scale = 2.0 * scale_min_ - new_scale;}
    // new energy
    double new_ene = scaleEnergy(new_scale);
    // accept or reject
    bool accept = doAccept(ene, new_ene, kbt);
    // in case of acceptance
    if(accept) {
      scale = new_scale;
      ene = new_ene;
      MCacc += 1.0;
    }
    // save best
    if(ene<ebest) {
      ebest = ene;
      scale_best = scale;
    }
  }
// calculate acceptance
  double accscale = MCacc / static_cast<double>(MCsteps);
// global communication
  if(!no_aver_ && nrep_>1) {
    if(replica_!=0) {
      scale_best = 0.0;
      ebest = 0.0;
      accscale = 0.0;
    }
    if(rank_==0) {
      multi_sim_comm.Sum(&scale_best, 1);
      multi_sim_comm.Sum(&ebest, 1);
      multi_sim_comm.Sum(&accscale, 1);
    }
  }
  // local communication
  if(rank_!=0) {
    scale_best = 0.0;
    ebest = 0.0;
    accscale = 0.0;
  }
  if(size_>1) {
    comm.Sum(&scale_best, 1);
    comm.Sum(&ebest, 1);
    comm.Sum(&accscale, 1);
  }
// set scale parameters
  getPntrToComponent("accscale")->set(accscale);
  getPntrToComponent("enescale")->set(ebest);
// return scale value
  return scale_best;
}

double EMMI::get_annealing(long long int step)
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

void EMMI::get_weights(double &weight, double &norm, double &neff)
{
  const double dnrep = static_cast<double>(nrep_);
  // calculate the weights either from BIAS
  if(do_reweight_) {
    std::vector<double> bias(nrep_,0);
    if(rank_==0) {
      bias[replica_] = getArgument(0);
      if(nrep_>1) multi_sim_comm.Sum(&bias[0], nrep_);
    }
    comm.Sum(&bias[0], nrep_);

    // accumulate weights
    if(!first_time_w_) {
      for(unsigned i=0; i<nrep_; ++i) {
        const double delta=bias[i]-average_weights_[i];
        // FIXME: multiplying by fractional decay here causes problems with numerical derivatives,
        // probably because we're making several calls to calculate(), causing accumulation
        // of epsilons. Maybe we can work on a temporary copy of the action instead?
        average_weights_[i]+=decay_w_*delta;
      }
    } else {
      first_time_w_ = false;
      for(unsigned i=0; i<nrep_; ++i) average_weights_[i] = bias[i];
    }

    // set average back into bias and set norm to one
    const double maxbias = *(std::max_element(average_weights_.begin(), average_weights_.end()));
    for(unsigned i=0; i<nrep_; ++i) bias[i] = std::exp((average_weights_[i]-maxbias)/kbt_);
    // set local weight, norm and weight variance
    weight = bias[replica_];
    double w2=0.;
    for(unsigned i=0; i<nrep_; ++i) {
      w2 += bias[i]*bias[i];
      norm += bias[i];
    }
    neff = norm*norm/w2;
    getPntrToComponent("weight")->set(weight/norm);
  } else {
    // or arithmetic ones
    neff = dnrep;
    weight = 1.0;
    norm = dnrep;
  }
  getPntrToComponent("neff")->set(neff);
}

void EMMI::calculate()
{

// calculate CV
  calculate_overlap();

  // rescale factor for ensemble average
  double weight = 0.;
  double neff = 0.;
  double norm = 0.;
  get_weights(weight, norm, neff);

  // in case of ensemble averaging, calculate average overlap
  if(!no_aver_ && nrep_>1) {
    // if master node, calculate average across replicas
    if(rank_==0) {
      for(unsigned i=0; i<ovmd_.size(); ++i) ovmd_ave_[i] = weight / norm * ovmd_[i];
      multi_sim_comm.Sum(&ovmd_ave_[0], ovmd_ave_.size());
    } else {
      for(unsigned i=0; i<ovmd_ave_.size(); ++i) ovmd_ave_[i] = 0.0;
    }
    // local communication
    if(size_>1) comm.Sum(&ovmd_ave_[0], ovmd_ave_.size());
  } else {
    for(unsigned i=0; i<ovmd_.size(); ++i) ovmd_ave_[i] = ovmd_[i];
  }

  // get time step
  long long int step = getStep();

  // do regression
  if(nregres_>0) {
    if(step%nregres_==0 && !getExchangeStep()) scale_ = doRegression();
    // set scale component
    getPntrToComponent("scale")->set(scale_);
  }

  // write model overlap to file
  if(ovstride_>0 && step%ovstride_==0) write_model_overlap(step);

  // clear energy and virial
  ene_ = 0.0;
  virial_.zero();

  // Gaussian noise
  if(noise_==0) calculate_Gauss();

  // Outliers noise
  if(noise_==1) calculate_Outliers();

  // Marginal noise
  if(noise_==2) calculate_Marginal();

  // get annealing rescale factor
  if(nanneal_>0) {
    anneal_ = get_annealing(step);
    getPntrToComponent("anneal")->set(anneal_);
  }

  // annealing rescale
  ene_ /= anneal_;

  std::vector<double> GMMid_der_av_(GMMid_der_.size(), 0.0);
  // in case of ensemble averaging
  if(!no_aver_ && nrep_>1) {
    // if master node, sum der_GMMid derivatives and ene
    if(rank_==0) {
      for(unsigned i=0; i<GMMid_der_.size(); ++i) {
        GMMid_der_av_[i] = (weight / norm) * GMMid_der_[i];
      }
      multi_sim_comm.Sum(&GMMid_der_av_[0], GMMid_der_av_.size());
      multi_sim_comm.Sum(&ene_, 1);
    } else {
      // set der_GMMid derivatives and energy to zero
      for(unsigned i=0; i<GMMid_der_av_.size(); ++i) GMMid_der_av_[i]=0.0;
      ene_ = 0.0;
    }
    // local communication
    if(size_>1) {
      comm.Sum(&GMMid_der_av_[0], GMMid_der_av_.size());
      comm.Sum(&ene_, 1);
    }
  } else {
    for (unsigned i = 0; i < GMMid_der_.size(); ++i) {
      GMMid_der_av_[i] = GMMid_der_[i];
    }
  }

  // clean temporary std::vector
  for(unsigned i=0; i<atom_der_.size(); ++i) atom_der_[i] = Vector(0,0,0);

  // get derivatives of bias with respect to atoms
  for(unsigned i=rank_; i<nl_.size(); i=i+size_) {
    // get indexes of data and model component
    unsigned id = nl_[i] / GMM_m_type_.size();
    unsigned im = nl_[i] % GMM_m_type_.size();
    // chain rule + replica normalization
    Vector tot_der = GMMid_der_av_[id] * ovmd_der_[i] * scale_ / anneal_;
    Vector pos;
    if(pbc_) pos = pbcDistance(GMM_d_m_[id], getPosition(im)) + GMM_d_m_[id];
    else     pos = getPosition(im);
    // increment derivatives and virial
    atom_der_[im] += tot_der;
    virial_ += Tensor(pos, -tot_der);
  }

  // communicate local derivatives and virial
  if(size_>1) {
    comm.Sum(&atom_der_[0][0], 3*atom_der_.size());
    comm.Sum(virial_);
  }

  // set derivatives, virial, and score
  for(unsigned i=0; i<atom_der_.size(); ++i) setAtomsDerivatives(getPntrToComponent("scoreb"), i, atom_der_[i]);
  setBoxDerivatives(getPntrToComponent("scoreb"), virial_);
  getPntrToComponent("scoreb")->set(ene_);

  if (do_reweight_) {
    double w_tmp = 0.;
    for (unsigned i = 0; i < ovmd_.size(); ++i) {
      w_tmp += (ovmd_[i] - ovmd_ave_[i]) * GMMid_der_[i];
    }
    w_tmp *= scale_ * (weight / norm) / kbt_ * decay_w_;

    setArgDerivatives(getPntrToComponent("scoreb"), w_tmp);
    getPntrToComponent("biasDer")->set(w_tmp);
  }

  // This part is needed only for Gaussian and Outliers noise models
  if(noise_!=2) {

    // do Montecarlo
    if(dsigma_[0]>0 && step%MCstride_==0 && !getExchangeStep()) doMonteCarlo();

    // print status
    if(step%statusstride_==0) print_status(step);

    // calculate acceptance ratio
    double acc = MCaccept_ / MCtrials_;

    // set value
    getPntrToComponent("acc")->set(acc);

  }

}

void EMMI::calculate_Gauss()
{
  // cycle on all the GMM groups
  for(unsigned i=0; i<GMM_d_grps_.size(); ++i) {
    double eneg = 0.0;
    // cycle on all the members of the group
    for(unsigned j=0; j<GMM_d_grps_[i].size(); ++j) {
      // id of the GMM component
      int GMMid = GMM_d_grps_[i][j];
      // calculate deviation
      double dev = ( scale_*ovmd_ave_[GMMid]-ovdd_[GMMid] ) / sigma_[i];
      // add to group energy
      eneg += 0.5 * dev * dev;
      // store derivative for later
      GMMid_der_[GMMid] = kbt_ * dev / sigma_[i];
    }
    // add to total energy along with normalizations and prior
    ene_ += kbt_ * ( eneg + (static_cast<double>(GMM_d_grps_[i].size())+prior_) * std::log(sigma_[i]) );
  }
}

void EMMI::calculate_Outliers()
{
  // cycle on all the GMM groups
  for(unsigned i=0; i<GMM_d_grps_.size(); ++i) {
    // cycle on all the members of the group
    double eneg = 0.0;
    for(unsigned j=0; j<GMM_d_grps_[i].size(); ++j) {
      // id of the GMM component
      int GMMid = GMM_d_grps_[i][j];
      // calculate deviation
      double dev = ( scale_*ovmd_ave_[GMMid]-ovdd_[GMMid] ) / sigma_[i];
      // add to group energy
      eneg += std::log1p( 0.5 * dev * dev );
      // store derivative for later
      GMMid_der_[GMMid] = kbt_ / ( 1.0 + 0.5 * dev * dev ) * dev / sigma_[i];
    }
    // add to total energy along with normalizations and prior
    ene_ += kbt_ * ( eneg + (static_cast<double>(GMM_d_grps_[i].size())+prior_) * std::log(sigma_[i]) );
  }
}

void EMMI::calculate_Marginal()
{
  // cycle on all the GMM groups
  for(unsigned i=0; i<GMM_d_grps_.size(); ++i) {
    // cycle on all the members of the group
    for(unsigned j=0; j<GMM_d_grps_[i].size(); ++j) {
      // id of the GMM component
      int GMMid = GMM_d_grps_[i][j];
      // calculate deviation
      double dev = ( scale_*ovmd_ave_[GMMid]-ovdd_[GMMid] );
      // calculate errf
      double errf = erf ( dev * inv_sqrt2_ / sigma_min_[i] );
      // add to group energy
      ene_ += -kbt_ * std::log ( 0.5 / dev * errf ) ;
      // store derivative for later
      GMMid_der_[GMMid] = - kbt_/errf*sqrt2_pi_*std::exp(-0.5*dev*dev/sigma_min_[i]/sigma_min_[i])/sigma_min_[i]+kbt_/dev;
    }
  }
}

}
}
