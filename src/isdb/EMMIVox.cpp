/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017,2018 The plumed team
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
#include "core/GenericMolInfo.h"
#include "core/ActionSet.h"
#include "tools/File.h"
#include "tools/OpenMP.h"
#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <numeric>
#include <ctime>
#include "tools/Random.h"
#include "simd_math_prims.h"
#ifdef __PLUMED_HAS_ARRAYFIRE
#include <arrayfire.h>
#include <af/util.h>
#endif

using namespace std;

namespace PLMD {
namespace isdb {

//+PLUMEDOC ISDB_COLVAR EMMIVOX
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
To use \ref EMMIVOX, the user should always add a \ref MOLINFO line and specify a pdb file of the system.

\note
To enhance sampling in single-structure refinement, one can use a Replica Exchange Method, such as Parallel Tempering.
In this case, the user should add the NO_AVER flag to the input line.

\note
\ref EMMIVOX can be used in combination with periodic and non-periodic systems. In the latter case, one should
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

# create EMMIVOX score
gmm: EMMIVOX NOPBC SIGMA_MIN=0.01 TEMP=300.0 NL_STRIDE=100 NL_CUTOFF=0.01 GMM_FILE=GMM_fit.dat ATOMS=protein-h

# translate into bias - apply every 2 steps
emr: BIASVALUE ARG=gmm.scoreb STRIDE=2

PRINT ARG=emr.* FILE=COLVAR STRIDE=500 FMT=%20.10f
\endplumedfile


*/
//+ENDPLUMEDOC

class EMMIVOX : public Colvar {

private:

// temperature in kbt
  double kbt_;
// model GMM - atom types
  vector<unsigned> GMM_m_type_;
// model GMM - list of atom sigmas - one per atom type
  vector<double> GMM_m_s0_;
  vector<Vector5d> GMM_m_s_;
// model GMM - list of atom weights - one per atom type
  vector<Vector5d> GMM_m_w_;
// model GMM - map between residue number and list of atoms
  map< unsigned, vector<unsigned> > GMM_m_resmap_;
// model GMM - list of residue ids
  vector<unsigned> GMM_m_res_;
// model GMM - list of neighboring voxels per atom
  vector< vector<unsigned> > GMM_m_nb_;
// model GMM - map between res id and bfactor
  map<unsigned,double> GMM_m_b_;
// model overlap
  vector<double> ovmd_;

// data GMM - means, sigma2 + beta option
  vector<Vector> GMM_d_m_;
  double         GMM_d_s_;
  vector<int>    GMM_d_beta_;
// data GMM - groups bookeeping
  vector < vector<int> > GMM_d_grps_;
// model GMM - list of neighboring atoms per voxel
  vector< vector<unsigned> > GMM_d_nb_;
// data GMM - overlap
  vector<double> ovdd_;

// derivatives
  vector<Vector> ovmd_der_;
  vector<Vector> atom_der_;
  vector<double> GMMid_der_;
// constants
  double inv_sqrt2_, sqrt2_pi_, inv_pi2_;
  vector<Vector5d> pref_;
  vector<Vector5d> invs2_;
  vector<Vector5d> cfact_;
// metainference
  unsigned nrep_;
  unsigned replica_;
  vector<double> sigma_min_;
  vector<double> ismin_;
// neighbor list
  double   nl_cutoff_;
  unsigned nl_stride_;
  double   ns_cutoff_;
  bool first_time_;
  vector<unsigned> nl_;
  vector<unsigned> ns_;
  vector<Vector> refpos_;
  vector<Vector> pos_;
// averaging
  bool no_aver_;
// correlation;
  bool do_corr_;
// Monte Carlo stuff
  Random   random_;
  // Scale Monte Carlo
  double scale_;
  double scale_min_;
  double scale_max_;
  double dscale_;
  double offset_;
  double doffset_;
  int    MCSstride_;
  double MCSaccept_;
  double MCStrials_;
// Bfact Monte Carlo
  double   dbfact_;
  double   bfactmin_;
  double   bfactmax_;
  int      MCBstride_;
  double   MCBaccept_;
  double   MCBtrials_;
  // status stuff
  unsigned int statusstride_;
  string       statusfilename_;
  OFile        statusfile_;
  bool         first_status_;
  // total energy and virial
  double ene_;
  Tensor virial_;
  // model overlap file
  unsigned int ovstride_;
  string       ovfilename_;
  // gpu stuff
  bool gpu_;
  int  deviceid_;
#ifdef __PLUMED_HAS_ARRAYFIRE
  // gpu stuff
  vector<float> pos_gpu_;
  vector<float> der_gpu_;
  vector<float> ovmd_gpu_;
  af::array pref_gpu;
  af::array invs2_gpu;
  af::array GMM_d_m_gpu;
  af::array ismin_gpu;
  af::array ovdd_gpu;
  // gpu nl stuff
  af::array pref_gpu_nl;
  af::array invs2_gpu_nl;
  af::array GMM_d_m_gpu_nl;
  af::array id_gpu;
  af::array im_gpu;
  af::array ov_k_sort, ov_k_keys;
  af::array d_k_sort, d_k_keys;
  // some constants
  float kbt, inv_sqrt2, sqrt2_pi;
#endif
//
// write file with model overlap
  void write_model_overlap(long int step);
// get median of vector
  double get_median(vector<double> &v);
// read and write status
  void read_status();
  void print_status(long int step);
// accept or reject
  bool doAccept(double oldE, double newE, double kbt);
// do MonteCarlo for Bfactor
  void doMonteCarloBfact();
// do MonteCarlo for scale
  void doMonteCarloScale();
// read error file
  vector<double> read_exp_errors(string errfile);
// calculate model GMM parameters
  vector<double> get_GMM_m(vector<AtomNumber> &atoms);
// read data file
  void get_exp_data(string datafile);
// auxiliary methods
  void prepare_gpu();
  void calculate_useful_stuff(double reso);
  void get_auxiliary_vectors();
  void push_auxiliary_gpu();
// calculate overlap between two Gaussians
  double get_overlap_der(const Vector &d_m, const Vector &m_m,
                         const Vector5d &pref, const Vector5d &invs2,
                         Vector &ov_der);
  double get_overlap(const Vector &d_m, const Vector &m_m, double d_s,
                     const Vector5d &cfact, const Vector5d &m_s, double bfact);
// update the neighbor list
  void update_neighbor_list();
// update the gpu
  void update_gpu();
// update the neighbor sphere
  void update_neighbor_sphere();
  bool do_neighbor_sphere();
// calculate on cpu/gpu
  void calculate_cpu();
  void calculate_gpu();
// calculate correlation
  double calculate_corr();
// Marginal noise
  double calculate_Marginal(double scale, double offset, vector<double> &GMMid_der);
  double calculate_Marginal(double scale, double offset);

public:
  static void registerKeywords( Keywords& keys );
  explicit EMMIVOX(const ActionOptions&);
// active methods:
  void prepare();
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(EMMIVOX,"EMMIVOX")

void EMMIVOX::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOMS","atoms for which we calculate the density map, typically all heavy atoms");
  keys.add("compulsory","DATA_FILE","file with the experimental data");
  keys.add("compulsory","NL_CUTOFF","The cutoff in distance for the neighbor list");
  keys.add("compulsory","NL_STRIDE","The frequency with which we are updating the neighbor list");
  keys.add("compulsory","NS_CUTOFF","The cutoff in distance for the outer neighbor sphere");
  keys.add("compulsory","SIGMA_MIN","minimum uncertainty");
  keys.add("compulsory","RESOLUTION", "Cryo-EM map resolution");
  keys.add("compulsory","VOXEL","Side of voxel grid");
  keys.add("compulsory","NORM_DENSITY","integral of the experimental density");
  keys.add("compulsory","WRITE_STRIDE","write the status to a file every N steps, this can be used for restart");
  keys.add("optional","DEVICEID","Identifier of the GPU to be used");
  keys.add("optional","DBFACT","MC step for bfactor");
  keys.add("optional","BFACT_MAX","Maximum value of bfactor");
  keys.add("optional","MCBFACT_STRIDE", "Bfactor Monte Carlo stride");
  keys.add("optional","ERR_FILE","file with experimental errors");
  keys.add("optional","STATUS_FILE","write a file with all the data useful for restart");
  keys.add("optional","SCALE_MIN","minimum scale");
  keys.add("optional","SCALE_MAX","maximum scale");
  keys.add("optional","DSCALE","maximum scale MC move");
  keys.add("optional","SCALE","scale factor");
  keys.add("optional","OFFSET","offset");
  keys.add("optional","DOFFSET","maximum offset MC move");
  keys.add("optional","MCSCALE_STRIDE", "scale factor Monte Carlo stride");
  keys.add("optional","TEMP","temperature");
  keys.add("optional","WRITE_OV_STRIDE","write model overlaps every N steps");
  keys.add("optional","WRITE_OV","write a file with model overlaps");
  keys.addFlag("NO_AVER",false,"don't do ensemble averaging in multi-replica mode");
  keys.addFlag("CORRELATION",false,"calculate correlation coefficient");
  keys.addFlag("GPU",false,"calculate EMMIVOX using ARRAYFIRE on an accelerator device");
  componentsAreNotOptional(keys);
  keys.addOutputComponent("scoreb","default","Bayesian score");
  keys.addOutputComponent("scale", "default","scale factor");
  keys.addOutputComponent("offset","default","offset");
  keys.addOutputComponent("accS",  "default","MC acceptance for scale");
  keys.addOutputComponent("accB",  "default", "Bfactor MC acceptance");
  keys.addOutputComponent("corr", "CORRELATION", "correlation coefficient");
}

EMMIVOX::EMMIVOX(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  inv_sqrt2_(0.707106781186548),
  sqrt2_pi_(0.797884560802865),
  inv_pi2_(0.050660591821169),
  first_time_(true), no_aver_(false), do_corr_(false),
  scale_(1.), dscale_(0.), offset_(0.), doffset_(0.),
  MCSstride_(1), MCSaccept_(0.), MCStrials_(0.),
  dbfact_(0.0), bfactmax_(4.0),
  MCBstride_(1), MCBaccept_(0.), MCBtrials_(0.),
  statusstride_(0), first_status_(true),
  ovstride_(0), gpu_(false), deviceid_(0)
{

  // list of atoms
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS", atoms);

  // file with experimental data
  string datafile;
  parse("DATA_FILE", datafile);

  // minimum value for error
  double sigma_min;
  parse("SIGMA_MIN", sigma_min);
  if(sigma_min<0) error("SIGMA_MIN should be greater or equal to zero");

  // status file parameters
  parse("WRITE_STRIDE", statusstride_);
  if(statusstride_<=0) error("you must specify a positive WRITE_STRIDE");
  parse("STATUS_FILE",  statusfilename_);
  if(statusfilename_=="") statusfilename_ = "MISTATUS"+getLabel();
  else                    statusfilename_ = statusfilename_+getLabel();

  // error file
  string errfile;
  parse("ERR_FILE", errfile);

  // integral of the experimetal density
  double norm_d;
  parse("NORM_DENSITY", norm_d);

  // temperature
  double temp=0.0;
  parse("TEMP",temp);
  // convert temp to kbt
  if(temp>0.0) kbt_=plumed.getAtoms().getKBoltzmann()*temp;
  else kbt_=plumed.getAtoms().getKbT();

  // scale MC
  parse("SCALE", scale_);
  parse("DSCALE",dscale_);
  parse("OFFSET",offset_);
  parse("DOFFSET",doffset_);
  parse("MCSCALE_STRIDE",MCSstride_);
  // other parameters
  if(dscale_>0.) {
    parse("SCALE_MIN",scale_min_);
    parse("SCALE_MAX",scale_max_);
    // checks
    if(MCSstride_<=0) error("you must specify a positive MCSCALE_STRIDE");
    if(scale_min_<=0.0) error("SCALE_MIN must be strictly positive");
    if(scale_max_<=scale_min_) error("SCALE_MAX must be greater than SCALE_MIN");
  }

  // Monte Carlo in B-factors
  parse("DBFACT",dbfact_);
  parse("BFACT_MAX",bfactmax_);
  parse("MCBFACT_STRIDE",MCBstride_);
  if(dbfact_<0) error("DBFACT should be greater or equal to zero");
  if(dbfact_>0 && MCBstride_<=0) error("you must specify a positive MCBFACT_STRIDE");
  if(dbfact_>0 && bfactmax_<=0)  error("you must specify a positive BFACT_MAX");

  // read map resolution
  double reso;
  parse("RESOLUTION", reso);
  if(reso<=0.) error("RESOLUTION should be strictly positive");

  // voxel size
  parse("VOXEL", GMM_d_s_);
  // transform into sigma^2 for overlap calculation
  GMM_d_s_ = pow(GMM_d_s_/4.0, 2.0);

  // neighbor list stuff
  parse("NL_CUTOFF",nl_cutoff_);
  if(nl_cutoff_<=0.0) error("NL_CUTOFF should be explicitly specified and positive");
  parse("NL_STRIDE",nl_stride_);
  if(nl_stride_<=0) error("NL_STRIDE should be explicitly specified and positive");
  parse("NS_CUTOFF",ns_cutoff_);
  if(ns_cutoff_<=nl_cutoff_) error("NS_CUTOFF should be greater than NL_CUTOFF");

  // averaging or not
  parseFlag("NO_AVER",no_aver_);

  // calculate correlation coefficient
  parseFlag("CORRELATION",do_corr_);

  // write overlap file
  parse("WRITE_OV_STRIDE", ovstride_);
  parse("WRITE_OV", ovfilename_);
  if(ovstride_>0 && ovfilename_=="") error("With WRITE_OV_STRIDE you must specify WRITE_OV");

  parseFlag("GPU",gpu_);
#ifndef  __PLUMED_HAS_ARRAYFIRE
  if(gpu_) error("To use the GPU mode PLUMED must be compiled with ARRAYFIRE");
#endif

  parse("DEVICEID",deviceid_);
#ifdef  __PLUMED_HAS_ARRAYFIRE
  af::setDevice(deviceid_);
  af::info();
#endif

  checkRead();

  // set parallel stuff
  unsigned mpisize=comm.Get_size();
  if(mpisize>1) error("EMMIVOX supports only OpenMP parallelization");

  // get number of replicas
  if(no_aver_) {
    nrep_ = 1;
  } else {
    nrep_ = multi_sim_comm.Get_size();
  }
  replica_ = multi_sim_comm.Get_rank();

  if(nrep_>1 && dbfact_>0) error("Bfactor sampling not supported with ensemble averaging");

  log.printf("  atoms involved : ");
  for(unsigned i=0; i<atoms.size(); ++i) log.printf("%d ",atoms[i].serial());
  log.printf("\n");
  log.printf("  experimental data file : %s\n", datafile.c_str());
  if(no_aver_) log.printf("  without ensemble averaging\n");
  if(gpu_) log.printf("  running on GPU on device: %d\n", deviceid_);
  log.printf("  neighbor list cutoff : %lf\n", nl_cutoff_);
  log.printf("  neighbor list stride : %u\n",  nl_stride_);
  log.printf("  neighbor sphere cutoff : %lf\n", ns_cutoff_);
  log.printf("  minimum uncertainty : %f\n",sigma_min);
  log.printf("  scale factor : %lf\n",scale_);
  log.printf("  offset : %lf\n",offset_);
  log.printf("  reading/writing to status file : %s\n",statusfilename_.c_str());
  log.printf("  with stride : %u\n",statusstride_);
  if(dscale_>0.0) {
    log.printf("  minimum scale : %lf\n", scale_min_);
    log.printf("  maximum scale : %lf\n", scale_max_);
    log.printf("  maximum scale MC move : %lf\n", dscale_);
    log.printf("  maximum offset MC move : %lf\n", doffset_);
    log.printf("  scale factor MC stride : %u\n", MCSstride_);
  }
  if(dbfact_>0) {
    log.printf("  max MC move in bfactor : %f\n",dbfact_);
    log.printf("  Bfactor MC stride : %u\n", MCBstride_);
  }
  if(errfile.size()>0) log.printf("  reading experimental errors from file : %s\n", errfile.c_str());
  log.printf("  temperature of the system in energy unit : %f\n",kbt_);
  log.printf("  number of replicas for averaging: %u\n",nrep_);
  log.printf("  id of the replica : %u\n",replica_);
  if(ovstride_>0) {
    log.printf("  stride for writing model overlaps : %u\n",ovstride_);
    log.printf("  file for writing model overlaps : %s\n", ovfilename_.c_str());
  }

  // calculate model GMM constant parameters
  vector<double> GMM_m_w = get_GMM_m(atoms);

  // prepare lists
  pos_.resize(GMM_m_w.size());
  refpos_.resize(GMM_m_w.size());

  // read data file
  get_exp_data(datafile);
  log.printf("  number of voxels : %u\n", static_cast<unsigned>(GMM_d_m_.size()));

  // normalize atom weight map
  double norm_m = accumulate(GMM_m_w.begin(),  GMM_m_w.end(),  0.0);
  // renormalization and constant factor
  for(unsigned i=0; i<GMM_m_w_.size(); ++i) {
    Vector5d cf;
    for(unsigned j=0; j<5; ++j) {
      GMM_m_w_[i][j] *= norm_d / norm_m;
      cf[j] = GMM_m_w_[i][j]/pow( 2.0*pi, 1.5 );
    }
    cfact_.push_back(cf);
  }

  // read experimental errors
  vector<double> exp_err;
  if(errfile.size()>0) exp_err = read_exp_errors(errfile);

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
    double ovdd_m = get_median(ovdd);
    double err_m  = get_median(err);
    // print out statistics
    log.printf("     # of members : %u\n", GMM_d_grps_[Gid].size());
    log.printf("     median overlap : %lf\n", ovdd_m);
    log.printf("     median error : %lf\n", err_m);
    // add minimum value of sigma for this group of GMMs
    sigma_min_.push_back(sqrt(err_m*err_m+sigma_min*ovdd_m*sigma_min*ovdd_m));
  }
  // populate ismin: cycle on all ovdd
  for(unsigned id=0; id<GMM_d_beta_.size(); ++id) {
    // id of the group
    unsigned ig = GMM_d_beta_[id];
    // and to ismin_
    ismin_.push_back(1.0 / sigma_min_[ig]);
  }

  // prepare gpu stuff
  if(gpu_) prepare_gpu();

  // calculate constant and auxiliary stuff
  calculate_useful_stuff(reso);

  // read status file if restarting
  if(getRestart()) read_status();

  // prepare auxiliary vectors
  get_auxiliary_vectors();

  // prepare data and derivative vectors
  ovmd_.resize(ovdd_.size());
  atom_der_.resize(GMM_m_type_.size());
  GMMid_der_.resize(ovdd_.size());

  // add components
  addComponentWithDerivatives("scoreb"); componentIsNotPeriodic("scoreb");
  addComponent("scale");                 componentIsNotPeriodic("scale");
  addComponent("offset");                componentIsNotPeriodic("offset");
  if(dbfact_>0)   {addComponent("accB"); componentIsNotPeriodic("accB");}
  if(dscale_>0.0) {addComponent("accS"); componentIsNotPeriodic("accS");}
  if(do_corr_)    {addComponent("corr"); componentIsNotPeriodic("corr");}

  // initialize random seed
  unsigned iseed = time(NULL)+replica_;
  random_.setSeed(-iseed);

  // request the atoms
  requestAtoms(atoms);

  // print bibliography
  log<<"  Bibliography "<<plumed.cite("Bonomi, Camilloni, Bioinformatics, 33, 3999 (2017)");
  log<<plumed.cite("Bonomi, Hanot, Greenberg, Sali, Nilges, Vendruscolo, Pellarin, Structure, 27, 175 (2019)");
  log<<plumed.cite("Bonomi, Pellarin, Vendruscolo, Biophys. J. 114, 1604 (2018)");
  if(!no_aver_ && nrep_>1)log<<plumed.cite("Bonomi, Camilloni, Cavalli, Vendruscolo, Sci. Adv. 2, e150117 (2016)");
  log<<"\n";
}

void EMMIVOX::prepare_gpu()
{
#ifdef __PLUMED_HAS_ARRAYFIRE
  // 1) convert constants to float
  kbt       = static_cast<float>(kbt_);
  inv_sqrt2 = static_cast<float>(inv_sqrt2_);
  sqrt2_pi  = static_cast<float>(sqrt2_pi_);
  // 2) create float version of ismin_
  vector<float> ismin_f;
  for(unsigned i=0; i<ismin_.size(); ++i) {
    ismin_f.push_back(static_cast<float>(ismin_[i]));
  }
  // create arrayfire [GMM_d_size, 1]
  ismin_gpu = af::array(ismin_f.size(), &ismin_f.front());
  // 3) create float version of ovdd_
  vector<float> ovdd_f;
  for(unsigned i=0; i<ovdd_.size(); ++i) {
    ovdd_f.push_back(static_cast<float>(ovdd_[i]));
  }
  // create arrayfire [GMM_d_size, 1]
  ovdd_gpu = af::array(ovdd_f.size(), &ovdd_f.front());
  // 4) store GMM_d_m_ as flattened vector of floats
  const unsigned GMM_d_size = GMM_d_m_.size();
  vector<float> GMM_d_m_gpu_(3*GMM_d_size);
  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for(unsigned i=0; i<GMM_d_size; ++i) {
    GMM_d_m_gpu_[i]              = static_cast<float>(GMM_d_m_[i][0]);
    GMM_d_m_gpu_[GMM_d_size+i]   = static_cast<float>(GMM_d_m_[i][1]);
    GMM_d_m_gpu_[2*GMM_d_size+i] = static_cast<float>(GMM_d_m_[i][2]);
  }
  // create array [GMM_d_size, 3]
  GMM_d_m_gpu = af::array(GMM_d_size, 3, &GMM_d_m_gpu_.front());
  // 5) dimensionate other vectors
  pos_gpu_.resize(3*GMM_m_type_.size());
  der_gpu_.resize(3*GMM_m_type_.size());
  ovmd_gpu_.resize(GMM_d_size);
#endif
}

void EMMIVOX::write_model_overlap(long int step)
{
  OFile ovfile;
  ovfile.link(*this);
  std::string num; Tools::convert(step,num);
  string name = ovfilename_+"-"+num;
  ovfile.open(name);
  ovfile.setHeavyFlush();
  ovfile.fmtField("%10.7e ");
// write overlaps
  for(unsigned i=0; i<ovmd_.size(); ++i) {
    ovfile.printField("Model", ovmd_[i]);
    ovfile.printField("ModelScaled", scale_ * ovmd_[i] + offset_);
    ovfile.printField("Data", ovdd_[i]);
    ovfile.printField();
  }
  ovfile.close();
}

double EMMIVOX::get_median(vector<double> &v)
{
// dimension of vector
  unsigned size = v.size();
// in case of only one entry
  if (size==1) {
    return v[0];
  } else {
    // reorder vector
    sort(v.begin(), v.end());
    // odd or even?
    if (size%2==0) {
      return (v[size/2-1]+v[size/2])/2.0;
    } else {
      return v[size/2];
    }
  }
}

void EMMIVOX::read_status()
{
  double MDtime;
// open file
  IFile *ifile = new IFile();
  ifile->link(*this);
  if(ifile->FileExist(statusfilename_)) {
    ifile->open(statusfilename_);
    while(ifile->scanField("MD_time", MDtime)) {
      // read scale and offset
      ifile->scanField("scale", scale_);
      ifile->scanField("offset", offset_);
      // read bfactors
      if(dbfact_>0) {
        for(unsigned i=0; i<GMM_m_res_.size(); ++i) {
          // convert i to string
          std::string num; Tools::convert(i,num);
          // read entries
          ifile->scanField("bfact"+num, GMM_m_b_[GMM_m_res_[i]]);
        }
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

void EMMIVOX::print_status(long int step)
{
// if first time open the file
  if(first_status_) {
    first_status_ = false;
    statusfile_.link(*this);
    statusfile_.open(statusfilename_);
    statusfile_.setHeavyFlush();
    statusfile_.fmtField("%10.7e ");
  }
// write fields
  double MDtime = static_cast<double>(step)*getTimeStep();
  statusfile_.printField("MD_time", MDtime);
  // write scale and offset
  statusfile_.printField("scale", scale_);
  statusfile_.printField("offset", offset_);
  // write bfactors only if doing fitting
  if(dbfact_>0) {
    for(unsigned i=0; i<GMM_m_res_.size(); ++i) {
      // convert i to string
      std::string num; Tools::convert(i,num);
      // print entry
      statusfile_.printField("bfact"+num, GMM_m_b_[GMM_m_res_[i]]);
    }
  }
  statusfile_.printField();
}

bool EMMIVOX::doAccept(double oldE, double newE, double kbt) {
  bool accept = false;
  // calculate delta energy
  double delta = ( newE - oldE ) / kbt;
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

vector<double> EMMIVOX::read_exp_errors(string errfile)
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
      for(int i=0; i<nexp; ++i) {
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

vector<double> EMMIVOX::get_GMM_m(vector<AtomNumber> &atoms)
{
  // list of weights - one per atom
  vector<double> GMM_m_w;

  auto* moldat=plumed.getActionSet().selectLatest<GenericMolInfo*>(this);
  // map of atom types to A and B coefficients of scattering factor
  // f(s) = A * exp(-B*s**2)
  // B is in Angstrom squared
  // map between an atom type and an index
  map<string, unsigned> type_map;
  type_map["C"]=0;
  type_map["O"]=1;
  type_map["N"]=2;
  type_map["S"]=3;
  // fill in sigma0 vector
  GMM_m_s0_.push_back(0.01*15.146);  // type 0
  GMM_m_s0_.push_back(0.01*8.59722); // type 1
  GMM_m_s0_.push_back(0.01*11.1116); // type 2
  GMM_m_s0_.push_back(0.01*15.8952); // type 3
  // fill in sigma vector
  GMM_m_s_.push_back(0.01*Vector5d(0.114,1.0825,5.4281,17.8811,51.1341));   // type 0
  GMM_m_s_.push_back(0.01*Vector5d(0.0652,0.6184,2.9449,9.6298,28.2194));   // type 1
  GMM_m_s_.push_back(0.01*Vector5d(0.0541,0.5165,2.8207,10.6297,34.3764));  // type 2
  GMM_m_s_.push_back(0.01*Vector5d(0.0838,0.7788,4.3462,15.5846,44.63655)); // type 3
  // fill in weight vector
  GMM_m_w_.push_back(Vector5d(0.0489,0.2091,0.7537,1.1420,0.3555)); // type 0
  GMM_m_w_.push_back(Vector5d(0.0365,0.1729,0.5805,0.8814,0.3121)); // type 1
  GMM_m_w_.push_back(Vector5d(0.0267,0.1328,0.5301,1.1020,0.4215)); // type 2
  GMM_m_w_.push_back(Vector5d(0.0915,0.4312,1.0847,2.4671,1.0852)); // type 3

  // check if MOLINFO line is present
  if( moldat ) {
    log<<"  MOLINFO DATA found with label " <<moldat->getLabel()<<", using proper atom names\n";
    for(unsigned i=0; i<atoms.size(); ++i) {
      // get atom name
      string name = moldat->getAtomName(atoms[i]);
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
        Vector5d w = GMM_m_w_[type_map[type_s]];
        GMM_m_w.push_back(w[0]+w[1]+w[2]+w[3]+w[4]);
        // get residue id
        unsigned ires = moldat->getResidueNumber(atoms[i]);
        // add to map and list
        GMM_m_resmap_[ires].push_back(i);
        GMM_m_res_.push_back(ires);
        // initialize Bfactor map
        GMM_m_b_[ires] = 0.0;
      } else {
        error("Wrong atom type "+type_s+" from atom name "+name+"\n");
      }
    }
  } else {
    error("MOLINFO DATA not found\n");
  }
  return GMM_m_w;
}

// read experimental data file in PLUMED format:
void EMMIVOX::get_exp_data(string datafile)
{
  Vector pos;
  double dens, beta;
  int idcomp;

// open file
  IFile *ifile = new IFile();
  if(ifile->FileExist(datafile)) {
    ifile->open(datafile);
    while(ifile->scanField("Id",idcomp)) {
      ifile->scanField("Pos_0",pos[0]);
      ifile->scanField("Pos_1",pos[1]);
      ifile->scanField("Pos_2",pos[2]);
      ifile->scanField("Beta",beta);
      ifile->scanField("Density",dens);
      // check beta
      if(beta<0) error("Beta must be positive!");
      // center of the Gaussian
      GMM_d_m_.push_back(pos);
      // beta
      GMM_d_beta_.push_back(beta);
      // experimental density
      ovdd_.push_back(dens);
      // new line
      ifile->scanField();
    }
    ifile->close();
  } else {
    error("Cannot find DATA_FILE "+datafile+"\n");
  }
  delete ifile;
  // now create a set from beta (unique set of values)
  set<int> bu(GMM_d_beta_.begin(), GMM_d_beta_.end());
  // now prepare the group vector
  GMM_d_grps_.resize(bu.size());
  // and fill it in
  for(unsigned i=0; i<GMM_d_beta_.size(); ++i) {
    if(GMM_d_beta_[i]>=static_cast<int>(GMM_d_grps_.size())) error("Check Beta values");
    GMM_d_grps_[GMM_d_beta_[i]].push_back(i);
  }
}

void EMMIVOX::calculate_useful_stuff(double reso)
{
  double bfactini = 0.0;
  // if doing Bfactor Monte Carlo
  if(dbfact_>0) {
    // average value of B
    double Bave = 0.0;
    for(unsigned i=0; i<GMM_m_type_.size(); ++i) {
      Bave += GMM_m_s0_[GMM_m_type_[i]];
    }
    Bave /= static_cast<double>(GMM_m_type_.size());
    // initialize B factor to reasonable value based on Gaussian width at half maximum height equal the resolution
    bfactini = 4.0 * ( 2.0 * pow(0.425*pi*reso,2) - Bave );
    // check for min and max
    bfactini = min(bfactmax_, max(bfactmin_, bfactini));
  }
  // set initial Bfactor
  for(map<unsigned,double>::iterator it=GMM_m_b_.begin(); it!=GMM_m_b_.end(); ++it) {
    it->second = bfactini;
  }
  log.printf("  experimental map resolution : %3.2f\n", reso);
  log.printf("  minimum Bfactor value       : %3.2f\n", bfactmin_);
  log.printf("  maximum Bfactor value       : %3.2f\n", bfactmax_);
  log.printf("  initial Bfactor value       : %3.2f\n", bfactini);
}

// prepare auxiliary vectors
void EMMIVOX::get_auxiliary_vectors()
{
// clear lists
  pref_.clear(); invs2_.clear();
// resize
  pref_.resize(GMM_m_res_.size()); invs2_.resize(GMM_m_res_.size());
  #pragma omp parallel num_threads(OpenMP::getNumThreads())
  {
    // private variables
    Vector5d pref, invs2;
    #pragma omp for
    // calculate constant quantities
    for(unsigned im=0; im<GMM_m_res_.size(); ++im) {
      // get atom type
      unsigned atype = GMM_m_type_[im];
      // get residue id
      unsigned ires = GMM_m_res_[im];
      // get bfactor
      double bfact = GMM_m_b_[ires];
      // sigma for 5 gaussians
      Vector5d m_s = GMM_m_s_[atype];
      // calculate constant quantities
      for(unsigned j=0; j<5; ++j) {
        double m_b = m_s[j] + bfact/4.0;
        invs2[j] = 1.0/(GMM_d_s_+inv_pi2_*m_b);
        pref[j]  = cfact_[atype][j] * pow(invs2[j],1.5);
      }
      // put into lists
      pref_[im] = pref;
      invs2_[im]= invs2;
    }
  }
  if(gpu_) push_auxiliary_gpu();
}

void EMMIVOX::push_auxiliary_gpu()
{
#ifdef __PLUMED_HAS_ARRAYFIRE
  // 1) create float version of pref_ and invs2_
  const unsigned GMM_m_size = GMM_m_res_.size();
  vector<float> pref(5*GMM_m_size);
  vector<float> invs2(5*GMM_m_size);
  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for(unsigned i=0; i<GMM_m_size; ++i) {
    for(unsigned j=0; j<5; ++j) {
      pref[GMM_m_size*j+i]  = static_cast<float>(pref_[i][j]);
      invs2[GMM_m_size*j+i] = static_cast<float>(invs2_[i][j]);
    }
  }
  // 2) initialize gpu arrays [GMM_m_size, 5]
  pref_gpu  = af::array(GMM_m_size, 5, &pref.front());
  invs2_gpu = af::array(GMM_m_size, 5, &invs2.front());
#endif
}

void EMMIVOX::doMonteCarloBfact()
{

// iterator over bfactor map
  map<unsigned,double>::iterator it;

// cycle over bfactor map
  for(it=GMM_m_b_.begin(); it!=GMM_m_b_.end(); ++it) {

    // residue id
    unsigned ires = it->first;
    // old bfactor
    double bfactold = it->second;

    // propose move in bfactor
    double bfactnew = bfactold + dbfact_ * ( 2.0 * random_.RandU01() - 1.0 );
    // check boundaries
    if(bfactnew > bfactmax_) {bfactnew = 2.0*bfactmax_ - bfactnew;}
    if(bfactnew < bfactmin_) {bfactnew = 2.0*bfactmin_ - bfactnew;}

    // useful quantities
    map<unsigned, double> deltaov;
    map<unsigned, double>::iterator itov;
    set<unsigned> ngbs;

    #pragma omp parallel num_threads(OpenMP::getNumThreads())
    {
      // private variables
      map<unsigned, double> deltaov_l;
      set<unsigned> ngbs_l;
      #pragma omp for
      // cycle over all the atoms belonging to residue ires
      for(unsigned ia=0; ia<GMM_m_resmap_[ires].size(); ++ia) {

        // get atom id
        unsigned im = GMM_m_resmap_[ires][ia];
        // get atom type
        unsigned atype = GMM_m_type_[im];
        // sigma for 5 Gaussians
        Vector5d m_s = GMM_m_s_[atype];
        // prefactors
        Vector5d cfact = cfact_[atype];
        // and position
        Vector pos = pos_[im];

        // cycle on all the components affected
        for(unsigned i=0; i<GMM_m_nb_[im].size(); ++i) {
          // voxel id
          unsigned id = GMM_m_nb_[im][i];
          // get contribution before change
          double dold=get_overlap(GMM_d_m_[id], pos, GMM_d_s_, cfact, m_s, bfactold);
          // get contribution after change
          double dnew=get_overlap(GMM_d_m_[id], pos, GMM_d_s_, cfact, m_s, bfactnew);
          // update delta overlap
          deltaov_l[id] += dnew-dold;
          // look for neighbors
          for(unsigned j=0; j<GMM_d_nb_[id].size(); ++j) {
            // atom index of potential neighbor
            unsigned in = GMM_d_nb_[id][j];
            // residue index of potential neighbor
            unsigned iresn = GMM_m_res_[in];
            // check if same residue
            if(ires==iresn) continue;
            // distance
            double dist = delta(pos,pos_[in]).modulo();
            // if closer than 0.5 nm, add residue to lists
            if(dist>0 && dist<0.5) ngbs_l.insert(iresn);
          }
        }
      }
      // add to global list
      #pragma omp critical
      for(itov=deltaov_l.begin(); itov!=deltaov_l.end(); ++itov) deltaov[itov->first] += itov->second;
      #pragma omp critical
      ngbs.insert(ngbs_l.begin(), ngbs_l.end());
    }

    // now calculate new and old score
    double old_ene = 0.0;
    double new_ene = 0.0;

    for(itov=deltaov.begin(); itov!=deltaov.end(); ++itov) {
      // id of the component
      unsigned id = itov->first;
      // new value
      double ovmdnew = ovmd_[id]+itov->second;
      // scores
      double devold = scale_*ovmd_[id]+offset_-ovdd_[id];
      double devnew = scale_*ovmdnew+offset_  -ovdd_[id];
      old_ene += -kbt_ * std::log( 0.5 / devold * erf ( devold * inv_sqrt2_ / sigma_min_[GMM_d_beta_[id]] ));
      new_ene += -kbt_ * std::log( 0.5 / devnew * erf ( devnew * inv_sqrt2_ / sigma_min_[GMM_d_beta_[id]] ));
    }

    // add restraint to keep Bfactor of close atoms close
    for(set<unsigned>::iterator is=ngbs.begin(); is!=ngbs.end(); ++is) {
      double gold = (bfactold-GMM_m_b_[*is])/0.2;
      double gnew = (bfactnew-GMM_m_b_[*is])/0.2;
      old_ene += 0.5 * kbt_ * gold * gold;
      new_ene += 0.5 * kbt_ * gnew * gnew;
    }

    // increment number of trials
    MCBtrials_ += 1.0;

    // accept or reject
    bool accept = doAccept(old_ene, new_ene, kbt_);

    // in case of acceptance
    if(accept) {
      // update acceptance rate
      MCBaccept_ += 1.0;
      // update bfactor
      it->second = bfactnew;
      // change all the ovmd_ affected
      for(itov=deltaov.begin(); itov!=deltaov.end(); ++itov) ovmd_[itov->first] += itov->second;
    }

  } // end cycle on bfactors

// update auxiliary lists
  get_auxiliary_vectors();
}

void EMMIVOX::doMonteCarloScale()
{
  // propose move in scale
  double ds = dscale_ * ( 2.0 * random_.RandU01() - 1.0 );
  double new_scale = scale_ + ds;
  // check boundaries
  if(new_scale > scale_max_) {new_scale = 2.0 * scale_max_ - new_scale;}
  if(new_scale < scale_min_) {new_scale = 2.0 * scale_min_ - new_scale;}
  // propose move in offset
  double doff = doffset_ * ( 2.0 * random_.RandU01() - 1.0 );
  double new_off = offset_ + doff;

  // communicate new_scale and new_off to other replicas
  if(!no_aver_ && nrep_>1) {
    if(replica_!=0) {new_scale = 0.0; new_off = 0.0;}
    multi_sim_comm.Sum(&new_scale, 1);
    multi_sim_comm.Sum(&new_off, 1);
  }

  // get new energy
  double new_ene = calculate_Marginal(new_scale, new_off);

  // with marginal, simply multiply by number of replicas!
  if(!no_aver_ && nrep_>1) new_ene *= static_cast<double>(nrep_);

  // accept or reject
  bool accept = doAccept(ene_, new_ene, kbt_);

  // increment number of trials
  MCStrials_ += 1.0;

  // communicate decision
  int do_update = 0;
  if(accept) do_update = 1;
  if(!no_aver_ && nrep_>1) {
    if(replica_!=0) do_update = 0;
    multi_sim_comm.Sum(&do_update, 1);
  }

  // in case of acceptance
  if(do_update) {
    scale_  = new_scale;
    offset_ = new_off;
    MCSaccept_ += 1.0;
    ene_ = new_ene;
  }
}

// get overlap and derivatives
double EMMIVOX::get_overlap_der(const Vector &d_m, const Vector &m_m,
                                const Vector5d &pref, const Vector5d &invs2,
                                Vector &ov_der)
{
  // initialize stuff
  double ov_tot = 0.0;
  // derivative accumulator
  double ov_der_p = 0.0;
  // calculate vector difference with/without pbc
  Vector md = delta(m_m, d_m);
  // norm squared
  double md2 = md[0]*md[0]+md[1]*md[1]+md[2]*md[2];
  // cycle on 5 Gaussians
  for(unsigned j=0; j<5; ++j) {
    // calculate exponent
    double ov = pref[j] * expapprox(-0.5 * md2 * invs2[j]);
    // update derivative prefix
    ov_der_p += ov * invs2[j];
    // increase total overlap
    ov_tot += ov;
  }
  // set derivative
  ov_der = ov_der_p * md;
  // return total overlap
  return ov_tot;
}

// get overlap
double EMMIVOX::get_overlap(const Vector &d_m, const Vector &m_m, double d_s,
                            const Vector5d &cfact, const Vector5d &m_s, double bfact)
{
  // calculate vector difference with/without pbc
  Vector md = delta(m_m, d_m);
  // norm squared
  double md2 = md[0]*md[0]+md[1]*md[1]+md[2]*md[2];
  // cycle on 5 Gaussians
  double ov_tot = 0.0;
  for(unsigned j=0; j<5; ++j) {
    // total value of b
    double m_b = m_s[j]+bfact/4.0;
    // calculate invs2
    double invs2 = 1.0/(d_s+inv_pi2_*m_b);
    // calculate exponent
    double ov = md2*invs2;
    // final calculation
    ov_tot += cfact[j] * pow(invs2, 1.5) * expapprox(-0.5*ov);
  }
  return ov_tot;
}

void EMMIVOX::update_neighbor_sphere()
{
  // dimension of atom vectors
  unsigned GMM_m_size = GMM_m_type_.size();

  // clear global list
  ns_.clear();

  // store reference positions
  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for(unsigned im=0; im<GMM_m_size; ++im) refpos_[im] = pos_[im];

  // cycle on GMM components - in parallel
  #pragma omp parallel num_threads(OpenMP::getNumThreads())
  {
    // private variables
    vector<unsigned> ns_l;
    Vector d_m;
    #pragma omp for
    for(unsigned id=0; id<ovdd_.size(); ++id) {
      // grid point
      d_m = GMM_d_m_[id];
      // cycle on all atoms
      for(unsigned im=0; im<GMM_m_size; ++im) {
        // calculate distance
        double dist = delta(refpos_[im], d_m).modulo();
        // add to local list
        if(dist<=ns_cutoff_) ns_l.push_back(id*GMM_m_size+im);
      }
    }
    // add to global list
    #pragma omp critical
    ns_.insert(ns_.end(), ns_l.begin(), ns_l.end());
  }
}

bool EMMIVOX::do_neighbor_sphere()
{
  vector<double> dist(refpos_.size(), 0.0);
  bool update = false;

// calculate displacement
  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for(unsigned im=0; im<refpos_.size(); ++im) dist[im] = delta(pos_[im],refpos_[im]).modulo();

// check if update or not
  double maxdist = *max_element(dist.begin(), dist.end());
  if(maxdist>=(ns_cutoff_-nl_cutoff_)) update=true;

// return if update or not
  return update;
}

void EMMIVOX::update_neighbor_list()
{
  // dimension of atom vectors
  const unsigned GMM_m_size = GMM_m_type_.size();

  // clear old neighbor list
  nl_.clear();

  // cycle on neighbour sphere - in parallel
  #pragma omp parallel num_threads(OpenMP::getNumThreads())
  {
    // private variables
    vector<unsigned> nl_l;
    unsigned id, im;
    #pragma omp for
    for(unsigned i=0; i<ns_.size(); ++i) {
      // get data (id) and atom (im) indexes
      id = ns_[i] / GMM_m_size;
      im = ns_[i] % GMM_m_size;
      // calculate distance
      double dist = delta(GMM_d_m_[id], pos_[im]).modulo();
      // add to local neighbour list
      if(dist<=nl_cutoff_) nl_l.push_back(ns_[i]);
    }
    // add to global list
    #pragma omp critical
    nl_.insert(nl_.end(), nl_l.begin(), nl_l.end());
  }
  // new dimension of neighbor list
  const unsigned nl_size = nl_.size();
  // now resize derivatives
  ovmd_der_.resize(nl_size);
  // in case of B-factors sampling
  if(dbfact_>0) {
    // now cycle over the neighbor list to creat a list of voxels per atom
    GMM_m_nb_.clear(); GMM_m_nb_.resize(GMM_m_size);
    GMM_d_nb_.clear(); GMM_d_nb_.resize(ovdd_.size());
    for(unsigned i=0; i<nl_size; ++i) {
      unsigned id = nl_[i] / GMM_m_size;
      unsigned im = nl_[i] % GMM_m_size;
      GMM_m_nb_[im].push_back(id);
      GMM_d_nb_[id].push_back(im);
    }
  }
  // in case, transfer data to gpu
  if(gpu_) update_gpu();
}

void EMMIVOX::update_gpu()
{
#ifdef __PLUMED_HAS_ARRAYFIRE
  // constants
  int GMM_m_size = GMM_m_type_.size();
  int nl_size = nl_.size();
  // communicate nl_ to GPU
  af::array nl_gpu = af::array(nl_size, &nl_.front());
  // build nl_id and nl_im arrays on the GPU
  id_gpu = nl_gpu / GMM_m_size;
  im_gpu = nl_gpu % GMM_m_size;
  // now we need to create pref_gpu_nl [nl_size, 5]
  pref_gpu_nl = pref_gpu(im_gpu, af::span);
  // and invs2_gpu_nl [nl_size, 5]
  invs2_gpu_nl = invs2_gpu(im_gpu, af::span);
  // and GMM_d_m_gpu_nl [nl_size, 3]
  GMM_d_m_gpu_nl = GMM_d_m_gpu(id_gpu, af::span);
  // finding sorting arrays
  af::array d_sort;
  af::seq s(nl_size);
  af::array s_arr = s;
  // this will be used for overlap derivatives
  af::sort(ov_k_keys, ov_k_sort, id_gpu, s_arr);
  // this will be used for final derivatives
  af::sort(d_k_keys, d_k_sort, im_gpu, s_arr);
#endif
}

void EMMIVOX::prepare()
{
  if(getExchangeStep()) first_time_=true;
}

// calculate fm+score on cpu
void EMMIVOX::calculate_cpu() {

  // number of atoms
  const unsigned GMM_m_size = GMM_m_type_.size();

  // clear overlap vector
  for(unsigned i=0; i<ovmd_.size(); ++i) ovmd_[i] = 0.0;

  // we have to cycle over all model and data GMM components in the neighbor list
  #pragma omp parallel num_threads(OpenMP::getNumThreads())
  {
    // private variables
    vector<double> v_(ovmd_.size(), 0.0);
    unsigned id, im;
    #pragma omp for
    for(unsigned i=0; i<nl_.size(); ++i) {
      // get data (id) and atom (im) indexes
      id = nl_[i] / GMM_m_size;
      im = nl_[i] % GMM_m_size;
      // add overlap with im component of model GMM
      v_[id] += get_overlap_der(GMM_d_m_[id],pos_[im],pref_[im],invs2_[im],ovmd_der_[i]);
    }
    #pragma omp critical
    for(unsigned i=0; i<ovmd_.size(); ++i) ovmd_[i] += v_[i];
  }

  // average overlap across replicas
  if(!no_aver_ && nrep_>1) {
    double escale = 1.0 / static_cast<double>(nrep_);
    multi_sim_comm.Sum(&ovmd_[0], ovmd_.size());
    for(unsigned i=0; i<ovmd_.size(); ++i) ovmd_[i] *= escale;
  }

  // calculate total energy and get derivatives
  ene_ = calculate_Marginal(scale_, offset_, GMMid_der_);

  // with marginal, simply multiply by number of replicas!
  if(!no_aver_ && nrep_>1) ene_ *= static_cast<double>(nrep_);

  // clear atom derivatives
  for(unsigned i=0; i<atom_der_.size(); ++i) atom_der_[i] = Vector(0,0,0);

  // calculate atom derivatives and virial
  virial_ = Tensor();
  // declare omp reduction for Tensors
  #pragma omp declare reduction( sumTensor : Tensor : omp_out += omp_in )

  #pragma omp parallel num_threads(OpenMP::getNumThreads())
  {
    // private stuff
    vector<Vector> atom_der(atom_der_.size(), Vector(0,0,0));
    unsigned id, im;
    Vector tot_der;
    #pragma omp for reduction (sumTensor : virial_)
    for(unsigned i=0; i<nl_.size(); ++i) {
      // get indexes of data and model component
      id = nl_[i] / GMM_m_size;
      im = nl_[i] % GMM_m_size;
      // chain rule + replica normalization
      tot_der = GMMid_der_[id] * ovmd_der_[i] * scale_;
      // increment atom derivatives
      atom_der[im] += tot_der;
      // and virial
      virial_ += Tensor(pos_[im], -tot_der);
    }
    #pragma omp critical
    for(unsigned i=0; i<atom_der_.size(); ++i) atom_der_[i] += atom_der[i];
  }
}

// calculate fm+score on gpu
void EMMIVOX::calculate_gpu()
{
#ifdef __PLUMED_HAS_ARRAYFIRE
  // number of atoms
  const unsigned GMM_m_size = GMM_m_type_.size();
  // number of data points
  const unsigned GMM_d_size = GMM_d_m_.size();
  // score
  const float scale  = static_cast<float>(scale_);
  const float offset = static_cast<float>(offset_);

  // fill positions in in parallel
  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for (unsigned i=0; i<GMM_m_size; ++i) {
    // fill vectors
    pos_gpu_[i]              = static_cast<float>(pos_[i][0]);
    pos_gpu_[GMM_m_size+i]   = static_cast<float>(pos_[i][1]);
    pos_gpu_[2*GMM_m_size+i] = static_cast<float>(pos_[i][2]);
  }
  // load positions into pos_gpu [GMM_m_size, 3]
  af::array pos_gpu = af::array(GMM_m_size, 3, &pos_gpu_.front());
  // create pos_gpu_nl [nl_size, 3]
  af::array pos_gpu_nl = pos_gpu(im_gpu, af::span);
  // calculate vector difference [nl_size, 3]
  af::array md = GMM_d_m_gpu_nl-pos_gpu_nl;
  // calculate norm squared by column [nl_size, 1]
  af::array md2 = sum(md*md,1);
  // tile to have 5 identical columns [nl_size, 5]
  af::array md2_5 = af::tile(md2, 1, 5);
  // calculate overlaps [nl_size, 5]
  af::array ov = pref_gpu_nl * af::exp(-0.5 * md2_5 * invs2_gpu_nl);
  // and derivatives [nl_size, 5]
  af::array ov_der_nl = invs2_gpu_nl * ov;
  // summation of overlaps over 5 columns [nl_size, 1]
  ov = sum(ov,1);
  // summation of derivatives over 5 rows [nl_size, 1]
  ov_der_nl = sum(ov_der_nl,1);
  // final derivative calculation [nl_size, 3]
  ov_der_nl = af::tile(ov_der_nl, 1, 3) * md;
  // now we have to sum up contributions from the same atom
  // sort by data id
  af::array ov_sort = ov(ov_k_sort);
  // sum equal keys (it works only if adjacent, that's why we sort before)
  af::array ov_k_sum, ov_sum;
  af::sumByKey(ov_k_sum, ov_sum, ov_k_keys, ov_sort);
  // there might be missing data points (i.e. data with no overlaps with atoms!)
  af::array ovmd_gpu = af::constant(0.0, GMM_d_size);
  // assign indexed [GMM_d_size, 1]
  ovmd_gpu(ov_k_sum) = ov_sum(af::span);

  // in case of metainference: average them across replicas
  if(!no_aver_ && nrep_>1) {
    // communicated ovmd_gpu [GMM_d_size, 1] to host
    ovmd_gpu.host(&ovmd_gpu_.front());
    const float escale = 1.0 / static_cast<float>(nrep_);
    multi_sim_comm.Sum(&ovmd_gpu_[0], ovmd_gpu_.size());
    for(unsigned i=0; i<ovmd_gpu_.size(); ++i) ovmd_gpu_[i] *= escale;
    // recommunicate back to device
    ovmd_gpu = af::array(ovmd_gpu_.size(), &ovmd_gpu_.front());
  }

  // calculate score
  // calculate deviation model/data [GMM_d_size, 1]
  af::array dev   = scale * ovmd_gpu + offset - ovdd_gpu;
  // check zero deviation and replace
  af::replace(dev,!(af::iszero(dev)),0.0000001);
  // error function [GMM_d_size, 1]
  af::array errf  = af::erf( dev * inv_sqrt2 * ismin_gpu );
  // energy [GMM_d_size, 1]
  af::array ene   = -kbt * af::log( 0.5 / dev * errf);
  // and derivatives [GMM_d_size, 1]
  af::array d_der = -kbt / errf * sqrt2_pi * af::exp( -0.5 * dev * dev * ismin_gpu * ismin_gpu ) * ismin_gpu + kbt / dev;
  // calculate total energy
  ene = sum(ene);
  // prepare array for derivatives wrt atoms [nl_size, 1]
  af::array der_gpu = d_der(id_gpu);
  // tile [nl_size, 3]
  der_gpu = af::tile(der_gpu, 1, 3);
  // multiple by ov_der_nl and scale [nl_size, 3]
  der_gpu = ov_der_nl * scale * der_gpu;
  // now we have to sum contributions for each atom
  // first, reshuffle der_gpu to sort by atom id [nl_size, 3]
  af::array der_gpu_s = der_gpu(d_k_sort, af::span);
  //  then sum contributions for atom id
  af::array d_k_sum, d_sum;
  af::sumByKey(d_k_sum, d_sum, d_k_keys, der_gpu_s, 0);
  // prepare final derivative vector [GMM_m_size, 3]
  af::array at_der = af::constant(0.0, GMM_m_size, 3);
  // and assign at the right places
  at_der(d_k_sum, af::span) = d_sum(af::span, af::span);

  // FINAL STUFF
  //
  // 1) communicate total energy
  float enef;
  ene.host(&enef);
  ene_ = static_cast<double>(enef);
  // with marginal, simply multiply by number of replicas!
  if(!no_aver_ && nrep_>1) ene_ *= static_cast<double>(nrep_);
  //
  // 2) communicate derivatives
  // der into der_gpu_
  at_der.host(&der_gpu_.front());
  // convert to double vectors
  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for(unsigned i=0; i<GMM_m_size; ++i) {
    atom_der_[i] = Vector(static_cast<double>(der_gpu_[i]),static_cast<double>(der_gpu_[GMM_m_size+i]),static_cast<double>(der_gpu_[2*GMM_m_size+i]));
  }
  //
  // 3) calculate virial
  virial_ = Tensor();
  // declare omp reduction for Tensors
  #pragma omp declare reduction( sumTensor : Tensor : omp_out += omp_in )

  #pragma omp parallel for num_threads(OpenMP::getNumThreads()) reduction (sumTensor : virial_)
  for(unsigned i=0; i<GMM_m_size; ++i) {
    virial_ += Tensor(pos_[i], -atom_der_[i]);
  }
  //
  // 4) communicate overlaps
  // these are needed only in certain situations
  long int step = getStep();
  bool do_comm = false;
  if(ovstride_>0 && step%ovstride_==0)  do_comm = true;
  if(dscale_>0   && step%MCSstride_==0) do_comm = true;
  if(dbfact_>0   && step%MCBstride_==0) do_comm = true;
  if(do_corr_) do_comm = true;
  if(do_comm) {
    // communicate
    ovmd_gpu.host(&ovmd_gpu_.front());
    // convert overlaps to double
    #pragma omp parallel for num_threads(OpenMP::getNumThreads())
    for(unsigned i=0; i<ovmd_.size(); ++i) {
      ovmd_[i] = static_cast<double>(ovmd_gpu_[i]);
    }
  }
#endif
}

void EMMIVOX::calculate()
{
  // store current positions
  pos_ = getPositions();

  // neighbor list update
  if(first_time_ || getExchangeStep() || getStep()%nl_stride_==0) {
    // check if time to update neighbor sphere
    bool update = false;
    if(first_time_ || getExchangeStep()) update = true;
    else update = do_neighbor_sphere();
    // update neighbor sphere
    if(update) update_neighbor_sphere();
    // update neighbor list
    update_neighbor_list();
    first_time_=false;
  }

  // calculate CV
  if(gpu_) calculate_gpu();
  else     calculate_cpu();

  // this part is common for GPUs and CPUs
  // set score, virial, and derivatives
  Value* score = getPntrToComponent("scoreb");
  score->set(ene_);
  setBoxDerivatives(score, virial_);
  #pragma omp parallel for
  for(unsigned i=0; i<atom_der_.size(); ++i) setAtomsDerivatives(score, i, atom_der_[i]);
  // set scale and offset value
  getPntrToComponent("scale")->set(scale_);
  getPntrToComponent("offset")->set(offset_);
  // calculate correlation coefficient
  if(do_corr_) {
    double cc = calculate_corr();
    getPntrToComponent("corr")->set(cc);
  }
  // PRINT STUFF to other files
  // get time step
  long int step = getStep();
  // print status
  if(step%statusstride_==0) print_status(step);
  // write model overlap to file
  if(ovstride_>0 && step%ovstride_==0) write_model_overlap(step);

  // SAMPLE OTHER PARAMETERS
  // Monte Carlo on scale
  if(dscale_>0) {
    // do Monte Carlo
    if(step%MCSstride_==0 && !getExchangeStep()) doMonteCarloScale();
    // calculate acceptance ratio
    double acc = MCSaccept_ / MCStrials_;
    // set acceptance value
    getPntrToComponent("accS")->set(acc);
  }
  // Monte Carlo on bfactors
  if(dbfact_>0) {
    // do Monte Carlo
    if(step%MCBstride_==0 && !getExchangeStep()) doMonteCarloBfact();
    // calculate acceptance ratio
    double acc = MCBaccept_ / MCBtrials_;
    // set value
    getPntrToComponent("accB")->set(acc);
  }

}

double EMMIVOX::calculate_corr()
{
// number of data points
  double nd = static_cast<double>(ovdd_.size());
// average ovmd_ and ovdd_
  double ave_md = std::accumulate(ovmd_.begin(), ovmd_.end(), 0.) / nd;
  double ave_dd = std::accumulate(ovdd_.begin(), ovdd_.end(), 0.) / nd;
// calculate corr
  double num = 0.;
  double den1 = 0.;
  double den2 = 0.;
  #pragma omp parallel for num_threads(OpenMP::getNumThreads()) reduction( + : num, den1, den2)
  for(unsigned i=0; i<ovdd_.size(); ++i) {
    double md = ovmd_[i]-ave_md;
    double dd = ovdd_[i]-ave_dd;
    num  += md*dd;
    den1 += md*md;
    den2 += dd*dd;
  }
// correlation coefficient between ovmd_ and ovdd_
  double cc = num / sqrt(den1*den2);
  return cc;
}

double EMMIVOX::calculate_Marginal(double scale, double offset, vector<double> &GMMid_der)
{
  double ene = 0.0;
  // cycle on all the overlaps
  #pragma omp parallel for num_threads(OpenMP::getNumThreads()) reduction( + : ene)
  for(unsigned id=0; id<ovdd_.size(); ++id) {
    // get ismin
    double ismin = ismin_[id];
    // calculate deviation
    double dev = ( scale * ovmd_[id] + offset - ovdd_[id] );
    // calculate errf
    double errf = erf ( dev * inv_sqrt2_ * ismin );
    // add to  energy
    ene += -kbt_ * std::log( 0.5 / dev * errf );
    // store derivative for later
    GMMid_der[id] = - kbt_/errf*sqrt2_pi_*exp(-0.5*dev*dev*ismin*ismin)*ismin+kbt_/dev;
  }
  // return total energy
  return ene;
}

double EMMIVOX::calculate_Marginal(double scale, double offset)
{
  double ene = 0.0;
  // cycle on all the overlaps
  #pragma omp parallel for num_threads(OpenMP::getNumThreads()) reduction( + : ene)
  for(unsigned id=0; id<ovdd_.size(); ++id) {
    // get ismin
    double ismin = ismin_[id];
    // calculate deviation
    double dev = ( scale * ovmd_[id] + offset - ovdd_[id] );
    // calculate errf
    double errf = erf ( dev * inv_sqrt2_ * ismin );
    // add to  energy
    ene += -kbt_ * std::log( 0.5 / dev * errf );
  }
  // return total energy
  return ene;
}

}
}
