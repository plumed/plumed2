/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2023 The plumed team
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

#ifdef __PLUMED_HAS_LIBTORCH
#include "colvar/Colvar.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Communicator.h"
#include "tools/Matrix.h"
#include "core/GenericMolInfo.h"
#include "core/ActionSet.h"
#include "tools/File.h"
#include "tools/OpenMP.h"
#include <string>
#include <cmath>
#include <map>
#include <numeric>
#include <ctime>
#include "tools/Random.h"

//this avoids nasty warnings on which we have no control
//Here I am using the pragma for gcc/clang, for visual studio you'll need something different
//This should be the correct direction,
//but the warning that it is blocking everithing is:
//  error: ISO C++11 requires at least one argument for the "..." in a variadic macro [-Werror]
//that is not named like [-Werror=shadow]
//so , at time of writing this I did not found any solution but deactivating warning in the makefile

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"

#include <torch/torch.h>
#include <torch/script.h>

#pragma GCC diagnostic push

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

namespace PLMD {
namespace isdb {

//+PLUMEDOC ISDB_COLVAR EMMIVOX
/*
Bayesian single-structure and ensemble refinement with cryo-EM maps.

This action implements the Bayesian approach for single-structure and ensemble refinement from cryo-EM maps introduced <a href="https://doi.org/10.1371/journal.pcbi.1012180">here</a>.
EMMIVox does not require fitting the cryo-EM map with a Gaussian Mixture Model, as done in [EMMI](EMMI.md), but uses directly the voxels in the deposited map.

When run in single-replica mode, this action allows atomistic, flexible refinement (and B-factors inference) of an individual structure into a density map.
A coarse-grained forward model can also be used in combination with the MARTINI force field.
Combined with a multi-replica framework (such as the -multi option in GROMACS), the user can model an ensemble of structures using
the Metainference approach that is introduced in the paper  cited below. The approach can be used to model continous dynamics of flexible regions as well as semi-ordered waters, lipids, and ions.

!!! note "installing libtorch"

    To use EMMIVOX, PLUMED must be linked against the LibTorch library using the instructions [on the module page](module_isdb.md).

## Examples

Complete tutorials for single-structure and ensemble refinement can be found <a href="https://github.com/COSBlab/EMMIVox">here</a>.

*/
//+ENDPLUMEDOC

class EMMIVOX : public Colvar {

private:

// temperature in kbt
  double kbt_;
// model - atom types
  std::vector<unsigned> Model_type_;
// model - list of atom sigmas - one per atom type
  std::vector<Vector5d> Model_s_;
// model - list of atom weights - one per atom type
  std::vector<Vector5d> Model_w_;
// model - map between residue/chain IDs and list of atoms
  std::map< std::pair<unsigned,std::string>, std::vector<unsigned> > Model_resmap_;
// model - list of residue/chain IDs per atom
  std::vector< std::pair<unsigned,std::string> > Model_res_;
// model - list of neighboring voxels per atom
  std::vector< std::vector<unsigned> > Model_nb_;
// model - map between residue/chain ID and bfactor
  std::map< std::pair<unsigned, std::string>, double> Model_b_;
// model - global list of residue/chain IDs
  std::vector< std::pair<unsigned,std::string> > Model_rlist_;
// model density
  std::vector<double> ovmd_;

// data map - voxel position
  std::vector<Vector> Map_m_;
// data map - density
  std::vector<double> ovdd_;
// data map - error
  std::vector<double> exp_err_;

// derivatives
  std::vector<Vector> ovmd_der_;
  std::vector<Vector> atom_der_;
  std::vector<double> score_der_;
// constants
  double inv_sqrt2_, sqrt2_pi_, inv_pi2_;
  std::vector<Vector5d> pref_;
  std::vector<Vector5d> invs2_;
  std::vector<Vector5d> cfact_;
  std::vector<double> cut_;
// metainference
  unsigned nrep_;
  unsigned replica_;
  std::vector<double> ismin_;
// neighbor list
  double nl_dist_cutoff_;
  double nl_gauss_cutoff_;
  unsigned nl_stride_;
  bool first_time_;
  std::vector< std::pair<unsigned,unsigned> > nl_;
  std::vector< std::pair<unsigned,unsigned> > ns_;
  std::vector<Vector> refpos_;
// averaging
  bool no_aver_;
// correlation;
  bool do_corr_;
// Monte Carlo stuff
  Random   random_;
  // Scale and Offset
  double scale_;
  double offset_;
// Bfact Monte Carlo
  double   dbfact_;
  double   bfactmin_;
  double   bfactmax_;
  double   bfactsig_;
  bool     bfactnoc_;
  bool     bfactread_;
  int      MCBstride_;
  double   MCBaccept_;
  double   MCBtrials_;
// residue neighbor list
  std::vector< std::vector<unsigned> > nl_res_;
  bool bfactemin_;
// Martini scattering factors
  bool martini_;
  // status stuff
  unsigned int statusstride_;
  std::string       statusfilename_;
  OFile        statusfile_;
  bool         first_status_;
  // total energy and virial
  double ene_;
  Tensor virial_;
  double eps_;
  // model density file
  unsigned int mapstride_;
  std::string       mapfilename_;
  // Libtorch stuff
  bool gpu_;
  torch::Tensor ovmd_gpu_;
  torch::Tensor ovmd_der_gpu_;
  torch::Tensor ismin_gpu_;
  torch::Tensor ovdd_gpu_;
  torch::Tensor Map_m_gpu_;
  torch::Tensor pref_gpu_;
  torch::Tensor invs2_gpu_;
  torch::Tensor nl_id_gpu_;
  torch::Tensor nl_im_gpu_;
  torch::Tensor pref_nl_gpu_;
  torch::Tensor invs2_nl_gpu_;
  torch::Tensor Map_m_nl_gpu_;
  torch::DeviceType device_t_;
//
// write file with model density
  void write_model_density(long int step);
// get median of vector
  double get_median(std::vector<double> v);
// read and write status
  void read_status();
  void print_status(long int step);
// accept or reject
  bool doAccept(double oldE, double newE, double kbt);
// vector of close residues
  void get_close_residues();
// do MonteCarlo for Bfactor
  void doMonteCarloBfact();
// calculate model parameters
  std::vector<double> get_Model_param(std::vector<AtomNumber> &atoms);
// read data file
  void get_exp_data(const std::string &datafile);
// auxiliary methods
  void prepare_gpu();
  void initialize_Bfactor(double reso);
  void get_auxiliary_vectors();
  void push_auxiliary_gpu();
// calculate overlap between two Gaussians
  double get_overlap(const Vector &d_m, const Vector &m_m,
                     const Vector5d &cfact, const Vector5d &m_s, double bfact);
// update the neighbor list
  void update_neighbor_list();
// update data on device
  void update_gpu();
// update the neighbor sphere
  void update_neighbor_sphere();
  bool do_neighbor_sphere();
// calculate forward model and score on device
  void calculate_fmod();
  void calculate_score();
// calculate correlation
  void calculate_corr();

public:
  static void registerKeywords( Keywords& keys );
  explicit EMMIVOX(const ActionOptions&);
// active methods:
  void prepare() override;
  void calculate() override;
};

PLUMED_REGISTER_ACTION(EMMIVOX,"EMMIVOX")

void EMMIVOX::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOMS","atoms used in the calculation of the density map, typically all heavy atoms");
  keys.add("compulsory","DATA_FILE","file with cryo-EM map");
  keys.add("compulsory","RESOLUTION", "cryo-EM map resolution");
  keys.add("compulsory","NORM_DENSITY","integral of experimental density");
  keys.add("compulsory","WRITE_STRIDE","stride for writing status file");
  keys.add("optional","NL_DIST_CUTOFF","neighbor list distance cutoff");
  keys.add("optional","NL_GAUSS_CUTOFF","neighbor list Gaussian sigma cutoff");
  keys.add("optional","NL_STRIDE","neighbor list update frequency");
  keys.add("optional","SIGMA_MIN","minimum density error");
  keys.add("optional","DBFACT","Bfactor MC step");
  keys.add("optional","BFACT_MAX","Bfactor maximum value");
  keys.add("optional","MCBFACT_STRIDE", "Bfactor MC stride");
  keys.add("optional","BFACT_SIGMA","Bfactor sigma prior");
  keys.add("optional","STATUS_FILE","write a file with all the data useful for restart");
  keys.add("optional","SCALE","scale factor");
  keys.add("optional","OFFSET","offset");
  keys.add("optional","TEMP","temperature");
  keys.add("optional","WRITE_MAP","file with model density");
  keys.add("optional","WRITE_MAP_STRIDE","stride for writing model density to file");
  keys.addFlag("NO_AVER",false,"no ensemble averaging in multi-replica mode");
  keys.addFlag("CORRELATION",false,"calculate correlation coefficient");
  keys.addFlag("GPU",false,"calculate EMMIVOX on GPU with Libtorch");
  keys.addFlag("BFACT_NOCHAIN",false,"Do not use chain ID for Bfactor MC");
  keys.addFlag("BFACT_READ",false,"Read Bfactor on RESTART (automatic with DBFACT>0)");
  keys.addFlag("BFACT_MINIMIZE",false,"Accept only moves that decrease energy");
  keys.addFlag("MARTINI",false,"Use Martini scattering factors");
  keys.addOutputComponent("scoreb","default","scalar","Bayesian score");
  keys.addOutputComponent("scale", "default","scalar","scale factor");
  keys.addOutputComponent("offset","default","scalar","offset");
  keys.addOutputComponent("accB",  "default","scalar", "Bfactor MC acceptance");
  keys.addOutputComponent("kbt",   "default","scalar", "temperature in energy unit");
  keys.addOutputComponent("corr",  "CORRELATION","scalar", "correlation coefficient");
  keys.addDOI("10.1126/sciadv.1501177");
}

EMMIVOX::EMMIVOX(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  nl_dist_cutoff_(1.0), nl_gauss_cutoff_(3.0), nl_stride_(50),
  first_time_(true), no_aver_(false), do_corr_(false),
  scale_(1.), offset_(0.),
  dbfact_(0.0), bfactmin_(0.05), bfactmax_(5.0),
  bfactsig_(0.1), bfactnoc_(false), bfactread_(false),
  MCBstride_(1), MCBaccept_(0.), MCBtrials_(0.), bfactemin_(false),
  martini_(false), statusstride_(0), first_status_(true),
  eps_(0.0001), mapstride_(0), gpu_(false) {
  // set constants
  inv_sqrt2_ = 1.0/sqrt(2.0);
  sqrt2_pi_  = sqrt(2.0 / M_PI);
  inv_pi2_   = 0.5 / M_PI / M_PI;

  // list of atoms
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS", atoms);

  // file with experimental cryo-EM map
  std::string datafile;
  parse("DATA_FILE", datafile);

  // neighbor list cutoffs
  parse("NL_DIST_CUTOFF",nl_dist_cutoff_);
  parse("NL_GAUSS_CUTOFF",nl_gauss_cutoff_);
  // checks
  if(nl_dist_cutoff_<=0. && nl_gauss_cutoff_<=0.) {
    error("You must specify either NL_DIST_CUTOFF or NL_GAUSS_CUTOFF or both");
  }
  if(nl_gauss_cutoff_<=0.) {
    nl_gauss_cutoff_ = 1.0e+10;
  }
  if(nl_dist_cutoff_<=0.) {
    nl_dist_cutoff_ = 1.0e+10;
  }
  // neighbor list update stride
  parse("NL_STRIDE",nl_stride_);
  if(nl_stride_<=0) {
    error("NL_STRIDE must be explicitly specified and positive");
  }

  // minimum value for error
  double sigma_min = 0.2;
  parse("SIGMA_MIN", sigma_min);
  if(sigma_min<0.) {
    error("SIGMA_MIN must be greater or equal to zero");
  }

  // status file parameters
  parse("WRITE_STRIDE", statusstride_);
  if(statusstride_<=0) {
    error("you must specify a positive WRITE_STRIDE");
  }
  parse("STATUS_FILE",  statusfilename_);
  if(statusfilename_=="") {
    statusfilename_ = "EMMIStatus"+getLabel();
  }

  // integral of the experimetal density
  double norm_d;
  parse("NORM_DENSITY", norm_d);

  // temperature
  kbt_ = getkBT();

  // scale and offset
  parse("SCALE", scale_);
  parse("OFFSET",offset_);

  // B-factors MC
  parse("DBFACT",dbfact_);
  // read Bfactors
  parseFlag("BFACT_READ",bfactread_);
  // do not use chains
  parseFlag("BFACT_NOCHAIN",bfactnoc_);
  // other parameters
  if(dbfact_>0.) {
    parse("MCBFACT_STRIDE",MCBstride_);
    parse("BFACT_MAX",bfactmax_);
    parse("BFACT_SIGMA",bfactsig_);
    parseFlag("BFACT_MINIMIZE",bfactemin_);
    // checks
    if(MCBstride_<=0) {
      error("you must specify a positive MCBFACT_STRIDE");
    }
    if(bfactmax_<=bfactmin_) {
      error("you must specify a positive BFACT_MAX");
    }
    if(MCBstride_%nl_stride_!=0) {
      error("MCBFACT_STRIDE must be multiple of NL_STRIDE");
    }
    if(bfactsig_<=0.) {
      error("you must specify a positive BFACT_SIGMA");
    }
  }

  // read map resolution
  double reso;
  parse("RESOLUTION", reso);
  if(reso<=0.) {
    error("RESOLUTION must be strictly positive");
  }

  // averaging or not
  parseFlag("NO_AVER",no_aver_);

  // calculate correlation coefficient
  parseFlag("CORRELATION",do_corr_);

  // write density file
  parse("WRITE_MAP_STRIDE", mapstride_);
  parse("WRITE_MAP", mapfilename_);
  if(mapstride_>0 && mapfilename_=="") {
    error("With WRITE_MAP_STRIDE you must specify WRITE_MAP");
  }

  // use GPU?
  parseFlag("GPU",gpu_);
  // set device
  if (gpu_ && torch::cuda::is_available()) {
    device_t_ = torch::kCUDA;
  } else {
    device_t_ = torch::kCPU;
    gpu_ = false;
  }

// Martini model
  parseFlag("MARTINI",martini_);

  // check read
  checkRead();

  // set parallel stuff
  unsigned mpisize=comm.Get_size();
  if(mpisize>1) {
    error("EMMIVOX supports only OpenMP parallelization");
  }

  // get number of replicas
  if(no_aver_) {
    nrep_ = 1;
  } else {
    nrep_ = multi_sim_comm.Get_size();
  }
  replica_ = multi_sim_comm.Get_rank();

  if(nrep_>1 && dbfact_>0) {
    error("Bfactor sampling not supported with ensemble averaging");
  }

  log.printf("  number of atoms involved : %u\n", atoms.size());
  log.printf("  experimental density map : %s\n", datafile.c_str());
  if(no_aver_) {
    log.printf("  without ensemble averaging\n");
  }
  if(gpu_) {
    log.printf("  running on GPU \n");
  } else {
    log.printf("  running on CPU \n");
  }
  if(nl_dist_cutoff_ <1.0e+10) {
    log.printf("  neighbor list distance cutoff : %lf\n", nl_dist_cutoff_);
  }
  if(nl_gauss_cutoff_<1.0e+10) {
    log.printf("  neighbor list Gaussian sigma cutoff : %lf\n", nl_gauss_cutoff_);
  }
  log.printf("  neighbor list update stride : %u\n",  nl_stride_);
  log.printf("  minimum density error : %f\n", sigma_min);
  log.printf("  scale factor : %lf\n", scale_);
  log.printf("  offset : %lf\n", offset_);
  log.printf("  reading/writing to status file : %s\n", statusfilename_.c_str());
  log.printf("  with stride : %u\n", statusstride_);
  if(dbfact_>0) {
    log.printf("  maximum Bfactor MC move : %f\n", dbfact_);
    log.printf("  stride MC move : %u\n", MCBstride_);
    log.printf("  using prior with sigma : %f\n", bfactsig_);
  }
  if(bfactread_) {
    log.printf("  reading Bfactors from file : %s\n", statusfilename_.c_str());
  }
  log.printf("  temperature of the system in energy unit : %f\n", kbt_);
  if(nrep_>1) {
    log.printf("  number of replicas for averaging: %u\n", nrep_);
    log.printf("  replica ID : %u\n", replica_);
  }
  if(mapstride_>0) {
    log.printf("  writing model density to file : %s\n", mapfilename_.c_str());
    log.printf("  with stride : %u\n", mapstride_);
  }
  if(martini_) {
    log.printf("  using Martini scattering factors\n");
  }

  // calculate model constant parameters
  std::vector<double> Model_w = get_Model_param(atoms);

  // read experimental map and errors
  get_exp_data(datafile);
  log.printf("  number of voxels : %u\n", static_cast<unsigned>(ovdd_.size()));

  // normalize atom weight map
  double norm_m = accumulate(Model_w.begin(),  Model_w.end(),  0.0);
  // renormalization and constant factor on atom types
  for(unsigned i=0; i<Model_w_.size(); ++i) {
    Vector5d cf;
    for(unsigned j=0; j<5; ++j) {
      Model_w_[i][j] *= norm_d / norm_m;
      cf[j] = Model_w_[i][j]/pow( 2.0*pi, 1.5 );
    }
    cfact_.push_back(cf);
  }

  // median density
  double ovdd_m = get_median(ovdd_);
  // median experimental error
  double err_m  = get_median(exp_err_);
  // minimum error
  double minerr = sigma_min*ovdd_m;
  // print out statistics
  log.printf("     median density : %lf\n", ovdd_m);
  log.printf("     minimum error  : %lf\n", minerr);
  log.printf("     median error   : %lf\n", err_m);
  // populate ismin: cycle on all voxels
  for(unsigned id=0; id<ovdd_.size(); ++id) {
    // define smin
    double smin = std::max(minerr, exp_err_[id]);
    // and to ismin_
    ismin_.push_back(1.0/smin);
  }

  // prepare gpu stuff: map centers, data, and error
  prepare_gpu();

  // initialize Bfactors
  initialize_Bfactor(reso);

  // read status file if restarting
  if(getRestart() || bfactread_) {
    read_status();
  }

  // prepare auxiliary vectors
  get_auxiliary_vectors();

  // prepare other vectors: data and derivatives
  ovmd_.resize(ovdd_.size());
  atom_der_.resize(Model_type_.size());
  score_der_.resize(ovdd_.size());

  // add components
  addComponentWithDerivatives("scoreb");
  componentIsNotPeriodic("scoreb");
  addComponent("scale");
  componentIsNotPeriodic("scale");
  addComponent("offset");
  componentIsNotPeriodic("offset");
  addComponent("kbt");
  componentIsNotPeriodic("kbt");
  if(dbfact_>0)   {
    addComponent("accB");
    componentIsNotPeriodic("accB");
  }
  if(do_corr_)    {
    addComponent("corr");
    componentIsNotPeriodic("corr");
  }

  // initialize random seed
  unsigned iseed = time(NULL)+replica_;
  random_.setSeed(-iseed);

  // request atoms
  requestAtoms(atoms);

  // print bibliography
  log<<"  Bibliography "<<plumed.cite("Bonomi, Camilloni, Bioinformatics, 33, 3999 (2017)");
  log<<plumed.cite("Hoff, Thomasen, Lindorff-Larsen, Bonomi, PLoS Comput. Biol. 20 (2024) e1012180");
  if(!no_aver_ && nrep_>1) {
    log<<plumed.cite("Bonomi, Camilloni, Cavalli, Vendruscolo, Sci. Adv. 2, e150117 (2016)");
  }
  log<<"\n";
}

void EMMIVOX::prepare_gpu() {
  // number of data points
  int nd = ovdd_.size();
  // 1) put ismin_ on device_t_
  ismin_gpu_ = torch::from_blob(ismin_.data(), {nd}, torch::kFloat64).to(torch::kFloat32).to(device_t_);
  // 2) put ovdd_ on device_t_
  ovdd_gpu_  = torch::from_blob(ovdd_.data(),  {nd}, torch::kFloat64).to(torch::kFloat32).to(device_t_);
  // 3) put Map_m_ on device_t_
  std::vector<double> Map_m_gpu(3*nd);
  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for(int i=0; i<nd; ++i) {
    Map_m_gpu[i]      = Map_m_[i][0];
    Map_m_gpu[i+nd]   = Map_m_[i][1];
    Map_m_gpu[i+2*nd] = Map_m_[i][2];
  }
  // libtorch tensor
  Map_m_gpu_ = torch::from_blob(Map_m_gpu.data(), {3,nd}, torch::kFloat64).clone().to(torch::kFloat32).to(device_t_);
}

void EMMIVOX::write_model_density(long int step) {
  OFile ovfile;
  ovfile.link(*this);
  std::string num;
  Tools::convert(step,num);
  std::string name = mapfilename_+"-"+num;
  ovfile.open(name);
  ovfile.setHeavyFlush();
  ovfile.fmtField("%10.7e ");
// write density
  for(unsigned i=0; i<ovmd_.size(); ++i) {
    ovfile.printField("Model", ovmd_[i]);
    ovfile.printField("ModelScaled", scale_ * ovmd_[i] + offset_);
    ovfile.printField("Data", ovdd_[i]);
    ovfile.printField();
  }
  ovfile.close();
}

double EMMIVOX::get_median(std::vector<double> v) {
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

void EMMIVOX::read_status() {
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
      // read bfactors if doing fitting of reading it at restart
      if(dbfact_>0 || bfactread_) {
        // cycle on residues
        for(unsigned ir=0; ir<Model_rlist_.size(); ++ir) {
          // key: pair of residue/chain IDs
          std::pair<unsigned,std::string> key = Model_rlist_[ir];
          // convert ires to std::string
          std::string num;
          Tools::convert(key.first,num);
          // read entry
          std::string ch = key.second;
          if(ch==" ") {
            ch="";
          }
          ifile->scanField("bf-"+num+":"+ch, Model_b_[key]);
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

void EMMIVOX::print_status(long int step) {
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
  // write bfactors only if doing fitting or reading bfactors
  if(dbfact_>0 || bfactread_) {
    // cycle on residues
    for(unsigned ir=0; ir<Model_rlist_.size(); ++ir) {
      // key: pair of residue/chain IDs
      std::pair<unsigned,std::string> key = Model_rlist_[ir];
      // bfactor from map
      double bf = Model_b_[key];
      // convert ires to std::string
      std::string num;
      Tools::convert(key.first,num);
      // print entry
      statusfile_.printField("bf-"+num+":"+key.second, bf);
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
    if( s < exp(-delta) ) {
      accept = true;
    }
  }
  return accept;
}

std::vector<double> EMMIVOX::get_Model_param(std::vector<AtomNumber> &atoms) {
  // check if MOLINFO line is present
  auto* moldat=plumed.getActionSet().selectLatest<GenericMolInfo*>(this);
  if(!moldat) {
    error("MOLINFO DATA not found\n");
  }
  log<<"  MOLINFO DATA found with label " <<moldat->getLabel()<<", using proper atom names\n";

  // list of weights - one per atom
  std::vector<double> Model_w;
  // 5-Gaussians parameters
  // map of atom types to A and B coefficients of scattering factor
  // f(s) = A * exp(-B*s**2)
  // B is in Angstrom squared
  // Elastic atomic scattering factors of electrons for neutral atoms
  // and s up to 6.0 A^-1: as implemented in PLUMED
  // map between an atom type and an index
  std::map<std::string, unsigned> type_map;
  // atomistic types
  type_map["C"]=0;
  type_map["O"]=1;
  type_map["N"]=2;
  type_map["S"]=3;
  type_map["P"]=4;
  type_map["F"]=5;
  type_map["NA"]=6;
  type_map["MG"]=7;
  type_map["CL"]=8;
  type_map["CA"]=9;
  type_map["K"]=10;
  type_map["ZN"]=11;
  // Martini types
  type_map["ALA_BB"]=12;
  type_map["ALA_SC1"]=13;
  type_map["CYS_BB"]=14;
  type_map["CYS_SC1"]=15;
  type_map["ASP_BB"]=16;
  type_map["ASP_SC1"]=17;
  type_map["GLU_BB"]=18;
  type_map["GLU_SC1"]=19;
  type_map["PHE_BB"]=20;
  type_map["PHE_SC1"]=21;
  type_map["PHE_SC2"]=22;
  type_map["PHE_SC3"]=23;
  type_map["GLY_BB"]=24;
  type_map["HIS_BB"]=25;
  type_map["HIS_SC1"]=26;
  type_map["HIS_SC2"]=27;
  type_map["HIS_SC3"]=28;
  type_map["ILE_BB"]=29;
  type_map["ILE_SC1"]=30;
  type_map["LYS_BB"]=31;
  type_map["LYS_SC1"]=32;
  type_map["LYS_SC2"]=33;
  type_map["LEU_BB"]=34;
  type_map["LEU_SC1"]=35;
  type_map["MET_BB"]=36;
  type_map["MET_SC1"]=37;
  type_map["ASN_BB"]=38;
  type_map["ASN_SC1"]=39;
  type_map["PRO_BB"]=40;
  type_map["PRO_SC1"]=41;
  type_map["GLN_BB"]=42;
  type_map["GLN_SC1"]=43;
  type_map["ARG_BB"]=44;
  type_map["ARG_SC1"]=45;
  type_map["ARG_SC2"]=46;
  type_map["SER_BB"]=47;
  type_map["SER_SC1"]=48;
  type_map["THR_BB"]=49;
  type_map["THR_SC1"]=50;
  type_map["VAL_BB"]=51;
  type_map["VAL_SC1"]=52;
  type_map["TRP_BB"]=53;
  type_map["TRP_SC1"]=54;
  type_map["TRP_SC2"]=55;
  type_map["TRP_SC3"]=56;
  type_map["TRP_SC4"]=57;
  type_map["TRP_SC5"]=58;
  type_map["TYR_BB"]=59;
  type_map["TYR_SC1"]=60;
  type_map["TYR_SC2"]=61;
  type_map["TYR_SC3"]=62;
  type_map["TYR_SC4"]=63;
  // fill in sigma vector for atoms
  Model_s_.push_back(0.01*Vector5d(0.114,1.0825,5.4281,17.8811,51.1341));   // C
  Model_s_.push_back(0.01*Vector5d(0.0652,0.6184,2.9449,9.6298,28.2194));   // O
  Model_s_.push_back(0.01*Vector5d(0.0541,0.5165,2.8207,10.6297,34.3764));  // N
  Model_s_.push_back(0.01*Vector5d(0.0838,0.7788,4.3462,15.5846,44.63655)); // S
  Model_s_.push_back(0.01*Vector5d(0.0977,0.9084,4.9654,18.5471,54.3648));  // P
  Model_s_.push_back(0.01*Vector5d(0.0613,0.5753,2.6858,8.8214,25.6668));   // F
  Model_s_.push_back(0.01*Vector5d(0.1684,1.7150,8.8386,50.8265,147.2073)); // NA
  Model_s_.push_back(0.01*Vector5d(0.1356,1.3579,6.9255,32.3165,92.1138));  // MG
  Model_s_.push_back(0.01*Vector5d(0.0694,0.6443,3.5351,12.5058,35.8633));  // CL
  Model_s_.push_back(0.01*Vector5d(0.1742,1.8329,8.8407,47.4583,134.9613)); // CA
  Model_s_.push_back(0.01*Vector5d(0.1660,1.6906,8.7447,46.7825,165.6923)); // K
  Model_s_.push_back(0.01*Vector5d(0.0876,0.8650,3.8612,18.8726,64.7016));  // ZN
  // fill in sigma vector for Martini beads
  Model_s_.push_back(0.01*Vector5d(22.000000,22.000000,22.000000,22.000000,22.000000)); // ALA_BB
  Model_s_.push_back(0.01*Vector5d(0.500000,0.500000,0.500000,0.500000,0.500000)); // ALA_SC1
  Model_s_.push_back(0.01*Vector5d(23.000000,23.000000,23.000000,23.000000,23.000000)); // CYS_BB
  Model_s_.push_back(0.01*Vector5d(8.500000,8.500000,8.500000,8.500000,8.500000)); // CYS_SC1
  Model_s_.push_back(0.01*Vector5d(23.000000,23.000000,23.000000,23.000000,23.000000)); // ASP_BB
  Model_s_.push_back(0.01*Vector5d(17.000000,17.000000,17.000000,17.000000,17.000000)); // ASP_SC1
  Model_s_.push_back(0.01*Vector5d(22.000000,22.000000,22.000000,22.000000,22.000000)); // GLU_BB
  Model_s_.push_back(0.01*Vector5d(24.000000,24.000000,24.000000,24.000000,24.000000)); // GLU_SC1
  Model_s_.push_back(0.01*Vector5d(23.000000,23.000000,23.000000,23.000000,23.000000)); // PHE_BB
  Model_s_.push_back(0.01*Vector5d(17.500000,17.500000,17.500000,17.500000,17.500000)); // PHE_SC1
  Model_s_.push_back(0.01*Vector5d(11.000000,11.000000,11.000000,11.000000,11.000000)); // PHE_SC2
  Model_s_.push_back(0.01*Vector5d(11.000000,11.000000,11.000000,11.000000,11.000000)); // PHE_SC3
  Model_s_.push_back(0.01*Vector5d(23.000000,23.000000,23.000000,23.000000,23.000000)); // GLY_BB
  Model_s_.push_back(0.01*Vector5d(23.000000,23.000000,23.000000,23.000000,23.000000)); // HIS_BB
  Model_s_.push_back(0.01*Vector5d(11.500000,11.500000,11.500000,11.500000,11.500000)); // HIS_SC1
  Model_s_.push_back(0.01*Vector5d(9.000000,9.000000,9.000000,9.000000,9.000000)); // HIS_SC2
  Model_s_.push_back(0.01*Vector5d(8.500000,8.500000,8.500000,8.500000,8.500000)); // HIS_SC3
  Model_s_.push_back(0.01*Vector5d(23.000000,23.000000,23.000000,23.000000,23.000000)); // ILE_BB
  Model_s_.push_back(0.01*Vector5d(25.500000,25.500000,25.500000,25.500000,25.500000)); // ILE_SC1
  Model_s_.push_back(0.01*Vector5d(22.000000,22.000000,22.000000,22.000000,22.000000)); // LYS_BB
  Model_s_.push_back(0.01*Vector5d(18.000000,18.000000,18.000000,18.000000,18.000000)); // LYS_SC1
  Model_s_.push_back(0.01*Vector5d(11.000000,11.000000,11.000000,11.000000,11.000000)); // LYS_SC2
  Model_s_.push_back(0.01*Vector5d(22.000000,22.000000,22.000000,22.000000,22.000000)); // LEU_BB
  Model_s_.push_back(0.01*Vector5d(21.500000,21.500000,21.500000,21.500000,21.500000)); // LEU_SC1
  Model_s_.push_back(0.01*Vector5d(22.000000,22.000000,22.000000,22.000000,22.000000)); // MET_BB
  Model_s_.push_back(0.01*Vector5d(22.500000,22.500000,22.500000,22.500000,22.500000)); // MET_SC1
  Model_s_.push_back(0.01*Vector5d(22.000000,22.000000,22.000000,22.000000,22.000000)); // ASN_BB
  Model_s_.push_back(0.01*Vector5d(18.500000,18.500000,18.500000,18.500000,18.500000)); // ASN_SC1
  Model_s_.push_back(0.01*Vector5d(23.500000,23.500000,23.500000,23.500000,23.500000)); // PRO_BB
  Model_s_.push_back(0.01*Vector5d(17.500000,17.500000,17.500000,17.500000,17.500000)); // PRO_SC1
  Model_s_.push_back(0.01*Vector5d(22.000000,22.000000,22.000000,22.000000,22.000000)); // GLN_BB
  Model_s_.push_back(0.01*Vector5d(24.500000,24.500000,24.500000,24.500000,24.500000)); // GLN_SC1
  Model_s_.push_back(0.01*Vector5d(23.000000,23.000000,23.000000,23.000000,23.000000)); // ARG_BB
  Model_s_.push_back(0.01*Vector5d(18.000000,18.000000,18.000000,18.000000,18.000000)); // ARG_SC1
  Model_s_.push_back(0.01*Vector5d(18.000000,18.000000,18.000000,18.000000,18.000000)); // ARG_SC2
  Model_s_.push_back(0.01*Vector5d(23.000000,23.000000,23.000000,23.000000,23.000000)); // SER_BB
  Model_s_.push_back(0.01*Vector5d(9.000000,9.000000,9.000000,9.000000,9.000000)); // SER_SC1
  Model_s_.push_back(0.01*Vector5d(23.000000,23.000000,23.000000,23.000000,23.000000)); // THR_BB
  Model_s_.push_back(0.01*Vector5d(17.000000,17.000000,17.000000,17.000000,17.000000)); // THR_SC1
  Model_s_.push_back(0.01*Vector5d(23.000000,23.000000,23.000000,23.000000,23.000000)); // VAL_BB
  Model_s_.push_back(0.01*Vector5d(18.000000,18.000000,18.000000,18.000000,18.000000)); // VAL_SC1
  Model_s_.push_back(0.01*Vector5d(23.000000,23.000000,23.000000,23.000000,23.000000)); // TRP_BB
  Model_s_.push_back(0.01*Vector5d(11.500000,11.500000,11.500000,11.500000,11.500000)); // TRP_SC1
  Model_s_.push_back(0.01*Vector5d(9.000000,9.000000,9.000000,9.000000,9.000000)); // TRP_SC2
  Model_s_.push_back(0.01*Vector5d(11.000000,11.000000,11.000000,11.000000,11.000000)); // TRP_SC3
  Model_s_.push_back(0.01*Vector5d(11.000000,11.000000,11.000000,11.000000,11.000000)); // TRP_SC4
  Model_s_.push_back(0.01*Vector5d(9.500000,9.500000,9.500000,9.500000,9.500000)); // TRP_SC5
  Model_s_.push_back(0.01*Vector5d(23.000000,23.000000,23.000000,23.000000,23.000000)); // TYR_BB
  Model_s_.push_back(0.01*Vector5d(12.000000,12.000000,12.000000,12.000000,12.000000)); // TYR_SC1
  Model_s_.push_back(0.01*Vector5d(11.000000,11.000000,11.000000,11.000000,11.000000)); // TYR_SC2
  Model_s_.push_back(0.01*Vector5d(11.000000,11.000000,11.000000,11.000000,11.000000)); // TYR_SC3
  Model_s_.push_back(0.01*Vector5d(8.500000,8.500000,8.500000,8.500000,8.500000)); // TYR_SC4
  // fill in weight vector for atoms
  Model_w_.push_back(Vector5d(0.0489,0.2091,0.7537,1.1420,0.3555)); // C
  Model_w_.push_back(Vector5d(0.0365,0.1729,0.5805,0.8814,0.3121)); // O
  Model_w_.push_back(Vector5d(0.0267,0.1328,0.5301,1.1020,0.4215)); // N
  Model_w_.push_back(Vector5d(0.0915,0.4312,1.0847,2.4671,1.0852)); // S
  Model_w_.push_back(Vector5d(0.1005,0.4615,1.0663,2.5854,1.2725)); // P
  Model_w_.push_back(Vector5d(0.0382,0.1822,0.5972,0.7707,0.2130)); // F
  Model_w_.push_back(Vector5d(0.1260,0.6442,0.8893,1.8197,1.2988)); // NA
  Model_w_.push_back(Vector5d(0.1130,0.5575,0.9046,2.1580,1.4735)); // MG
  Model_w_.push_back(Vector5d(0.0799,0.3891,1.0037,2.3332,1.0507)); // CL
  Model_w_.push_back(Vector5d(0.2355,0.9916,2.3959,3.7252,2.5647)); // CA
  Model_w_.push_back(Vector5d(0.2149,0.8703,2.4999,2.3591,3.0318)); // K
  Model_w_.push_back(Vector5d(0.1780,0.8096,1.6744,1.9499,1.4495)); // ZN
  // fill in weight vector for Martini beads
  Model_w_.push_back(Vector5d(1.800000,1.800000,1.800000,1.800000,1.800000)); // ALA_BB
  Model_w_.push_back(Vector5d(0.100000,0.100000,0.100000,0.100000,0.100000)); // ALA_SC1
  Model_w_.push_back(Vector5d(1.900000,1.900000,1.900000,1.900000,1.900000)); // CYS_BB
  Model_w_.push_back(Vector5d(1.100000,1.100000,1.100000,1.100000,1.100000)); // CYS_SC1
  Model_w_.push_back(Vector5d(1.900000,1.900000,1.900000,1.900000,1.900000)); // ASP_BB
  Model_w_.push_back(Vector5d(1.700000,1.700000,1.700000,1.700000,1.700000)); // ASP_SC1
  Model_w_.push_back(Vector5d(1.800000,1.800000,1.800000,1.800000,1.800000)); // GLU_BB
  Model_w_.push_back(Vector5d(2.300000,2.300000,2.300000,2.300000,2.300000)); // GLU_SC1
  Model_w_.push_back(Vector5d(1.900000,1.900000,1.900000,1.900000,1.900000)); // PHE_BB
  Model_w_.push_back(Vector5d(1.400000,1.400000,1.400000,1.400000,1.400000)); // PHE_SC1
  Model_w_.push_back(Vector5d(0.900000,0.900000,0.900000,0.900000,0.900000)); // PHE_SC2
  Model_w_.push_back(Vector5d(0.900000,0.900000,0.900000,0.900000,0.900000)); // PHE_SC3
  Model_w_.push_back(Vector5d(1.900000,1.900000,1.900000,1.900000,1.900000)); // GLY_BB
  Model_w_.push_back(Vector5d(1.900000,1.900000,1.900000,1.900000,1.900000)); // HIS_BB
  Model_w_.push_back(Vector5d(0.900000,0.900000,0.900000,0.900000,0.900000)); // HIS_SC1
  Model_w_.push_back(Vector5d(0.800000,0.800000,0.800000,0.800000,0.800000)); // HIS_SC2
  Model_w_.push_back(Vector5d(0.800000,0.800000,0.800000,0.800000,0.800000)); // HIS_SC3
  Model_w_.push_back(Vector5d(1.900000,1.900000,1.900000,1.900000,1.900000)); // ILE_BB
  Model_w_.push_back(Vector5d(2.000000,2.000000,2.000000,2.000000,2.000000)); // ILE_SC1
  Model_w_.push_back(Vector5d(1.800000,1.800000,1.800000,1.800000,1.800000)); // LYS_BB
  Model_w_.push_back(Vector5d(1.400000,1.400000,1.400000,1.400000,1.400000)); // LYS_SC1
  Model_w_.push_back(Vector5d(0.900000,0.900000,0.900000,0.900000,0.900000)); // LYS_SC2
  Model_w_.push_back(Vector5d(1.800000,1.800000,1.800000,1.800000,1.800000)); // LEU_BB
  Model_w_.push_back(Vector5d(1.900000,1.900000,1.900000,1.900000,1.900000)); // LEU_SC1
  Model_w_.push_back(Vector5d(1.800000,1.800000,1.800000,1.800000,1.800000)); // MET_BB
  Model_w_.push_back(Vector5d(2.300000,2.300000,2.300000,2.300000,2.300000)); // MET_SC1
  Model_w_.push_back(Vector5d(1.800000,1.800000,1.800000,1.800000,1.800000)); // ASN_BB
  Model_w_.push_back(Vector5d(1.800000,1.800000,1.800000,1.800000,1.800000)); // ASN_SC1
  Model_w_.push_back(Vector5d(1.900000,1.900000,1.900000,1.900000,1.900000)); // PRO_BB
  Model_w_.push_back(Vector5d(1.400000,1.400000,1.400000,1.400000,1.400000)); // PRO_SC1
  Model_w_.push_back(Vector5d(1.800000,1.800000,1.800000,1.800000,1.800000)); // GLN_BB
  Model_w_.push_back(Vector5d(2.300000,2.300000,2.300000,2.300000,2.300000)); // GLN_SC1
  Model_w_.push_back(Vector5d(1.900000,1.900000,1.900000,1.900000,1.900000)); // ARG_BB
  Model_w_.push_back(Vector5d(1.400000,1.400000,1.400000,1.400000,1.400000)); // ARG_SC1
  Model_w_.push_back(Vector5d(1.800000,1.800000,1.800000,1.800000,1.800000)); // ARG_SC2
  Model_w_.push_back(Vector5d(1.900000,1.900000,1.900000,1.900000,1.900000)); // SER_BB
  Model_w_.push_back(Vector5d(0.800000,0.800000,0.800000,0.800000,0.800000)); // SER_SC1
  Model_w_.push_back(Vector5d(1.900000,1.900000,1.900000,1.900000,1.900000)); // THR_BB
  Model_w_.push_back(Vector5d(1.400000,1.400000,1.400000,1.400000,1.400000)); // THR_SC1
  Model_w_.push_back(Vector5d(1.900000,1.900000,1.900000,1.900000,1.900000)); // VAL_BB
  Model_w_.push_back(Vector5d(1.400000,1.400000,1.400000,1.400000,1.400000)); // VAL_SC1
  Model_w_.push_back(Vector5d(1.900000,1.900000,1.900000,1.900000,1.900000)); // TRP_BB
  Model_w_.push_back(Vector5d(0.900000,0.900000,0.900000,0.900000,0.900000)); // TRP_SC1
  Model_w_.push_back(Vector5d(0.800000,0.800000,0.800000,0.800000,0.800000)); // TRP_SC2
  Model_w_.push_back(Vector5d(0.900000,0.900000,0.900000,0.900000,0.900000)); // TRP_SC3
  Model_w_.push_back(Vector5d(0.900000,0.900000,0.900000,0.900000,0.900000)); // TRP_SC4
  Model_w_.push_back(Vector5d(0.800000,0.800000,0.800000,0.800000,0.800000)); // TRP_SC5
  Model_w_.push_back(Vector5d(1.900000,1.900000,1.900000,1.900000,1.900000)); // TYR_BB
  Model_w_.push_back(Vector5d(0.900000,0.900000,0.900000,0.900000,0.900000)); // TYR_SC1
  Model_w_.push_back(Vector5d(0.900000,0.900000,0.900000,0.900000,0.900000)); // TYR_SC2
  Model_w_.push_back(Vector5d(0.900000,0.900000,0.900000,0.900000,0.900000)); // TYR_SC3
  Model_w_.push_back(Vector5d(0.800000,0.800000,0.800000,0.800000,0.800000)); // TYR_SC4
  // cycle on atoms
  for(unsigned i=0; i<atoms.size(); ++i) {
    // get atom name
    std::string name = moldat->getAtomName(atoms[i]);
    // get residue name
    std::string resname = moldat->getResidueName(atoms[i]);
    // type of atoms/bead
    std::string type_s;
    // Martini model
    if(martini_) {
      type_s = resname+"_"+name;
      // Atomistic model
    } else {
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
      // convert to std::string
      type_s = std::string(1,type);
      // special cases
      if(name=="SOD" || name=="NA" || name =="Na") {
        type_s = "NA";
      }
      if(name=="MG"  || name=="Mg") {
        type_s = "MG";
      }
      if(name=="CLA" || name=="CL" || name =="Cl") {
        type_s = "CL";
      }
      if((resname=="CAL" || resname=="CA") && (name=="CAL" || name=="CA" || name =="C0")) {
        type_s = "CA";
      }
      if(name=="POT" || name=="K") {
        type_s = "K";
      }
      if(name=="ZN"  || name=="Zn") {
        type_s = "ZN";
      }
    }
    // check if key in map
    if(type_map.find(type_s) != type_map.end()) {
      // save atom type
      Model_type_.push_back(type_map[type_s]);
      // this will be normalized in the final density
      Vector5d w = Model_w_[type_map[type_s]];
      Model_w.push_back(w[0]+w[1]+w[2]+w[3]+w[4]);
      // get residue id
      unsigned ires = moldat->getResidueNumber(atoms[i]);
      // and chain
      std::string c ("*");
      if(!bfactnoc_) {
        c = moldat->getChainID(atoms[i]);
      }
      // define pair residue/chain IDs
      std::pair<unsigned,std::string> key = std::make_pair(ires,c);
      // add to map between residue/chain and list of atoms
      Model_resmap_[key].push_back(i);
      // and global list of residue/chain per atom
      Model_res_.push_back(key);
      // initialize Bfactor map
      Model_b_[key] = 0.0;
    } else {
      error("Wrong atom type "+type_s+" from atom name "+name+"\n");
    }
  }
  // create ordered vector of residue-chain IDs
  for(unsigned i=0; i<Model_res_.size(); ++i) {
    std::pair<unsigned,std::string> key = Model_res_[i];
    // search in Model_rlist_
    if(find(Model_rlist_.begin(), Model_rlist_.end(), key) == Model_rlist_.end()) {
      Model_rlist_.push_back(key);
    }
  }
  // return weights
  return Model_w;
}

// read experimental data file in PLUMED format:
void EMMIVOX::get_exp_data(const std::string &datafile) {
  Vector pos;
  double dens, err;
  int idcomp;

// open file
  IFile *ifile = new IFile();
  if(ifile->FileExist(datafile)) {
    ifile->open(datafile);
    while(ifile->scanField("Id",idcomp)) {
      ifile->scanField("Pos_0",pos[0]);
      ifile->scanField("Pos_1",pos[1]);
      ifile->scanField("Pos_2",pos[2]);
      ifile->scanField("Density",dens);
      ifile->scanField("Error",err);
      // voxel center
      Map_m_.push_back(pos);
      // experimental density
      ovdd_.push_back(dens);
      // error
      exp_err_.push_back(err);
      // new line
      ifile->scanField();
    }
    ifile->close();
  } else {
    error("Cannot find DATA_FILE "+datafile+"\n");
  }
  delete ifile;
}

void EMMIVOX::initialize_Bfactor(double reso) {
  double bfactini = 0.0;
  // if doing Bfactor Monte Carlo
  if(dbfact_>0) {
    // initialize B factor based on empirical relation between resolution and average bfactor
    // calculated on ~8000 cryo-EM data with resolution < 5 Ang
    // Bfact = A*reso**2+B; with A=6.95408 B=-2.45697/100.0 nm^2
    bfactini = 6.95408*reso*reso - 0.01*2.45697;
    // check for min and max
    bfactini = std::min(bfactmax_, std::max(bfactmin_, bfactini));
  }
  // set initial Bfactor
  for(std::map< std::pair<unsigned,std::string>, double>::iterator it=Model_b_.begin(); it!=Model_b_.end(); ++it) {
    it->second = bfactini;
  }
  log.printf("  experimental map resolution : %3.2f\n", reso);
  // if doing Bfactor Monte Carlo
  if(dbfact_>0) {
    log.printf("  minimum Bfactor value : %3.2f\n", bfactmin_);
    log.printf("  maximum Bfactor value : %3.2f\n", bfactmax_);
    log.printf("  initial Bfactor value : %3.2f\n", bfactini);
  }
}

// prepare auxiliary vectors
void EMMIVOX::get_auxiliary_vectors() {
// number of atoms
  unsigned natoms = Model_res_.size();
// clear lists
  pref_.clear();
  invs2_.clear();
  cut_.clear();
// resize
  pref_.resize(natoms);
  invs2_.resize(natoms);
  cut_.resize(natoms);
// cycle on all atoms
  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for(unsigned im=0; im<natoms; ++im) {
    // get atom type
    unsigned atype = Model_type_[im];
    // get residue/chain IDs
    std::pair<unsigned,std::string> key = Model_res_[im];
    // get bfactor
    double bfact = Model_b_[key];
    // sigma for 5 gaussians
    Vector5d m_s = Model_s_[atype];
    // calculate constant quantities
    Vector5d pref, invs2;
    // calculate cutoff
    double n = 0.0;
    double d = 0.0;
    for(unsigned j=0; j<5; ++j) {
      // total value of b
      double m_b = m_s[j] + bfact/4.0;
      // calculate invs2
      invs2[j] = 1.0/(inv_pi2_*m_b);
      // prefactor
      pref[j]  = cfact_[atype][j] * pow(invs2[j],1.5);
      // cutoff
      n += pref[j] / invs2[j];
      d += pref[j];
    }
    // put into global lists
    pref_[im]  = pref;
    invs2_[im] = invs2;
    cut_[im] = std::min(nl_dist_cutoff_, sqrt(n/d)*nl_gauss_cutoff_);
  }
  // push to GPU
  push_auxiliary_gpu();
}

void EMMIVOX::push_auxiliary_gpu() {
  // 1) create vector of pref_ and invs2_
  int natoms = Model_type_.size();
  std::vector<double> pref(5*natoms), invs2(5*natoms);
  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for(int i=0; i<natoms; ++i) {
    for(int j=0; j<5; ++j) {
      pref[i+j*natoms]  = pref_[i][j];
      invs2[i+j*natoms] = invs2_[i][j];
    }
  }
  // 2) initialize gpu tensors
  pref_gpu_  = torch::from_blob(pref.data(),  {5,natoms}, torch::kFloat64).clone().to(torch::kFloat32).to(device_t_);
  invs2_gpu_ = torch::from_blob(invs2.data(), {5,natoms}, torch::kFloat64).clone().to(torch::kFloat32).to(device_t_);
}

void EMMIVOX::get_close_residues() {
  // clear neighbor list
  nl_res_.clear();
  nl_res_.resize(Model_rlist_.size());

  // loop in parallel
  #pragma omp parallel num_threads(OpenMP::getNumThreads())
  {
    // private variable
    std::vector< std::vector<unsigned> > nl_res_l(Model_rlist_.size());
    // cycle on residues/chains #1
    #pragma omp for
    for(unsigned i=0; i<Model_rlist_.size()-1; ++i) {

      // key1: pair of residue/chain IDs
      std::pair<unsigned,std::string> key1 = Model_rlist_[i];

      // cycle over residues/chains #2
      for(unsigned j=i+1; j<Model_rlist_.size(); ++j) {

        // key2: pair of residue/chain IDs
        std::pair<unsigned,std::string> key2 = Model_rlist_[j];

        // set flag neighbor
        bool neigh = false;

        // cycle over all the atoms belonging to key1
        for(unsigned im1=0; im1<Model_resmap_[key1].size(); ++im1) {
          // get atom position #1
          Vector pos1 = getPosition(Model_resmap_[key1][im1]);
          // cycle over all the atoms belonging to key2
          for(unsigned im2=0; im2<Model_resmap_[key2].size(); ++im2) {
            // get atom position #2
            Vector pos2 = getPosition(Model_resmap_[key2][im2]);
            // if closer than 0.5 nm, then residues key1 and key2 are neighbors
            if(delta(pos1,pos2).modulo()<0.5) {
              // set neighbors
              neigh = true;
              // and exit
              break;
            }
          }
          // check if neighbor already found
          if(neigh) {
            break;
          }
        }

        // if neighbors, add to local list
        if(neigh) {
          nl_res_l[i].push_back(j);
          nl_res_l[j].push_back(i);
        }
      }
    }
    // add to global list
    #pragma omp critical
    {
      for(unsigned i=0; i<nl_res_.size(); ++i) {
        nl_res_[i].insert(nl_res_[i].end(), nl_res_l[i].begin(), nl_res_l[i].end());
      }
    }
  }
}

void EMMIVOX::doMonteCarloBfact() {
// update residue neighbor list
  get_close_residues();

// cycle over residues/chains
  for(unsigned ir=0; ir<Model_rlist_.size(); ++ir) {

    // key: pair of residue/chain IDs
    std::pair<unsigned,std::string> key = Model_rlist_[ir];
    // old bfactor
    double bfactold = Model_b_[key];

    // propose move in bfactor
    double bfactnew = bfactold + dbfact_ * ( 2.0 * random_.RandU01() - 1.0 );
    // check boundaries
    if(bfactnew > bfactmax_) {
      bfactnew = 2.0*bfactmax_ - bfactnew;
    }
    if(bfactnew < bfactmin_) {
      bfactnew = 2.0*bfactmin_ - bfactnew;
    }

    // useful quantities
    std::map<unsigned, double> deltaov;

    #pragma omp parallel num_threads(OpenMP::getNumThreads())
    {
      // private variables
      std::map<unsigned, double> deltaov_l;
      #pragma omp for
      // cycle over all the atoms belonging to key (residue/chain)
      for(unsigned ia=0; ia<Model_resmap_[key].size(); ++ia) {

        // get atom id
        unsigned im = Model_resmap_[key][ia];
        // get atom type
        unsigned atype = Model_type_[im];
        // sigma for 5 Gaussians
        Vector5d m_s = Model_s_[atype];
        // prefactors
        Vector5d cfact = cfact_[atype];
        // and position
        Vector pos = getPosition(im);

        // cycle on all the neighboring voxels affected by a change in Bfactor
        for(unsigned i=0; i<Model_nb_[im].size(); ++i) {
          // voxel id
          unsigned id = Model_nb_[im][i];
          // get contribution to density in id before change
          double dold = get_overlap(Map_m_[id], pos, cfact, m_s, bfactold);
          // get contribution after change
          double dnew = get_overlap(Map_m_[id], pos, cfact, m_s, bfactnew);
          // update delta density
          deltaov_l[id] += dnew-dold;
        }
      }
      // add to global list
      #pragma omp critical
      {
        for(std::map<unsigned,double>::iterator itov=deltaov_l.begin(); itov!=deltaov_l.end(); ++itov) {
          deltaov[itov->first] += itov->second;
        }
      }
    }

    // now calculate new and old score
    double old_ene = 0.0;
    double new_ene = 0.0;

    // cycle on all affected voxels
    for(std::map<unsigned,double>::iterator itov=deltaov.begin(); itov!=deltaov.end(); ++itov) {
      // id of the component
      unsigned id = itov->first;
      // new value
      double ovmdnew = ovmd_[id]+itov->second;
      // deviations
      double devold = scale_ * ovmd_[id] + offset_ - ovdd_[id];
      double devnew = scale_ * ovmdnew   + offset_ - ovdd_[id];
      // inverse of sigma_min
      double ismin = ismin_[id];
      // scores
      if(devold==0.0) {
        old_ene += -kbt_ * std::log( 0.5 * sqrt2_pi_ * ismin );
      } else {
        old_ene += -kbt_ * std::log( 0.5 / devold * erf ( devold * inv_sqrt2_ * ismin ));
      }
      if(devnew==0.0) {
        new_ene += -kbt_ * std::log( 0.5 * sqrt2_pi_ * ismin );
      } else {
        new_ene += -kbt_ * std::log( 0.5 / devnew * erf ( devnew * inv_sqrt2_ * ismin ));
      }
    }

    // list of neighboring residues
    std::vector<unsigned> close = nl_res_[ir];
    // add restraint to keep Bfactors of close residues close
    for(unsigned i=0; i<close.size(); ++i) {
      // residue/chain IDs of neighbor
      std::pair<unsigned,std::string> keyn = Model_rlist_[close[i]];
      // deviations
      double devold = bfactold - Model_b_[keyn];
      double devnew = bfactnew - Model_b_[keyn];
      // inverse of sigma_min
      double ismin = 1.0 / bfactsig_;
      // scores
      if(devold==0.0) {
        old_ene += -kbt_ * std::log( 0.5 * sqrt2_pi_ * ismin );
      } else {
        old_ene += -kbt_ * std::log( 0.5 / devold * erf ( devold * inv_sqrt2_ * ismin ));
      }
      if(devnew==0.0) {
        new_ene += -kbt_ * std::log( 0.5 * sqrt2_pi_ * ismin );
      } else {
        new_ene += -kbt_ * std::log( 0.5 / devnew * erf ( devnew * inv_sqrt2_ * ismin ));
      }
    }

    // increment number of trials
    MCBtrials_ += 1.0;

    // accept or reject
    bool accept = false;
    if(bfactemin_) {
      if(new_ene < old_ene) {
        accept = true;
      }
    } else {
      accept = doAccept(old_ene, new_ene, kbt_);
    }

    // in case of acceptance
    if(accept) {
      // update acceptance rate
      MCBaccept_ += 1.0;
      // update bfactor
      Model_b_[key] = bfactnew;
      // change all the ovmd_ affected
      for(std::map<unsigned,double>::iterator itov=deltaov.begin(); itov!=deltaov.end(); ++itov) {
        ovmd_[itov->first] += itov->second;
      }
    }

  } // end cycle on bfactors

// update auxiliary lists (to update pref_gpu_, invs2_gpu_, and cut_ on CPU/GPU)
  get_auxiliary_vectors();
// update neighbor list (new cut_ + update pref_nl_gpu_ and invs2_nl_gpu_ on GPU)
  update_neighbor_list();
// recalculate fmod (to update derivatives)
  calculate_fmod();
}

// get overlap
double EMMIVOX::get_overlap(const Vector &d_m, const Vector &m_m,
                            const Vector5d &cfact, const Vector5d &m_s, double bfact) {
  // calculate vector difference
  Vector md = delta(m_m, d_m);
  // norm squared
  double md2 = md[0]*md[0]+md[1]*md[1]+md[2]*md[2];
  // cycle on 5 Gaussians
  double ov_tot = 0.0;
  for(unsigned j=0; j<5; ++j) {
    // total value of b
    double m_b = m_s[j]+bfact/4.0;
    // calculate invs2
    double invs2 = 1.0/(inv_pi2_*m_b);
    // final calculation
    ov_tot += cfact[j] * pow(invs2, 1.5) * std::exp(-0.5 * md2 * invs2);
  }
  return ov_tot;
}

void EMMIVOX::update_neighbor_sphere() {
  // number of atoms
  unsigned natoms = Model_type_.size();
  // clear neighbor sphere
  ns_.clear();
  // store reference positions
  refpos_ = getPositions();

  // cycle on voxels - in parallel
  #pragma omp parallel num_threads(OpenMP::getNumThreads())
  {
    // private variables
    std::vector< std::pair<unsigned,unsigned> > ns_l;
    #pragma omp for
    for(unsigned id=0; id<ovdd_.size(); ++id) {
      // grid point
      Vector d_m = Map_m_[id];
      // cycle on atoms
      for(unsigned im=0; im<natoms; ++im) {
        // calculate distance
        double dist = delta(getPosition(im), d_m).modulo();
        // add to local list
        if(dist<=2.0*cut_[im]) {
          ns_l.push_back(std::make_pair(id,im));
        }
      }
    }
    // add to global list
    #pragma omp critical
    ns_.insert(ns_.end(), ns_l.begin(), ns_l.end());
  }
}

bool EMMIVOX::do_neighbor_sphere() {
  std::vector<double> dist(getPositions().size());
  bool update = false;

// calculate displacement
  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for(unsigned im=0; im<dist.size(); ++im) {
    dist[im] = delta(getPosition(im),refpos_[im]).modulo()/cut_[im];
  }

// check if update or not
  double maxdist = *max_element(dist.begin(), dist.end());
  if(maxdist>=1.0) {
    update=true;
  }

// return if update or not
  return update;
}

void EMMIVOX::update_neighbor_list() {
  // number of atoms
  unsigned natoms = Model_type_.size();
  // clear neighbor list
  nl_.clear();

  // cycle on neighbour sphere - in parallel
  #pragma omp parallel num_threads(OpenMP::getNumThreads())
  {
    // private variables
    std::vector< std::pair<unsigned,unsigned> > nl_l;
    #pragma omp for
    for(unsigned long long i=0; i<ns_.size(); ++i) {
      // calculate distance
      double dist = delta(Map_m_[ns_[i].first], getPosition(ns_[i].second)).modulo();
      // add to local neighbour list
      if(dist<=cut_[ns_[i].second]) {
        nl_l.push_back(ns_[i]);
      }
    }
    // add to global list
    #pragma omp critical
    nl_.insert(nl_.end(), nl_l.begin(), nl_l.end());
  }

  // new dimension of neighbor list
  unsigned long long nl_size = nl_.size();
  // now resize derivatives
  ovmd_der_.resize(nl_size);

  // in case of B-factors sampling - at the right step
  if(dbfact_>0 && getStep()%MCBstride_==0) {
    // clear vectors
    Model_nb_.clear();
    Model_nb_.resize(natoms);
    // cycle over the neighbor list to creat a list of voxels per atom
    #pragma omp parallel num_threads(OpenMP::getNumThreads())
    {
      // private variables
      std::vector< std::vector<unsigned> > Model_nb_l(natoms);
      #pragma omp for
      for(unsigned long long i=0; i<nl_size; ++i) {
        Model_nb_l[nl_[i].second].push_back(nl_[i].first);
      }
      // add to global list
      #pragma omp critical
      {
        for(unsigned i=0; i<natoms; ++i) {
          Model_nb_[i].insert(Model_nb_[i].end(), Model_nb_l[i].begin(), Model_nb_l[i].end());
        }
      }
    }
  }

  // transfer data to gpu
  update_gpu();
}

void EMMIVOX::update_gpu() {
  // dimension of neighbor list
  unsigned nl_size = nl_.size();
  // create useful vectors
  std::vector<int> nl_id(nl_size), nl_im(nl_size);
  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for(unsigned long long i=0; i<nl_size; ++i) {
    nl_id[i] = static_cast<int>(nl_[i].first);
    nl_im[i] = static_cast<int>(nl_[i].second);
  }
  // create tensors on device
  nl_id_gpu_ = torch::from_blob(nl_id.data(), {nl_size}, torch::kInt32).clone().to(device_t_);
  nl_im_gpu_ = torch::from_blob(nl_im.data(), {nl_size}, torch::kInt32).clone().to(device_t_);
  // now we need to create pref_nl_gpu_ [5,nl_size]
  pref_nl_gpu_  = torch::index_select(pref_gpu_,1,nl_im_gpu_);
  // and invs2_nl_gpu_ [5,nl_size]
  invs2_nl_gpu_ = torch::index_select(invs2_gpu_,1,nl_im_gpu_);
  // and Map_m_nl_gpu_ [3,nl_size]
  Map_m_nl_gpu_ = torch::index_select(Map_m_gpu_,1,nl_id_gpu_);
}

void EMMIVOX::prepare() {
  if(getExchangeStep()) {
    first_time_=true;
  }
}

// calculate forward model on gpu
void EMMIVOX::calculate_fmod() {
  // number of atoms
  int natoms = Model_type_.size();
  // number of data points
  int nd = ovdd_.size();

  // fill positions in in parallel
  std::vector<double> posg(3*natoms);
  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for (int i=0; i<natoms; ++i) {
    // fill vectors
    posg[i]          = getPosition(i)[0];
    posg[i+natoms]   = getPosition(i)[1];
    posg[i+2*natoms] = getPosition(i)[2];
  }
  // transfer positions to pos_gpu [3,natoms]
  torch::Tensor pos_gpu = torch::from_blob(posg.data(), {3,natoms}, torch::kFloat64).to(torch::kFloat32).to(device_t_);
  // create pos_nl_gpu_ [3,nl_size]
  torch::Tensor pos_nl_gpu = torch::index_select(pos_gpu,1,nl_im_gpu_);
  // calculate vector difference [3,nl_size]
  torch::Tensor md = Map_m_nl_gpu_ - pos_nl_gpu;
  // calculate norm squared by column [1,nl_size]
  torch::Tensor md2 = torch::sum(md*md,0);
  // calculate density [5,nl_size]
  torch::Tensor ov = pref_nl_gpu_ * torch::exp(-0.5 * md2 * invs2_nl_gpu_);
  // and derivatives [5,nl_size]
  ovmd_der_gpu_ = invs2_nl_gpu_ * ov;
  // sum density over 5 columns [1,nl_size]
  ov = torch::sum(ov,0);
  // sum contributions from the same atom
  auto options = torch::TensorOptions().device(device_t_).dtype(torch::kFloat32);
  ovmd_gpu_ = torch::zeros({nd}, options);
  ovmd_gpu_.index_add_(0, nl_id_gpu_, ov);
  // sum derivatives over 5 rows [1,nl_size] and multiply by md [3,nl_size]
  ovmd_der_gpu_ = md * torch::sum(ovmd_der_gpu_,0);

  // in case of metainference: average them across replicas
  if(!no_aver_ && nrep_>1) {
    // communicate ovmd_gpu_ to CPU [1, nd]
    torch::Tensor ovmd_cpu = ovmd_gpu_.detach().to(torch::kCPU).to(torch::kFloat64);
    // and put them in ovmd_
    ovmd_ = std::vector<double>(ovmd_cpu.data_ptr<double>(), ovmd_cpu.data_ptr<double>() + ovmd_cpu.numel());
    // sum across replicas
    multi_sim_comm.Sum(&ovmd_[0], nd);
    // and divide by number of replicas
    double escale = 1.0 / static_cast<double>(nrep_);
    for(int i=0; i<nd; ++i) {
      ovmd_[i] *= escale;
    }
    // put back on device
    ovmd_gpu_ = torch::from_blob(ovmd_.data(), {nd}, torch::kFloat64).to(torch::kFloat32).to(device_t_);
  }

  // communicate back model density
  // this is needed only in certain situations
  long int step = getStep();
  bool do_comm = false;
  if(mapstride_>0 && step%mapstride_==0) {
    do_comm = true;
  }
  if(dbfact_>0    && step%MCBstride_==0) {
    do_comm = true;
  }
  if(do_corr_) {
    do_comm = true;
  }
  // in case of metainference: already communicated
  if(!no_aver_ && nrep_>1) {
    do_comm = false;
  }
  if(do_comm) {
    // communicate ovmd_gpu_ to CPU [1, nd]
    torch::Tensor ovmd_cpu = ovmd_gpu_.detach().to(torch::kCPU).to(torch::kFloat64);
    // and put them in ovmd_
    ovmd_ = std::vector<double>(ovmd_cpu.data_ptr<double>(), ovmd_cpu.data_ptr<double>() + ovmd_cpu.numel());
  }
}

// calculate score
void EMMIVOX::calculate_score() {
  // number of atoms
  int natoms = Model_type_.size();

  // calculate deviation model/data [1, nd]
  torch::Tensor dev = scale_ * ovmd_gpu_ + offset_ - ovdd_gpu_;
  // error function [1, nd]
  torch::Tensor errf = torch::erf( dev * inv_sqrt2_ * ismin_gpu_ );
  // take care of dev = zero
  torch::Tensor zeros_d = torch::ne(dev, 0.0);
  // redefine dev
  dev = dev * zeros_d + eps_ * torch::logical_not(zeros_d);
  // take care of errf = zero
  torch::Tensor zeros_e = torch::ne(errf, 0.0);
  // redefine errf
  errf = errf * zeros_e + eps_ * torch::logical_not(zeros_e);
  // logical AND: both dev and errf different from zero
  torch::Tensor zeros = torch::logical_and(zeros_d, zeros_e);
  // energy - with limit dev going to zero
  torch::Tensor ene = 0.5 * ( errf / dev * zeros + torch::logical_not(zeros) * sqrt2_pi_ *  ismin_gpu_);
  // logarithm and sum
  ene = -kbt_ * torch::sum(torch::log(ene));
  // and derivatives [1, nd]
  torch::Tensor d_der = -kbt_ * zeros * ( sqrt2_pi_ * torch::exp( -0.5 * dev * dev * ismin_gpu_ * ismin_gpu_ ) * ismin_gpu_ / errf - 1.0 / dev );
  // tensor for derivatives wrt atoms [1, nl_size]
  torch::Tensor der_gpu = torch::index_select(d_der,0,nl_id_gpu_);
  // multiply by ovmd_der_gpu_ and scale [3, nl_size]
  der_gpu = ovmd_der_gpu_ * scale_ * der_gpu;
  // sum contributions for each atom
  auto options = torch::TensorOptions().device(device_t_).dtype(torch::kFloat32);
  torch::Tensor atoms_der_gpu = torch::zeros({3,natoms}, options);
  atoms_der_gpu.index_add_(1, nl_im_gpu_, der_gpu);

  // FINAL STUFF
  //
  // 1) communicate total energy to CPU
  torch::Tensor ene_cpu = ene.detach().to(torch::kCPU).to(torch::kFloat64);
  ene_ = *ene_cpu.data_ptr<double>();
  // with marginal, simply multiply by number of replicas!
  if(!no_aver_ && nrep_>1) {
    ene_ *= static_cast<double>(nrep_);
  }
  //
  // 2) communicate derivatives to CPU
  torch::Tensor atom_der_cpu = atoms_der_gpu.detach().to(torch::kCPU).to(torch::kFloat64);
  // convert to std::vector<double>
  std::vector<double> atom_der = std::vector<double>(atom_der_cpu.data_ptr<double>(), atom_der_cpu.data_ptr<double>() + atom_der_cpu.numel());
  // and put in atom_der_
  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for(int i=0; i<natoms; ++i) {
    atom_der_[i] = Vector(atom_der[i],atom_der[i+natoms],atom_der[i+2*natoms]);
  }
  //
  // 3) calculate virial on CPU
  Tensor virial;
  // declare omp reduction for Tensors
  #pragma omp declare reduction( sumTensor : Tensor : omp_out += omp_in )

  #pragma omp parallel for num_threads(OpenMP::getNumThreads()) reduction (sumTensor : virial)
  for(int i=0; i<natoms; ++i) {
    virial += Tensor(getPosition(i), -atom_der_[i]);
  }
  // store virial
  virial_ = virial;
}

void EMMIVOX::calculate() {
  // get time step
  long int step = getStep();

  // set temperature value
  getPntrToComponent("kbt")->set(kbt_);

  // neighbor list update
  if(first_time_ || getExchangeStep() || step%nl_stride_==0) {
    // check if time to update neighbor sphere
    bool update = false;
    if(first_time_ || getExchangeStep()) {
      update = true;
    } else {
      update = do_neighbor_sphere();
    }
    // update neighbor sphere
    if(update) {
      update_neighbor_sphere();
    }
    // update neighbor list
    update_neighbor_list();
    // set flag
    first_time_=false;
  }

  // calculate forward model
  calculate_fmod();

  // Monte Carlo on bfactors
  if(dbfact_>0) {
    double acc = 0.0;
    // do Monte Carlo
    if(step%MCBstride_==0 && !getExchangeStep() && step>0) {
      doMonteCarloBfact();
    }
    // calculate acceptance ratio
    if(MCBtrials_>0) {
      acc = MCBaccept_ / MCBtrials_;
    }
    // set value
    getPntrToComponent("accB")->set(acc);
  }

  // calculate score
  calculate_score();

  // set score, virial, and derivatives
  Value* score = getPntrToComponent("scoreb");
  score->set(ene_);
  setBoxDerivatives(score, virial_);
  #pragma omp parallel for
  for(unsigned i=0; i<atom_der_.size(); ++i) {
    setAtomsDerivatives(score, i, atom_der_[i]);
  }
  // set scale and offset value
  getPntrToComponent("scale")->set(scale_);
  getPntrToComponent("offset")->set(offset_);
  // calculate correlation coefficient
  if(do_corr_) {
    calculate_corr();
  }
  // PRINT other quantities to files
  // - status file
  if(step%statusstride_==0) {
    print_status(step);
  }
  // - density file
  if(mapstride_>0 && step%mapstride_==0) {
    write_model_density(step);
  }
}

void EMMIVOX::calculate_corr() {
// number of data points
  double nd = static_cast<double>(ovdd_.size());
// average ovmd_ and ovdd_
  double ave_md = std::accumulate(ovmd_.begin(), ovmd_.end(), 0.) / nd;
  double ave_dd = std::accumulate(ovdd_.begin(), ovdd_.end(), 0.) / nd;
// calculate correlation
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
// correlation coefficient
  double cc = num / sqrt(den1*den2);
// set plumed
  getPntrToComponent("corr")->set(cc);
}


}
}

#endif
