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
#include "core/SetupMolInfo.h"
#include "core/ActionSet.h"
#include "tools/File.h"
#include "tools/OpenMP.h"

#include <string>
#include <cmath>
#include <map>
#include <numeric>
#include <ctime>
#include "tools/Random.h"

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
  vector<double> sigma_;
  vector<double> sigma_min_;
  vector<double> sigma_max_;
  vector<double> dsigma_;
// neighbor list
  double   nl_cutoff_;
  unsigned nl_stride_;
  double   ns_cutoff_;
  bool first_time_;
  vector<unsigned> nl_;
  vector<unsigned> ns_;
  vector<Vector> refpos_;
// averaging
  bool no_aver_;
// Monte Carlo stuff
  int      MCstride_;
  double   MCaccept_;
  double   MCtrials_;
  Random   random_;
// Bfact Monte Carlo
  int      MCBstride_;
  double   MCBaccept_;
  double   MCBtrials_;
  double   dbfact_;
  double   bfactmin_;
  double   bfactmax_;
  bool     readbf_;
  // status stuff
  unsigned int statusstride_;
  string       statusfilename_;
  OFile        statusfile_;
  bool         first_status_;
  // Scale Monte Carlo
  double scale_;
  double scale_min_;
  double scale_max_;
  double dscale_;
  double offset_;
  double doffset_;
  double MCSaccept_;
  double MCStrials_;
  // prior exponent
  double prior_;
  // noise type
  unsigned noise_;
  // total energy
  double ene_;
  // model overlap file
  unsigned int ovstride_;
  string       ovfilename_;

// write file with model overlap
  void write_model_overlap(long int step);
// get median of vector
  double get_median(vector<double> &v);
// read and write status
  void read_status();
  void print_status(long int step);
// accept or reject
  bool doAccept(double oldE, double newE, double kbt);
// do MonteCarlo
  void doMonteCarlo(vector<double> &eneg);
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
  void calculate_useful_stuff(double reso);
  void get_auxiliary_vectors();
// calculate overlap between two Gaussians
  double get_overlap_der(const Vector &d_m, const Vector &m_m,
                         const Vector5d &pref, const Vector5d &invs2,
                         Vector &ov_der);
  double get_overlap(const Vector &d_m, const Vector &m_m, double d_s,
                     const Vector5d &cfact, const Vector5d &m_s, double bfact);
// update the neighbor list
  void update_neighbor_list();
// update the neighbor sphere
  void update_neighbor_sphere();
  bool do_neighbor_sphere();
// calculate overlap
  void calculate_overlap();
// Gaussian noise
  double calculate_Gauss_group(unsigned igroup, double sigma,
                               double scale, double offset, vector<double> &GMMid_der);
  double calculate_Gauss(vector<double> &eneg, vector<double> &sigma,
                         double scale, double offset, vector<double> &GMMid_der);
// Outliers noise
  double calculate_Outliers_group(unsigned igroup, double sigma,
                                  double scale, double offset, vector<double> &GMMid_der);
  double calculate_Outliers(vector<double> &eneg, vector<double> &sigma,
                            double scale, double offset, vector<double> &GMMid_der);
// Marginal noise
  double calculate_Marginal(double scale, double offset, vector<double> &GMMid_der);

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
  keys.add("compulsory","NOISETYPE","functional form of the noise (GAUSS, OUTLIERS, MARGINAL)");
  keys.add("compulsory","NORM_DENSITY","integral of the experimental density");
  keys.add("compulsory","WRITE_STRIDE","write the status to a file every N steps, this can be used for restart");
  keys.add("optional","SIGMA0","initial value of the uncertainty");
  keys.add("optional","DSIGMA","MC step for uncertainties");
  keys.add("optional","MC_STRIDE", "Monte Carlo stride");
  keys.add("optional","DBFACT","MC step for bfactor");
  keys.add("optional","BFACT_MAX","Maximum value of bfactor");
  keys.add("optional","MCBFACT_STRIDE", "Bfactor Monte Carlo stride");
  keys.addFlag("READ_BFACT",false,"read Bfactor from status file at restart");
  keys.add("optional","ERR_FILE","file with experimental errors");
  keys.add("optional","STATUS_FILE","write a file with all the data useful for restart");
  keys.add("optional","SCALE_MIN","minimum scale");
  keys.add("optional","SCALE_MAX","maximum scale");
  keys.add("optional","DSCALE","maximum scale MC move");
  keys.add("optional","SCALE","scale factor");
  keys.add("optional","OFFSET","offset");
  keys.add("optional","DOFFSET","maximum offset MC move");
  keys.add("optional","TEMP","temperature");
  keys.add("optional","PRIOR", "exponent of uncertainty prior");
  keys.add("optional","WRITE_OV_STRIDE","write model overlaps every N steps");
  keys.add("optional","WRITE_OV","write a file with model overlaps");
  keys.addFlag("NO_AVER",false,"don't do ensemble averaging in multi-replica mode");
  componentsAreNotOptional(keys);
  keys.addOutputComponent("scoreb","default","Bayesian score");
  keys.addOutputComponent("acc",   "NOISETYPE","MC acceptance for uncertainty");
  keys.addOutputComponent("accB",  "default", "Bfactor MC acceptance");
  keys.addOutputComponent("scale", "default","scale factor");
  keys.addOutputComponent("offset","default","offset");
  keys.addOutputComponent("accscale", "default","MC acceptance for scale");
}

EMMIVOX::EMMIVOX(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  inv_sqrt2_(0.707106781186548),
  sqrt2_pi_(0.797884560802865),
  inv_pi2_(0.050660591821169),
  first_time_(true), no_aver_(false),
  MCstride_(1), MCaccept_(0.), MCtrials_(0.),
  MCBstride_(1), MCBaccept_(0.), MCBtrials_(0.),
  dbfact_(0.0), bfactmax_(4.0), readbf_(false),
  statusstride_(0), first_status_(true),
  scale_(1.), dscale_(0.), offset_(0.), doffset_(0.),
  MCSaccept_(0.), MCStrials_(0.),
  prior_(1.), ovstride_(0)
{

  // list of atoms
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS", atoms);

  // file with experimental data
  string datafile;
  parse("DATA_FILE", datafile);

  // type of data noise
  string noise;
  parse("NOISETYPE",noise);
  if      (noise=="GAUSS")   noise_ = 0;
  else if(noise=="OUTLIERS") noise_ = 1;
  else if(noise=="MARGINAL") noise_ = 2;
  else error("Unknown noise type!");

  // minimum value for error
  double sigma_min;
  parse("SIGMA_MIN", sigma_min);
  if(sigma_min<0) error("SIGMA_MIN should be greater or equal to zero");

  // Monte Carlo in B-factors
  parse("DBFACT", dbfact_);
  parse("BFACT_MAX", bfactmax_);
  parse("MCBFACT_STRIDE", MCBstride_);
  parseFlag("READ_BFACT", readbf_);
  if(dbfact_<0) error("DBFACT should be greater or equal to zero");
  if(dbfact_>0 && MCBstride_<=0) error("you must specify a positive MCBFACT_STRIDE");
  if(dbfact_>0 && bfactmax_<=0)  error("you must specify a positive BFACT_MAX");

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
  }

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

  // exponent of uncertainty prior
  parse("PRIOR",prior_);

  // scale MC
  parse("SCALE", scale_);
  parse("DSCALE",dscale_);
  parse("OFFSET",offset_);
  parse("DOFFSET",doffset_);
  // other parameters
  if(dscale_>0.) {
    parse("SCALE_MIN",scale_min_);
    parse("SCALE_MAX",scale_max_);
    // checks
    if(scale_min_<=0.0) error("SCALE_MIN must be strictly positive");
    if(scale_max_<=scale_min_) error("SCALE_MAX must be greater than SCALE_MIN");
  }

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

  // write overlap file
  parse("WRITE_OV_STRIDE", ovstride_);
  parse("WRITE_OV", ovfilename_);
  if(ovstride_>0 && ovfilename_=="") error("With WRITE_OV_STRIDE you must specify WRITE_OV");

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
  log.printf("  type of data noise : %s\n", noise.c_str());
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
  }
  if(noise_!=2) {
    log.printf("  initial value of the uncertainty : %f\n",sigma_ini);
    log.printf("  max MC move in uncertainty : %f\n",dsigma);
    log.printf("  MC stride : %u\n", MCstride_);
  }
  if(dbfact_>0) {
    log.printf("  max MC move in bfactor : %f\n",dbfact_);
    log.printf("  Bfactor MC stride : %u\n", MCBstride_);
  }
  if(errfile.size()>0) log.printf("  reading experimental errors from file : %s\n", errfile.c_str());
  log.printf("  temperature of the system in energy unit : %f\n",kbt_);
  log.printf("  prior exponent : %f\n",prior_);
  log.printf("  number of replicas for averaging: %u\n",nrep_);
  log.printf("  id of the replica : %u\n",replica_);
  if(ovstride_>0) {
    log.printf("  stride for writing model overlaps : %u\n",ovstride_);
    log.printf("  file for writing model overlaps : %s\n", ovfilename_.c_str());
  }

  // calculate model GMM constant parameters
  vector<double> GMM_m_w = get_GMM_m(atoms);

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
  if(dbfact_>0) {addComponent("accB");   componentIsNotPeriodic("accB");}
  if(noise_!=2) {if(dsigma_[0]>0) {addComponent("acc");    componentIsNotPeriodic("acc");}}
  if(dscale_>0.0) {addComponent("accscale");  componentIsNotPeriodic("accscale");}

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
  for(int i=0; i<ovmd_.size(); ++i) {
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
  double bf;
// open file
  IFile *ifile = new IFile();
  ifile->link(*this);
  if(ifile->FileExist(statusfilename_)) {
    ifile->open(statusfilename_);
    while(ifile->scanField("MD_time", MDtime)) {
      // read sigma only if not marginal noise
      if(noise_!=2) {
        for(unsigned i=0; i<sigma_.size(); ++i) {
          // convert i to string
          std::string num; Tools::convert(i,num);
          // read entries
          ifile->scanField("s"+num, sigma_[i]);
        }
      }
      // read scale and offset
      ifile->scanField("scale", scale_);
      ifile->scanField("offset", offset_);
      // read bfactors
      for(unsigned i=0; i<GMM_m_res_.size(); ++i) {
        // convert i to string
        std::string num; Tools::convert(i,num);
        // read entries
        if(readbf_) ifile->scanField("bfact"+num, GMM_m_b_[GMM_m_res_[i]]);
        else        ifile->scanField("bfact"+num, bf);
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
  // write sigma only if not marginal noise
  if(noise_!=2) {
    for(unsigned i=0; i<sigma_.size(); ++i) {
      // convert i to string
      std::string num; Tools::convert(i,num);
      // print entry
      statusfile_.printField("s"+num, sigma_[i]);
    }
  }
  // always write scale and offset
  statusfile_.printField("scale", scale_);
  statusfile_.printField("offset", offset_);
  // always write bfactors
  for(unsigned i=0; i<GMM_m_res_.size(); ++i) {
    // convert i to string
    std::string num; Tools::convert(i,num);
    // print entry
    statusfile_.printField("bfact"+num, GMM_m_b_[GMM_m_res_[i]]);
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

vector<double> EMMIVOX::get_GMM_m(vector<AtomNumber> &atoms)
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
        Vector5d w = GMM_m_w_[type_map[type_s]];
        GMM_m_w.push_back(w[0]+w[1]+w[2]+w[3]+w[4]);
        // get residue id
        unsigned ires = moldat[0]->getResidueNumber(atoms[i]);
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
    if(GMM_d_beta_[i]>=GMM_d_grps_.size()) error("Check Beta values");
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
}


void EMMIVOX::doMonteCarlo(vector<double> &eneg)
{
  // local acceptance
  double MCaccept = 0.0;
  // new derivatives per GMM component
  vector<double> new_GMMid_der(GMMid_der_.size(), 0.0);

  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(MCaccept)
  {
    #pragma omp for reduction( + : MCaccept)
    // cycle over all the groups in parallel
    for(unsigned i=0; i<sigma_.size(); ++i) {

      // generate random move
      double new_s = sigma_[i] + dsigma_[i] * ( 2.0 * random_.RandU01() - 1.0 );
      // check boundaries
      if(new_s > sigma_max_[i]) {new_s = 2.0 * sigma_max_[i] - new_s;}
      if(new_s < sigma_min_[i]) {new_s = 2.0 * sigma_min_[i] - new_s;}

      // calculate new group energy
      double new_ene = 0.0;

      // in case of Gaussian noise
      if(noise_==0) new_ene = calculate_Gauss_group(i, new_s, scale_, offset_, new_GMMid_der);

      // in case of Outliers noise
      if(noise_==1) new_ene = calculate_Outliers_group(i, new_s, scale_, offset_, new_GMMid_der);

      // accept or reject
      bool accept = doAccept(eneg[i], new_ene, kbt_);

      // in case of acceptance
      if(accept) {
        sigma_[i] = new_s;
        eneg[i]   = new_ene;
        MCaccept += 1.0;
        // derivatives of the GMM components in group i
        for(unsigned j=0; j<GMM_d_grps_[i].size(); ++j) GMMid_der_[GMM_d_grps_[i][j]] = new_GMMid_der[GMM_d_grps_[i][j]];
      }
    } // end of cycle over sigmas
  }
  // recalculate total energy from updated energy groups
  ene_ = accumulate(eneg.begin(), eneg.end(), 0.0);
  // update trials
  MCtrials_ += static_cast<double>(sigma_.size());
  // increase acceptance
  MCaccept_ += MCaccept;
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
        Vector pos = getPosition(im);

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
            double dist = delta(pos,getPosition(in)).modulo();
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

    if(noise_==0) {
      for(itov=deltaov.begin(); itov!=deltaov.end(); ++itov) {
        // id of the component
        unsigned id = itov->first;
        // new value
        double ovmdnew = ovmd_[id]+itov->second;
        // scores
        double devold = ( scale_*ovmd_[id]+offset_-ovdd_[id] ) / sigma_[GMM_d_beta_[id]];
        double devnew = ( scale_*ovmdnew+offset_  -ovdd_[id] ) / sigma_[GMM_d_beta_[id]];
        old_ene += 0.5 * kbt_ * devold * devold;
        new_ene += 0.5 * kbt_ * devnew * devnew;
      }
    }
    if(noise_==1) {
      for(itov=deltaov.begin(); itov!=deltaov.end(); ++itov) {
        // id of the component
        unsigned id = itov->first;
        // new value
        double ovmdnew = ovmd_[id]+itov->second;
        // scores
        double devold = ( scale_*ovmd_[id]+offset_-ovdd_[id] ) / sigma_[GMM_d_beta_[id]];
        double devnew = ( scale_*ovmdnew+offset_  -ovdd_[id] ) / sigma_[GMM_d_beta_[id]];
        old_ene += kbt_ * std::log( 1.0 + 0.5 * devold * devold );
        new_ene += kbt_ * std::log( 1.0 + 0.5 * devnew * devnew );
      }
    }
    if(noise_==2) {
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

  // new total energy
  double new_ene;
  // prepare new energy per group - not used anymore
  vector<double> new_eneg(GMM_d_grps_.size(), 0.0);
  // and derivatives
  vector<double> new_GMMid_der(GMMid_der_.size(), 0.0);

  // in case of Gaussian noise
  if(noise_==0) new_ene = calculate_Gauss(new_eneg, sigma_, new_scale, new_off, new_GMMid_der);

  // in case of Outliers noise
  if(noise_==1) new_ene = calculate_Outliers(new_eneg, sigma_, new_scale, new_off, new_GMMid_der);

  // in case of Marginal noise
  if(noise_==2) new_ene = calculate_Marginal(new_scale, new_off, new_GMMid_der);

  // in case sum new energy across replicas
  if(!no_aver_ && nrep_>1) multi_sim_comm.Sum(&new_ene, 1);

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
    scale_ = new_scale;
    offset_ = new_off;
    ene_ = new_ene;
    MCSaccept_ += 1.0;
    GMMid_der_ = new_GMMid_der;
    // in case sum derivatives across replicas
    if(!no_aver_ && nrep_>1) multi_sim_comm.Sum(&GMMid_der_[0], GMMid_der_.size());
  }
}

// get overlap and derivatives
double EMMIVOX::get_overlap_der(const Vector &d_m, const Vector &m_m,
                                const Vector5d &pref, const Vector5d &invs2,
                                Vector &ov_der)
{
  // initialize stuff
  double ov_tot = 0.0;
  // clear derivatives
  ov_der = Vector(0,0,0);
  // calculate vector difference with/without pbc
  Vector md = delta(m_m, d_m);
  // cycle on 5 Gaussians
  for(unsigned j=0; j<5; ++j) {
    // calculate exponent
    double ov = (md[0]*md[0]+md[1]*md[1]+md[2]*md[2])*invs2[j];
    // final calculation
    ov = pref[j] * exp(-0.5*ov);
    // derivatives
    ov_der += ov * Vector(md[0]*invs2[j],md[1]*invs2[j],md[2]*invs2[j]);
    // increase total overlap
    ov_tot += ov;
  }
  return ov_tot;
}

// get overlap
double EMMIVOX::get_overlap(const Vector &d_m, const Vector &m_m, double d_s,
                            const Vector5d &cfact, const Vector5d &m_s, double bfact)
{
  // calculate vector difference with/without pbc
  Vector md = delta(m_m, d_m);
  // cycle on 5 Gaussians
  double ov_tot = 0.0;
  for(unsigned j=0; j<5; ++j) {
    // total value of b
    double m_b = m_s[j]+bfact/4.0;
    // calculate invs2
    double invs2 = 1.0/(d_s+inv_pi2_*m_b);
    // calculate exponent
    double ov = (md[0]*md[0]+md[1]*md[1]+md[2]*md[2])*invs2;
    // final calculation
    ov_tot += cfact[j] * pow(invs2, 1.5) * exp(-0.5*ov);
  }
  return ov_tot;
}

void EMMIVOX::update_neighbor_sphere()
{
  // dimension of atom vectors
  unsigned GMM_m_size = GMM_m_type_.size();

  // clear global list and reference positions
  ns_.clear(); refpos_.clear();
  // allocate reference positions
  refpos_.resize(GMM_m_size);

  // store reference positions
  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for(unsigned im=0; im<GMM_m_size; ++im) refpos_[im] = getPosition(im);

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
  for(unsigned im=0; im<refpos_.size(); ++im) dist[im] = delta(getPosition(im),refpos_[im]).modulo();

// check if update or not
  double maxdist = *max_element(dist.begin(), dist.end());
  if(maxdist>=(ns_cutoff_-nl_cutoff_)) update=true;

// return if update or not
  return update;
}

void EMMIVOX::update_neighbor_list()
{
  // dimension of atom vectors
  unsigned GMM_m_size = GMM_m_type_.size();

  // clear global list
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
      double dist = delta(GMM_d_m_[id], getPosition(im)).modulo();
      // add to local neighbour list
      if(dist<=nl_cutoff_) nl_l.push_back(ns_[i]);
    }
    // add to global list
    #pragma omp critical
    nl_.insert(nl_.end(), nl_l.begin(), nl_l.end());
  }
  // get size
  unsigned tot_size = nl_.size();
  // now resize derivatives
  ovmd_der_.resize(tot_size);
  // in case of B-factors sampling
  if(dbfact_>0) {
    // now cycle over the neighbor list to creat a list of voxels per atom
    GMM_m_nb_.clear(); GMM_m_nb_.resize(GMM_m_size);
    GMM_d_nb_.clear(); GMM_d_nb_.resize(ovdd_.size());
    unsigned id, im;
    for(unsigned i=0; i<tot_size; ++i) {
      id = nl_[i] / GMM_m_size;
      im = nl_[i] % GMM_m_size;
      GMM_m_nb_[im].push_back(id);
      GMM_d_nb_[id].push_back(im);
    }
  }
}

void EMMIVOX::prepare()
{
  if(getExchangeStep()) first_time_=true;
}

// overlap calculator
void EMMIVOX::calculate_overlap() {

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

  // clear overlap vector
  for(unsigned i=0; i<ovmd_.size(); ++i) ovmd_[i] = 0.0;

  // we have to cycle over all model and data GMM components in the neighbor list
  unsigned GMM_m_size = GMM_m_type_.size();
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
      v_[id] += get_overlap_der(GMM_d_m_[id],getPosition(im),pref_[im],invs2_[im],ovmd_der_[i]);
    }
    #pragma omp critical
    for(unsigned i=0; i<ovmd_.size(); ++i) ovmd_[i] += v_[i];
  }

}

void EMMIVOX::calculate()
{

// calculate CV
  calculate_overlap();

  // rescale factor for ensemble average
  double escale = 1.0 / static_cast<double>(nrep_);

  // in case of ensemble averaging, calculate average overlap
  if(!no_aver_ && nrep_>1) {
    multi_sim_comm.Sum(&ovmd_[0], ovmd_.size());
    for(unsigned i=0; i<ovmd_.size(); ++i) ovmd_[i] *= escale;
  }

  // get time step
  long int step = getStep();

  // Monte Carlo on b factors - this will change the overlaps
  if(dbfact_>0) {
    // do Monte Carlo
    if(step%MCBstride_==0 && !getExchangeStep()) doMonteCarloBfact();
    // calculate acceptance ratio
    double acc = MCBaccept_ / MCBtrials_;
    // set value
    getPntrToComponent("accB")->set(acc);
  }

  // write model overlap to file
  if(ovstride_>0 && step%ovstride_==0) write_model_overlap(step);

  // prepare energy per group
  vector<double> eneg(GMM_d_grps_.size(), 0.0);

  // calculate total energy and energies per group
  // Gaussian noise
  if(noise_==0) ene_ = calculate_Gauss(eneg, sigma_, scale_, offset_, GMMid_der_);

  // Outliers noise
  if(noise_==1) ene_ = calculate_Outliers(eneg, sigma_, scale_, offset_, GMMid_der_);

  // Marginal noise
  if(noise_==2) ene_ = calculate_Marginal(scale_, offset_, GMMid_der_);

  // Move SIGMAs
  // This part is needed only for Gaussian and Outliers noise models
  if(noise_!=2) {
    // do Monte Carlo
    if(dsigma_[0]>0 && step%MCstride_==0 && !getExchangeStep()) doMonteCarlo(eneg);
    if(dsigma_[0]>0) {
      // calculate acceptance ratio
      double acc = MCaccept_ / MCtrials_;
      // set value
      getPntrToComponent("acc")->set(acc);
    }
  }

  // up to now, energy, energy per group and derivatives are per replica
  // in case of ensemble averaging: sum across replicas
  if(!no_aver_ && nrep_>1) {
    multi_sim_comm.Sum(&GMMid_der_[0], GMMid_der_.size());
    multi_sim_comm.Sum(&ene_, 1);
  }

  // Monte Carlo on scale
  if(dscale_>0) {
    // do Monte Carlo
    if(step%MCstride_==0 && !getExchangeStep()) doMonteCarloScale();
    // calculate acceptance ratio
    double acc = MCSaccept_ / MCStrials_;
    // set acceptance value
    getPntrToComponent("accscale")->set(acc);
  }
  // set scale and offset value
  getPntrToComponent("scale")->set(scale_);
  getPntrToComponent("offset")->set(offset_);

  // print status
  if(step%statusstride_==0) print_status(step);

  // clear atom derivatives
  for(unsigned i=0; i<atom_der_.size(); ++i) atom_der_[i] = Vector(0,0,0);

  // calculate atom derivatives and virial
  Tensor virial;
  // declare omp reduction for Tensors
  #pragma omp declare reduction( sumTensor : Tensor : omp_out += omp_in )

  #pragma omp parallel num_threads(OpenMP::getNumThreads())
  {
    // private stuff
    vector<Vector> atom_der(atom_der_.size(), Vector(0,0,0));
    unsigned id, im;
    Vector tot_der;
    #pragma omp for reduction (sumTensor : virial)
    for(unsigned i=0; i<nl_.size(); ++i) {
      // get indexes of data and model component
      id = nl_[i] / GMM_m_type_.size();
      im = nl_[i] % GMM_m_type_.size();
      // chain rule + replica normalization
      tot_der = GMMid_der_[id] * ovmd_der_[i] * escale * scale_;
      // increment atom derivatives
      atom_der[im] += tot_der;
      // and virial
      virial += Tensor(getPosition(im), -tot_der);
    }
    #pragma omp critical
    for(unsigned i=0; i<atom_der_.size(); ++i) atom_der_[i] += atom_der[i];
  }

  #pragma omp parallel for
  for(unsigned i=0; i<atom_der_.size(); ++i) setAtomsDerivatives(getPntrToComponent("scoreb"), i, atom_der_[i]);

  // set score and virial
  getPntrToComponent("scoreb")->set(ene_);
  setBoxDerivatives(getPntrToComponent("scoreb"), virial);
}

double EMMIVOX::calculate_Gauss_group(unsigned igroup, double sigma,
                                      double scale, double offset, vector<double> &GMMid_der)
{
  double eneg = 0.0;
  // cycle on all the members of the group
  for(unsigned j=0; j<GMM_d_grps_[igroup].size(); ++j) {
    // id of the GMM component
    int id = GMM_d_grps_[igroup][j];
    // calculate deviation
    double dev = ( scale * ovmd_[id] + offset - ovdd_[id] ) / sigma;
    // add to group energy
    eneg += dev * dev;
    // store derivative for later
    GMMid_der[id] = kbt_ * dev / sigma;
  }
  // add normalizations and prior
  eneg = kbt_ * ( 0.5 * eneg + (static_cast<double>(GMM_d_grps_[igroup].size())+prior_) * std::log(sigma) );
  // return energy per group
  return eneg;
}

double EMMIVOX::calculate_Gauss(vector<double> &eneg, vector<double> &sigma,
                                double scale, double offset, vector<double> &GMMid_der)
{
  double ene = 0.0;
  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(ene)
  {
    // cycle on all the GMM groups in parallel
    #pragma omp for reduction( + : ene)
    for(unsigned i=0; i<GMM_d_grps_.size(); ++i) {
      // get group contribution
      eneg[i] = calculate_Gauss_group(i, sigma[i], scale, offset, GMMid_der);
      // add to total energy
      ene += eneg[i];
    }
  }
  // return total energy
  return ene;
}

double EMMIVOX::calculate_Outliers_group(unsigned igroup,
    double sigma, double scale, double offset, vector<double> &GMMid_der)
{
  double eneg = 0.0;
  // cycle on all the members of the group
  for(unsigned j=0; j<GMM_d_grps_[igroup].size(); ++j) {
    // id of the GMM component
    int id = GMM_d_grps_[igroup][j];
    // calculate deviation
    double dev = ( scale * ovmd_[id] + offset - ovdd_[id] ) / sigma;
    // add to group energy
    eneg += std::log( 1.0 + 0.5 * dev * dev );
    // store derivative for later
    GMMid_der[id] = kbt_ / ( 1.0 + 0.5 * dev * dev ) * dev / sigma;
  }
  // add normalizations and prior
  eneg = kbt_ * ( eneg + (static_cast<double>(GMM_d_grps_[igroup].size())+prior_) * std::log(sigma) );
  // return energy per group
  return eneg;
}

double EMMIVOX::calculate_Outliers(vector<double> &eneg, vector<double> &sigma,
                                   double scale, double offset, vector<double> &GMMid_der)
{
  double ene = 0.0;
  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(ene)
  {
    // cycle on all the GMM groups in parallel
    #pragma omp for reduction( + : ene)
    for(unsigned i=0; i<GMM_d_grps_.size(); ++i) {
      // get group contribution
      eneg[i] = calculate_Outliers_group(i, sigma[i], scale, offset, GMMid_der);
      // add to total energy
      ene += eneg[i];
    }
  }
  // return total energy
  return ene;
}

double EMMIVOX::calculate_Marginal(double scale, double offset, vector<double> &GMMid_der)
{
  double ene = 0.0;
  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(ene)
  {
    // cycle on all the GMM groups
    #pragma omp for reduction( + : ene)
    for(unsigned i=0; i<GMM_d_grps_.size(); ++i) {
      // cycle on all the members of the group
      for(unsigned j=0; j<GMM_d_grps_[i].size(); ++j) {
        // id of the GMM component
        int id = GMM_d_grps_[i][j];
        // calculate deviation
        double dev = ( scale * ovmd_[id] + offset - ovdd_[id] );
        // calculate errf
        double errf = erf ( dev * inv_sqrt2_ / sigma_min_[i] );
        // add to  energy
        ene += -kbt_ * std::log ( 0.5 / dev * errf ) ;
        // store derivative for later
        GMMid_der[id] = - kbt_/errf*sqrt2_pi_*exp(-0.5*dev*dev/sigma_min_[i]/sigma_min_[i])/sigma_min_[i]+kbt_/dev;
      }
    }
  }
  // return total energy
  return ene;
}

}
}
