/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2020 of Glen Hocky

The FISST module is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The FISST module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "bias/Bias.h"
#include "core/ActionRegister.h"
#include "core/Atoms.h"
#include "core/PlumedMain.h"
#include "tools/File.h"
#include "tools/Matrix.h"
#include "tools/Random.h"
#include "legendre_rule_fast.h"

#include <iostream>


using namespace PLMD;
using namespace bias;

//namespace is lowercase to match
//module names being all lowercase

namespace PLMD {
namespace fisst {

//+PLUMEDOC FISSTMOD_BIAS FISST
/*
Compute and apply the optimal linear force on an observable to enhance sampling of conformational distributions over a range of applied forces.

This method is described in \cite Hartmann-FISST-2019

If the system's Hamiltonian is given by:
\f[
    H(\vec{p},\vec{q}) = \sum_{j} \frac{p_j^2}{2m_j} + U(\vec{q}),
\f]

This bias modifies the Hamiltonian to be:
\f[
  H'(\vec{p},\vec{q}) = H(\vec{p},\vec{q}) - \bar{F} Q
\f]

where for CV \f$Q\f$, a coupling constant \f${\bar{F}}\f$ is determined
adaptively according to the FISST algorithm.

Specifically,
\f[
\bar{F}(Q)=\frac{ \int_{F_{min}}^{F_{max}} e^{\beta F Q(\vec{q})} \omega(F) F dF}{\int_{F_{min}}^{F_{max}} e^{\beta F Q(\vec{q})} \omega(F) dF},
\f]

where \f$\vec{q}\f$ are the molecular coordinates of the system, and \f$w(F)\f$ is a weighting function that is learned on the fly for each force by the FISST algorithm (starting from an initial weight distribution, uniform by default).

The target for \f$w(F)=1/Z_q(F)\f$, where
\f[
    Z_q(F) \equiv \int d\vec{q} e^{-\beta U(\vec{q}) + \beta F Q(\vec{q})}.
\f]

FISST also computes and writes Observable Weights \f$W_F(\vec{q}_t)\f$ for a molecular configuration at time \f$t\f$, so that averages of other quantities \f$A(\vec{q})\f$ can be reconstructed later at different force values (over a trajectory with \f$T\f$ samples):
\f[
    \langle A \rangle_F = \frac{1}{T} \sum_t W_F(\vec{q}_t) A(\vec{q}_t).
\f]


\par Examples

In the following example, an adaptive restraint is learned to bias the distance between two atoms in a system, for a force range of 0-100 pN.

\plumedfile
UNITS LENGTH=A TIME=fs ENERGY=kcal/mol

b1: GROUP ATOMS=1
b2: GROUP ATOMS=12

dend: DISTANCE ATOMS=b1,b2

#The conversion factor is 69.4786 pN = 1 kcal/mol/Angstrom

#0 pN to 100 pN
f: FISST MIN_FORCE=0 MAX_FORCE=1.44 PERIOD=100 NINTERPOLATE=31 ARG=dend OUT_RESTART=pull.restart.txt OUT_OBSERVABLE=pull.observable.txt OBSERVABLE_FREQ=1000

PRINT ARG=dend,f.dend_fbar,f.bias,f.force2 FILE=pull.colvar.txt STRIDE=1000
\endplumedfile


*/
//+ENDPLUMEDOC


class FISST : public Bias {


private:
  /*We will get this and store it once, since on-the-fly changing number of CVs will be fatal*/
  const unsigned int ncvs_;
  std::vector<double> center_;
  std::vector<double> current_avg_force_;

  std::vector<double> forces_;
  std::vector<double> force_weight_;
  std::vector<double> gauss_weight_;
  std::vector<double> partition_estimate_;
  std::vector<double> observable_weight_;

  std::string in_restart_name_;
  std::string out_restart_name_;
  std::string out_observable_name_;
  std::string fmt_;
  std::string initial_weight_dist_;
  OFile out_restart_;
  OFile out_observable_;
  IFile in_restart_;
  bool b_freeze_;
  bool b_adaptive_;
  bool b_restart_;
  bool b_write_restart_;
  bool b_write_observable_;
  bool b_first_restart_sample_;
  int period_;
  int reset_period_;
  int observable_freq_;
  int n_interpolation_;
  int n_samples_;
  double kbt_;
  double beta_;
  //change min_force and max_force to vectors if going to do more than one cv
  double max_force_;
  double min_force_;
  double initial_weight_rate_;
  double threshold_;
  Random rand_;


  Value* value_force2_;
  void readInRestart();
  void NormalizeForceWeights();
  /*setup output restart*/
  void setupOutRestart();
  void setupOutObservable();
  /*write output restart*/
  void writeOutRestart();
  void writeOutObservable();
  void update_statistics();
  void update_bias();
  void apply_bias();
  void compute_observable_weight();

public:
  explicit FISST(const ActionOptions&);
  void calculate();
  void update();
  void turnOnDerivatives();
  static void registerKeywords(Keywords& keys);
  ~FISST();
};

PLUMED_REGISTER_ACTION(FISST,"FISST")

void FISST::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","PERIOD","Steps corresponding to the learning rate");
  keys.add("optional","RESET_PERIOD","Reset the learning statistics every time this number of steps comes around.");
  keys.add("compulsory","NINTERPOLATE","Number of grid points on which to do interpolation.");
  keys.add("compulsory","MIN_FORCE","Minimum force (per CV) to use for sampling. Units: [Energy]/[CV]  (can be negative).");
  keys.add("compulsory","MAX_FORCE","Maximum force (per CV) to use for sampling.");
  keys.add("compulsory","CENTER","0","The CV value at which the applied bias energy will be zero");
  keys.add("optional","KBT","The system temperature in units of KB*T. If not provided will be taken from MD code (if available)");

  keys.add("optional","INITIAL_WEIGHT_DIST","Starting distribution for the force weights (options: UNIFORM, EXP, GAUSS).");
  keys.add("optional","INITIAL_WEIGHT_RATE","Rate of decay for exponential and gaussian distributions. W(F)~exp(-r |F|^d).");

  keys.add("optional","RESTART_FMT","the format that should be used to output real numbers in FISST restarts.");
  keys.add("optional","OUT_RESTART","Output file for all information needed to continue FISST simulation."
           "If you have the RESTART directive set (global or for FISST), this file will be appended to."
           "Note that the header will be printed again if appending.");
  keys.add("optional","IN_RESTART","Read this file to continue an FISST simulation. "
           "If same as OUT_RESTART and you have not set the RESTART directive, the file will be backed-up and overwritten with new output."
           "If you do have the RESTART flag set and it is the same name as OUT_RESTART, this file will be appended.");
  keys.add("optional","OUT_OBSERVABLE","Output file putting weights needed to compute observables at different force values."
           "If you have the RESTART directive set (global or for FISST), this file will be appended to. "
           "Note that the header will be printed again if appending.");
  keys.add("optional","OBSERVABLE_FREQ","How often to write out observable weights (default=period).");
  keys.addFlag("FREEZE",false,"Fix bias weights at current level (only used for restarting).");
  keys.use("RESTART");
  keys.addOutputComponent("force2","default","squared value of force from the bias.");
  keys.addOutputComponent("_fbar","default", "For each named CV biased, there will be a corresponding output CV_fbar storing the current linear bias prefactor.");
}

FISST::FISST(const ActionOptions&ao):
  PLUMED_BIAS_INIT(ao),
  ncvs_(getNumberOfArguments()),
  current_avg_force_(ncvs_,0.0),
  center_(ncvs_,0.0),
  //change min_force and max_force to vectors if going to do more than one cv
  min_force_(0.0),
  max_force_(0.0),
  in_restart_name_(""),
  out_restart_name_(""),
  out_observable_name_(""),
  fmt_("%e"),
  b_freeze_(false),
  b_restart_(false),
  b_write_restart_(false),
  b_write_observable_(false),
  b_first_restart_sample_(true),
  n_interpolation_(0),
  n_samples_(0),
  initial_weight_rate_(0),
  initial_weight_dist_("UNIFORM"),
  period_(0),
  reset_period_(0),
  observable_freq_(0),
  kbt_(0.0),
  value_force2_(NULL)
{
  if(ncvs_==0)
    error("Must specify at least one CV with ARG");

  //temporary
  if(ncvs_>1)
    error("FISST only supports using one CV right now");

  addComponent("force2");
  componentIsNotPeriodic("force2");
  value_force2_ = getPntrToComponent("force2");

  for(unsigned int i = 0; i<ncvs_; i++) {
    std::string comp = getPntrToArgument(i)->getName() + "_fbar";
    addComponent(comp);
    componentIsNotPeriodic(comp);
  }

  parseVector("CENTER",center_);
  //change min_force and max_force to vectors if going to do more than one cv
  parse("MIN_FORCE",min_force_);
  parse("MAX_FORCE",max_force_);
  parse("PERIOD",period_);
  parse("RESET_PERIOD",reset_period_);
  parse("INITIAL_WEIGHT_DIST",initial_weight_dist_);
  parse("INITIAL_WEIGHT_RATE",initial_weight_rate_);
  parse("OBSERVABLE_FREQ",observable_freq_);
  parse("NINTERPOLATE",n_interpolation_);
  parseFlag("FREEZE",b_freeze_);
  parse("KBT",kbt_);
  parse("RESTART_FMT", fmt_);
  fmt_ = " " + fmt_;//add space since parse strips them
  parse("OUT_RESTART",out_restart_name_);
  parse("OUT_OBSERVABLE",out_observable_name_);
  parse("IN_RESTART",in_restart_name_);
  checkRead();

  if(center_.size() != ncvs_)
    error("Must have same number of CENTER arguments as ARG arguments");

  if(in_restart_name_ != "") {
    b_restart_ = true;
    log.printf("  reading simulation information from file: %s\n",in_restart_name_.c_str());
    readInRestart();
  } else {

    if(! kbt_ > 0.0)
      kbt_ = plumed.getAtoms().getKbT();

    //in driver, this results in kbt of 0
    if(kbt_ == 0) {
      error("  Unable to determine valid kBT. "
            "Could be because you are runnning from driver or MD didn't give temperature.\n"
            "Consider setting temperature manually with the KBT keyword.");
    }

    log.printf("  kBT = %f\n",kbt_);
    log.printf("  Updating with a time scale of %i steps\n",period_);

    log.printf("  Using centers for CVs of:");
    for(unsigned int i = 0; i< ncvs_; i++) {
      log.printf(" %f ",center_[i]);
    }
    log.printf("\n");
    observable_weight_.resize(n_interpolation_);
    for(unsigned int i = 0; i<n_interpolation_; i++) observable_weight_[i] = 1.0;

    forces_.resize(n_interpolation_);
    force_weight_.resize(n_interpolation_);
    //using code from the MIST project
    gauss_weight_.resize(n_interpolation_);
    legendre_compute_glr(n_interpolation_, &forces_[0], &gauss_weight_[0]);
    rescale(min_force_, max_force_, n_interpolation_, &forces_[0], &gauss_weight_[0]);

    log.printf("Using weight distribution %s with rate %f\n",initial_weight_dist_.c_str(),initial_weight_rate_);
    if(initial_weight_dist_ == "UNIFORM" ) {
      for(unsigned int i = 0; i<n_interpolation_; i++) force_weight_[i] = 1.0;
    }
    else if (initial_weight_dist_ == "EXP" ) {
      for(unsigned int i = 0; i<n_interpolation_; i++) force_weight_[i] = exp(-fabs(forces_[i])*initial_weight_rate_);
    }
    else if (initial_weight_dist_ == "GAUSS" ) {
      for(unsigned int i = 0; i<n_interpolation_; i++) force_weight_[i] = exp(-pow(forces_[i],2)*initial_weight_rate_);
    }
    else {
      error("  Specified weight distribution is not from the allowed list.");

    }

    partition_estimate_.resize(n_interpolation_);
    NormalizeForceWeights();
    double sum = 0.0;
    for(unsigned int i = 0; i<n_interpolation_; i++) {
      //setting partition estimate as 1/w_i
      partition_estimate_[i] = 1/force_weight_[i];
      log.printf("force/gauss weight/force_weight: %i %f %f %f\n",i,forces_[i],gauss_weight_[i],force_weight_[i]);
      sum+=gauss_weight_[i]*force_weight_[i];
    }
    log.printf("--Sum_i w_i g_i: %f\n",sum);

  }

  //set inverse temperature
  beta_ = 1/kbt_;

  if(b_freeze_ && b_restart_) {
    log.printf("  freezing weights read in from the restart file\n");
  }

  if(out_restart_name_.length()>0) {
    log.printf("  writing restart information every %i steps to file %s with format %s\n",abs(period_),out_restart_name_.c_str(), fmt_.c_str());
    b_write_restart_ = true;
    setupOutRestart();
  }
  if(out_observable_name_.length()>0) {
    if(observable_freq_==0) observable_freq_ = period_;
    log.printf("  writing observable information every %i steps to file %s with format %s\n",observable_freq_,out_observable_name_.c_str(), fmt_.c_str());
    b_write_observable_ = true;
    setupOutObservable();
  }

  //add citation later:
  //log<<"  Bibliography "<<plumed.cite("")<<"\n";
}

void FISST::NormalizeForceWeights() {
  double denom = 0.0;

  for(unsigned i=0; i<n_interpolation_; i++)
    denom += gauss_weight_[i] * force_weight_[i];

  for(unsigned i=0; i<n_interpolation_; i++)
    force_weight_[i] /= denom;
}

void FISST::readInRestart() {
  in_restart_.open(in_restart_name_);

  if(in_restart_.FieldExist("kbt")) {
    in_restart_.scanField("kbt",kbt_);
  } else { error("No field 'kbt' in restart file"); }
  log.printf("  with kBT = %f\n",kbt_);

  if(in_restart_.FieldExist("period")) {
    in_restart_.scanField("period",period_);
  } else { error("No field 'period' in restart file"); }
  log.printf("  Updating every %i steps\n",period_);

//this one can be optional
  if(in_restart_.FieldExist("reset_period")) {
    in_restart_.scanField("reset_period",reset_period_);
  }
  log.printf("  Resetting statistics every %i steps\n",reset_period_);

  if(in_restart_.FieldExist("n_interpolation")) {
    in_restart_.scanField("n_interpolation",n_interpolation_);
  } else { error("No field 'n_interpolation' in restart file"); }

  if(in_restart_.FieldExist("min_force")) {
    in_restart_.scanField("min_force",min_force_);
  } else { error("No field 'min_force' in restart file"); }
  if(in_restart_.FieldExist("max_force")) {
    in_restart_.scanField("max_force",max_force_);
  } else { error("No field 'max_force' in restart file"); }
  log.printf("  with forces from min_force=%e to max_force=%e over %i bins\n",min_force_,max_force_,n_interpolation_);

  unsigned int N = 0;
  std::string cv_name;
  double tmp, time;

  while(in_restart_.scanField("time",time)) {
    in_restart_.scanField("nsamples",n_samples_);

    observable_weight_.resize(n_interpolation_);
    partition_estimate_.resize(n_interpolation_);
    force_weight_.resize(n_interpolation_);
    gauss_weight_.resize(n_interpolation_);
    forces_.resize(n_interpolation_);

    for(unsigned int i = 0; i<ncvs_; ++i) {
      cv_name = getPntrToArgument(i)->getName();
      in_restart_.scanField(cv_name,tmp);
      for(unsigned int j =0; j<n_interpolation_; ++j) {
        in_restart_.scanField(cv_name + "_f"+std::to_string(j),forces_[j]);
        in_restart_.scanField(cv_name + "_g"+std::to_string(j),gauss_weight_[j]);
        in_restart_.scanField(cv_name + "_w"+std::to_string(j),force_weight_[j]);
        in_restart_.scanField(cv_name + "_z"+std::to_string(j),partition_estimate_[j]);
      }
    }
    N++;

    in_restart_.scanField();
  }

  double sum = 0.0;
  for(unsigned int j =0; j<n_interpolation_; ++j) {
    //clear observable weight, which will be set later
    observable_weight_[j] = 1.0;

    //setting partition estimate as 1/w_i
    log.printf("force/gauss weight/force_weight: %i %e %e %e\n",j,forces_[j],gauss_weight_[j],force_weight_[j]);
    sum+=gauss_weight_[j]*force_weight_[j];
  }
  log.printf("--Sum_i w_i g_i: %f\n",sum);

  in_restart_.close();
}

void FISST::setupOutObservable() {
  out_observable_.link(*this);
  out_observable_.fmtField(fmt_);
  out_observable_.open(out_observable_name_);
  out_observable_.setHeavyFlush();

  out_observable_.addConstantField("kbt").printField("kbt",kbt_);
  out_observable_.addConstantField("n_interpolation").printField("n_interpolation",n_interpolation_);
  out_observable_.addConstantField("period").printField("period",period_);
  out_observable_.addConstantField("min_force").printField("min_force",min_force_);
  out_observable_.addConstantField("max_force").printField("max_force",max_force_);
}

void FISST::setupOutRestart() {
  out_restart_.link(*this);
  out_restart_.fmtField(fmt_);
  out_restart_.open(out_restart_name_);
  out_restart_.setHeavyFlush();

  out_restart_.addConstantField("kbt").printField("kbt",kbt_);
  out_restart_.addConstantField("n_interpolation").printField("n_interpolation",n_interpolation_);
  out_restart_.addConstantField("period").printField("period",period_);
  if(reset_period_>0) out_restart_.addConstantField("reset_period").printField("reset_period",reset_period_);
  out_restart_.addConstantField("min_force").printField("min_force",min_force_);
  out_restart_.addConstantField("max_force").printField("max_force",max_force_);
}

void FISST::writeOutRestart() {
  std::string cv_name;
  out_restart_.printField("time",getTimeStep()*getStep());
  out_restart_.printField("nsamples",n_samples_);

  for(unsigned int i = 0; i<ncvs_; ++i) {
    cv_name = getPntrToArgument(i)->getName();
    double Q_i = difference(i, center_[i], getArgument(i));
    out_restart_.printField(cv_name,Q_i);
    for(int j = 0; j < n_interpolation_; j++ ) {
//have to update this for multiple cvs
      out_restart_.printField(cv_name + "_f"+std::to_string(j),forces_[j]);
      out_restart_.printField(cv_name + "_g"+std::to_string(j),gauss_weight_[j]);
      out_restart_.printField(cv_name + "_w"+std::to_string(j),force_weight_[j]);
      out_restart_.printField(cv_name + "_z"+std::to_string(j),partition_estimate_[j]);
    }
  }
  out_restart_.printField();
}

void FISST::writeOutObservable() {
  std::string cv_name;
  out_observable_.printField("time",getTimeStep()*getStep());
  out_observable_.printField("nsamples",n_samples_);

  for(unsigned int i = 0; i<ncvs_; ++i) {
    cv_name = getPntrToArgument(i)->getName();
    double Q_i = difference(i, center_[i], getArgument(i));
    out_observable_.printField(cv_name,Q_i);
    out_observable_.printField(cv_name + "_fbar",current_avg_force_[i]);
    for(int j = 0; j < n_interpolation_; j++ ) {
//have to update this for multiple cvs
      out_observable_.printField(cv_name + "_f"+std::to_string(j),forces_[j]);
      out_observable_.printField(cv_name + "_ow"+std::to_string(j),observable_weight_[j]);
    }
  }
  out_observable_.printField();
}


void FISST::calculate() {
  if(getStep() == 0 ) {
    if(b_write_restart_) writeOutRestart();
    if(b_write_observable_) writeOutObservable();
  }

  if(! b_freeze_) {
    if(b_restart_ && b_first_restart_sample_) {
      //dont' update statistics if restarting and first sample
      b_first_restart_sample_ = false;
    }
    else {
      update_statistics();
    }
  }
  update_bias();
  apply_bias();

  //check about writing restart file
  if(getStep()>0 && getStep()%period_==0) {
    if(b_write_restart_) writeOutRestart();
  }
  if(getStep()>0 && getStep()%observable_freq_==0) {
    if(b_write_observable_) {
      compute_observable_weight();
      writeOutObservable();
    }
  }
}


void FISST::apply_bias() {
  //Compute linear force as in "restraint"
  double ene = 0, totf2 = 0, cv, m, f;

  for(unsigned int i = 0; i < ncvs_; ++i) {
    cv = difference(i, center_[i], getArgument(i));
    double fbar = current_avg_force_[i];
    ene -= fbar*cv;
    setOutputForce(i,fbar);
    totf2 += fbar*fbar;

    std::string fbar_name_ = getPntrToArgument(i)->getName() + "_fbar";
    Value* fbar_ = getPntrToComponent(fbar_name_);
    fbar_->set(fbar);
  };

  setBias(ene);
  value_force2_->set(totf2);
  //log.flush();
}

void FISST::update_statistics()  {
//get stride is for multiple time stepping
  double dt=getTimeStep()*getStride();
  double h = dt/(period_*getTimeStep());
  double fbar_denum_integral = 0.0;

  int step = getStep();
  if(reset_period_>0 && step>0 && step%reset_period_==0) {
    n_samples_=1;
  }
  else {
    n_samples_++;
  }
  double d_n_samples = (double)n_samples_;

  for(unsigned int i = 0; i < ncvs_; ++i) {
    double Q_i = difference(i, center_[i], getArgument(i));
    for(unsigned int j=0; j<n_interpolation_; j++)
    {
      //if multiple cvs, these need to be updated to have 2 columns
      double f_j = forces_[j];
      double w_j = force_weight_[j];
      double g_j = gauss_weight_[j];

      fbar_denum_integral += g_j * w_j * exp(beta_*f_j * Q_i);
    }

    for(unsigned int j=0; j<n_interpolation_; j++)
    {
      double f_j = forces_[j];
      double sample_weight = exp(beta_*f_j * Q_i) / fbar_denum_integral;

      partition_estimate_[j] = sample_weight/d_n_samples + partition_estimate_[j]*(d_n_samples-1)/(d_n_samples);

      double w_jn = force_weight_[j];
      double z_jn = partition_estimate_[j];

      double w_jp1 = (1.0 - h) * w_jn + h / z_jn;
      force_weight_[j] = w_jp1;
    }
  }

  // make sure that the weights are normalised
  NormalizeForceWeights();
}


void FISST::update_bias()
{
  for(unsigned int i = 0; i < ncvs_; ++i) {
    double Q_i = difference(i, center_[i], getArgument(i));
    double fbar_num_integral = 0.0;
    double fbar_denum_integral = 0.0;

    for(unsigned int j=0; j<n_interpolation_; j++ ) {
      double f_j = forces_[j];
      double w_j = force_weight_[j];
      double g_j = gauss_weight_[j];

      fbar_num_integral += g_j * f_j * w_j * exp(beta_*f_j*Q_i);
      fbar_denum_integral += g_j * w_j * exp(beta_*f_j*Q_i);
    }

    current_avg_force_[i] = fbar_num_integral/fbar_denum_integral;
  }
}

void FISST::compute_observable_weight() {
  double obs_num = (max_force_ - min_force_);

  for(unsigned int i = 0; i < ncvs_; ++i) {
    double Q_i = difference(i, center_[i], getArgument(i));

    for(unsigned int j=0; j<n_interpolation_; j++ ) {
      double z_j = partition_estimate_[j];
      double f_j = forces_[j];
      double denum_integral = 0.0;

      for( unsigned int k=0; k<n_interpolation_; k++ ) {
        double f_k = forces_[k];
        double w_k = force_weight_[k];
        double g_k = gauss_weight_[k];

        denum_integral += g_k * w_k * exp(beta_*(f_k-f_j)*Q_i);
      }
      observable_weight_[j] = obs_num/(denum_integral*z_j);
    }
  }
}



void FISST::update() {
  //pass
}

FISST::~FISST() {
  out_restart_.close();
  out_observable_.close();
}

void FISST::turnOnDerivatives() {
  // do nothing
  // this is to avoid errors triggered when a bias is used as a CV
  // (This is done in ExtendedLagrangian.cpp)
}


}
}//close the 2 namespaces
