/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2017 of Glen Hocky and Andrew White

The eds module is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The eds module is distributed in the hope that it will be useful,
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



#include <iostream>


using namespace PLMD;
using namespace bias;

//namespace is lowercase to match
//module names being all lowercase

namespace PLMD {
namespace eds {

//+PLUMEDOC EDSMOD_BIAS EDS
/*
Add a linear bias on a set of observables.

This force is the same as the linear part of the bias in \ref
RESTRAINT, but this bias has the ability to compute prefactors
adaptively using the scheme of White and Voth \cite white2014efficient
in order to match target observable values for a set of CVs. You can
see a tutorial on EDS specifically for biasing coordination number at
<a
href="http://thewhitelab.org/Blog/tutorial/2017/05/10/lammps-coordination-number-tutorial/">
Andrew White's webpage</a>.

The addition to the potential is of the form
\f[
  \sum_i \frac{\alpha_i}{s_i} x_i
\f]

where for CV \f$x_i\f$, a coupling constant \f${\alpha}_i\f$ is determined
adaptively or set by the user to match a target value for
\f$x_i\f$. \f$s_i\f$ is a scale parameter, which by default is set to
the target value. It may also be set separately.

\warning
It is not possible to set the target value of the observable
to zero with the default value of \f$s_i\f$ as this will cause a
divide-by-zero error. Instead, set \f$s_i=1\f$ or modify the CV so the
desired target value is no longer zero.

Notice that a similar method is available as \ref MAXENT, although with different features and using a different optimization algorithm.

\par Examples

The following input for a harmonic oscillator of two beads will
adaptively find a linear bias to change the mean and variance to the
target values. The PRINT line shows how to access the value of the
coupling constants.

\plumedfile
dist: DISTANCE ATOMS=1,2
# this is the squared of the distance
dist2: COMBINE ARG=dist POWERS=2 PERIODIC=NO

#bias mean and variance
eds: EDS ARG=dist,dist2 CENTER=2.0,1.0 PERIOD=50000 TEMP=1.0
PRINT ARG=dist,dist2,eds.dist_coupling,eds.dist2_coupling,eds.bias,eds.force2 FILE=colvars.dat STRIDE=100
\endplumedfile

Rather than trying to find the coupling constants adaptively, one can ramp up to a constant value.
\plumedfile
#ramp couplings from 0,0 to -1,1 over 50000 steps
eds: EDS ARG=dist,dist2 CENTER=2.0,1.0 FIXED=-1,1 RAMP PERIOD=50000 TEMP=1.0

#same as above, except starting at -0.5,0.5 rather than default of 0,0
eds: EDS ARG=dist,dist2 CENTER=2.0,1.0 FIXED=-1,1 INIT=-0.5,0.5 RAMP PERIOD=50000 TEMP=1.0
\endplumedfile

A restart file can be added to dump information needed to restart/continue simulation using these parameters every PERIOD.
\plumedfile
#add the option to write to a restart file
eds: EDS ARG=dist,dist2 CENTER=2.0,1.0 PERIOD=50000 TEMP=1.0 OUT_RESTART=restart.dat
\endplumedfile

Read in a previous restart file. Adding RESTART flag makes output append
\plumedfile
eds: EDS ARG=dist,dist2 CENTER=2.0,1.0 PERIOD=50000 TEMP=1.0 IN_RESTART=restart.dat RESTART
\endplumedfile

Read in a previous restart file and freeze the bias at the final level from the previous simulation
\plumedfile
eds: EDS ARG=dist,dist2 CENTER=2.0,1.0 TEMP=1.0 IN_RESTART=restart.dat FREEZE
\endplumedfile

Read in a previous restart file and freeze the bias at the mean from the previous simulation
\plumedfile
eds: EDS ARG=dist,dist2 CENTER=2.0,1.0 TEMP=1.0 IN_RESTART=restart.dat FREEZE MEAN
\endplumedfile

Read in a previous restart file and continue the bias, but use the mean from the previous run as the starting point
\plumedfile
eds: EDS ARG=dist,dist2 CENTER=2.0,1.0 PERIOD=50000 TEMP=1.0 IN_RESTART=restart.dat MEAN
\endplumedfile


*/
//+ENDPLUMEDOC

class EDS : public Bias {


private:
  /*We will get this and store it once, since on-the-fly changing number of CVs will be fatal*/
  const unsigned int ncvs_;
  std::vector<double> center_;
  std::vector<Value*> center_values_;
  std::vector<double> scale_;
  std::vector<double> current_coupling_;
  std::vector<double> set_coupling_;
  std::vector<double> target_coupling_;
  std::vector<double> max_coupling_range_;
  std::vector<double> max_coupling_grad_;
  std::vector<double> coupling_rate_;
  std::vector<double> coupling_accum_;
  std::vector<double> means_;
  std::vector<double> ssds_;
  std::vector<double> step_size_;
  std::vector<Value*> out_coupling_;
  Matrix<double> covar_;
  std::string in_restart_name_;
  std::string out_restart_name_;
  std::string fmt_;
  OFile out_restart_;
  IFile in_restart_;
  bool b_c_values_;
  bool b_adaptive_;
  bool b_freeze_;
  bool b_equil_;
  bool b_ramp_;
  bool b_covar_;
  bool b_restart_;
  bool b_write_restart_;
  bool b_hard_c_range_;
  int seed_;
  int update_period_;
  int avg_coupling_count_;
  int update_calls_;
  double kbt_;
  double c_range_increase_f_;
  double multi_prop_;
  Random rand_;
  Value* value_force2_;

  /*read input restart. b_mean sets if we use mean or final value for freeze*/
  void readInRestart(const bool b_mean);
  /*setup output restart*/
  void setupOutRestart();
  /*write output restart*/
  void writeOutRestart();
  void update_statistics();
  void calc_covar_step_size();
  void calc_ssd_step_size();
  void reset_statistics();
  void update_bias();
  void apply_bias();

public:
  explicit EDS(const ActionOptions&);
  void calculate();
  void update();
  void turnOnDerivatives();
  static void registerKeywords(Keywords& keys);
  ~EDS();
};

PLUMED_REGISTER_ACTION(EDS,"EDS")

void EDS::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("optional","CENTER","The desired centers (equilibrium values) which will be sought during the adaptive linear biasing. This is for fixed values");
  keys.add("optional","CENTER_ARG","The desired centers (equilibrium values) which will be sought during the adaptive linear biasing. "
           "CENTER_ARG is for calculated centers, e.g. from a CV or analysis. ");

  keys.add("optional","PERIOD","Steps over which to adjust bias for adaptive or ramping");

  keys.add("compulsory","RANGE","3.0","The largest magnitude of the force constant which one expects (in kBT) for each CV based");
  keys.add("compulsory","SEED","0","Seed for random order of changing bias");
  keys.add("compulsory","INIT","0","Starting value for coupling constant");
  keys.add("compulsory","FIXED","0","Fixed target values for coupling constant. Non-adaptive.");
  keys.add("optional","BIAS_SCALE","A divisor to set the units of the bias. "
           "If not set, this will be the experimental value by default (as is done in White and Voth 2014).");
  keys.add("optional","TEMP","The system temperature. If not provided will be taken from MD code (if available)");
  keys.add("optional","MULTI_PROP","What proportion of dimensions to update at each step. "
           "Must be in interval [1,0), where 1 indicates all and any other indicates a stochastic update. "
           "If not set, default is 1 / N, where N is the number of CVs. ");
  keys.add("optional","RESTART_FMT","the format that should be used to output real numbers in EDS restarts");
  keys.add("optional","OUT_RESTART","Output file for all information needed to continue EDS simulation. "
           "If you have the RESTART directive set (global or for EDS), this file will be appended to. "
           "Note that the header will be printed again if appending.");
  keys.add("optional","IN_RESTART","Read this file to continue an EDS simulation. "
           "If same as OUT_RESTART and you have not set the RESTART directive, the file will be backed-up and overwritten with new output. "
           "If you do have the RESTART flag set and it is the same name as OUT_RESTART, this file will be appended.");

  keys.addFlag("RAMP",false,"Slowly increase bias constant to a fixed value");
  keys.addFlag("COVAR",false,"Utilize the covariance matrix when updating the bias. Default Off, but may be enabled due to other options");
  keys.addFlag("FREEZE",false,"Fix bias at current level (only used for restarting).");
  keys.addFlag("MEAN",false,"Instead of using final bias level from restart, use average. Can only be used in conjunction with FREEZE");

  keys.use("RESTART");

  keys.addOutputComponent("force2","default","squared value of force from the bias");
  keys.addOutputComponent("_coupling","default", "For each named CV biased, there will be a corresponding output CV_coupling storing the current linear bias prefactor.");
}

EDS::EDS(const ActionOptions&ao):
  PLUMED_BIAS_INIT(ao),
  ncvs_(getNumberOfArguments()),
  scale_(ncvs_,0.0),
  current_coupling_(ncvs_,0.0),
  set_coupling_(ncvs_,0.0),
  target_coupling_(ncvs_,0.0),
  max_coupling_range_(ncvs_,3.0),
  max_coupling_grad_(ncvs_,0.0),
  coupling_rate_(ncvs_,1.0),
  coupling_accum_(ncvs_,0.0),
  means_(ncvs_,0.0),
  step_size_(ncvs_,0.0),
  out_coupling_(ncvs_,NULL),
  in_restart_name_(""),
  out_restart_name_(""),
  fmt_("%f"),
  b_adaptive_(true),
  b_freeze_(false),
  b_equil_(true),
  b_ramp_(false),
  b_covar_(false),
  b_restart_(false),
  b_write_restart_(false),
  b_hard_c_range_(false),
  seed_(0),
  update_period_(0),
  avg_coupling_count_(1),
  update_calls_(0),
  kbt_(0.0),
  c_range_increase_f_(1.25),
  multi_prop_(-1.0),
  value_force2_(NULL)
{
  double temp=-1.0;
  bool b_mean=false;

  addComponent("force2");
  componentIsNotPeriodic("force2");
  value_force2_ = getPntrToComponent("force2");

  for(unsigned int i = 0; i<ncvs_; i++) {
    std::string comp = getPntrToArgument(i)->getName() + "_coupling";
    addComponent(comp);
    componentIsNotPeriodic(comp);
    out_coupling_[i]=getPntrToComponent(comp);
  }

  parseVector("CENTER",center_);
  parseArgumentList("CENTER_ARG",center_values_);
  parseVector("BIAS_SCALE", scale_);
  parseVector("RANGE",max_coupling_range_);
  parseVector("FIXED",target_coupling_);
  parseVector("INIT",set_coupling_);
  parse("PERIOD",update_period_);
  parse("TEMP",temp);
  parse("SEED",seed_);
  parse("MULTI_PROP",multi_prop_);
  parse("RESTART_FMT", fmt_);
  fmt_ = " " + fmt_;//add space since parse strips them
  parse("OUT_RESTART",out_restart_name_);
  parseFlag("RAMP",b_ramp_);
  parseFlag("FREEZE",b_freeze_);
  parseFlag("MEAN",b_mean);
  parseFlag("COVAR",b_covar_);
  parse("IN_RESTART",in_restart_name_);
  checkRead();

  /*
   * Things that are different when usnig changing centers:
   * 1. Scale
   * 2. The log file
   * 3. Reading Restarts
   */

  if(center_.size() == 0) {
    if(center_values_.size() == 0)
      error("Must set either CENTER or CENTER_ARG");
    else if(center_values_.size() != ncvs_)
      error("CENTER_ARG must contain the same number of variables as ARG");
    b_c_values_ = true;
    center_.resize(ncvs_);
    log.printf("  EDS will use possibly varying centers\n");
  } else {
    if(center_.size() != ncvs_)
      error("Must have same number of CENTER arguments as ARG arguments");
    else if(center_values_.size() != 0)
      error("You can only set CENTER or CENTER_ARG. Not both");
    b_c_values_ = false;
    log.printf("  EDS will use fixed centers\n");
  }



  log.printf("  setting scaling:");
  if(scale_.size() > 0  && scale_.size() < ncvs_) {
    error("the number of BIAS_SCALE values be the same as number of CVs");
  } else if(scale_.size() == 0 && b_c_values_) {
    log.printf(" Setting SCALE to be 1 for all CVs\n");
    scale_.resize(ncvs_);
    for(unsigned int i = 0; i < ncvs_; ++i)
      scale_[i] = 1;
  } else if(scale_.size() == 0 && !b_c_values_) {
    log.printf(" (default) ");

    scale_.resize(ncvs_);
    for(unsigned int i = 0; i < scale_.size(); i++) {
      if(center_[i]==0)
        error("BIAS_SCALE parameter has been set to CENTER value of 0 (as is default). This will divide by 0, so giving up. See doc for EDS bias");
      scale_[i] = center_[i];
    }
  } else {
    for(unsigned int i = 0; i < scale_.size(); i++)
      log.printf(" %f",scale_[i]);
  }
  log.printf("\n");


  if(b_covar_) {
    log.printf("  EDS will utilize covariance matrix for update steps\n");
    covar_.resize(ncvs_, ncvs_);
  } else {
    log.printf("  EDS will utilize variance for update steps\n");
    ssds_.resize(ncvs_);
  }


  if (b_mean == true and b_freeze_ == false) {
    error("EDS keyworkd MEAN can only be used along with keyword FREEZE");
  }

  if(in_restart_name_ != "") {
    b_restart_ = true;
    log.printf("  reading simulation information from file: %s\n",in_restart_name_.c_str());
    readInRestart(b_mean);
  } else {

    if(temp>=0.0) kbt_=plumed.getAtoms().getKBoltzmann()*temp;
    else kbt_ = plumed.getAtoms().getKbT();

    //in driver, this results in kbt of 0
    if(kbt_ == 0) {
      error("  Unable to determine valid kBT. "
            "Could be because you are runnning from driver or MD didn't give temperature.\n"
            "Consider setting temperature manually with the TEMP keyword.");
      kbt_ = 1;
    }

    log.printf("  kBT = %f\n",kbt_);
    log.printf("  Updating every %i steps\n",update_period_);

    if(!b_c_values_) {
      log.printf("  with centers:");
      for(unsigned int i = 0; i< ncvs_; i++) {
        log.printf(" %f ",center_[i]);
      }
    } else {
      log.printf("  with actions centers:");
      for(unsigned int i = 0; i< ncvs_; i++) {
        log.printf(" %s ",center_values_[i]->getName().c_str());
        //add dependency on these actions
        addDependency(center_values_[i]->getPntrToAction());
      }
    }

    log.printf("\n  with initial ranges / rates:\n");
    for(unsigned int i = 0; i<max_coupling_range_.size(); i++) {
      //this is just an empirical guess. Bigger range, bigger grads. Less frequent updates, bigger changes
      max_coupling_range_[i]*=kbt_;
      max_coupling_grad_[i] = max_coupling_range_[i]*update_period_/100.;
      log.printf("    %f / %f\n",max_coupling_range_[i],max_coupling_grad_[i]);
    }

    if(seed_>0) {
      log.printf("  setting random seed = %i",seed_);
      rand_.setSeed(seed_);
    }

    for(unsigned int i = 0; i<ncvs_; ++i) if(target_coupling_[i]!=0.0) b_adaptive_=false;

    if(!b_adaptive_) {
      if(b_ramp_) {
        log.printf("  ramping up coupling constants over %i steps\n",update_period_);
      }

      log.printf("  with starting coupling constants");
      for(unsigned int i = 0; i<set_coupling_.size(); i++) log.printf(" %f",set_coupling_[i]);
      log.printf("\n");
      log.printf("  and final coupling constants");
      for(unsigned int i = 0; i<target_coupling_.size(); i++) log.printf(" %f",target_coupling_[i]);
      log.printf("\n");
    }

    //now do setup
    if(b_ramp_) {
      update_period_*=-1;
    }

    for(unsigned int i = 0; i<set_coupling_.size(); i++) current_coupling_[i] = set_coupling_[i];

    // if b_adaptive_, then first half will be used for equilibrating and second half for statistics
    if(update_period_>0) {
      update_period_ /= 2;
    }


  }

  if(b_freeze_) {
    b_adaptive_=false;
    update_period_ = 0;
    if (b_mean) {
      log.printf("  freezing bias at the average level from the restart file\n");
    } else {
      log.printf("  freezing bias at current level\n");
    }
  }

  if(multi_prop_ == -1.0) {
    log.printf("  Will update each dimension stochastically with probability 1 / number of CVs\n");
    multi_prop_ = 1.0 / ncvs_;
  } else if(multi_prop_ > 0 && multi_prop_ <= 1.0) {
    log.printf("  Will update each dimension stochastically with probability %f\n", multi_prop_);
  } else {
    error("  MULTI_PROP must be between 0 and 1\n");
  }

  if(out_restart_name_.length()>0) {
    log.printf("  writing restart information every %i steps to file %s with format %s\n",abs(update_period_),out_restart_name_.c_str(), fmt_.c_str());
    b_write_restart_ = true;
    setupOutRestart();
  }

  log<<"  Bibliography "<<plumed.cite("White and Voth, J. Chem. Theory Comput. 10 (8), 3023-3030 (2014)")<<"\n";
}

void EDS::readInRestart(const bool b_mean) {
  int adaptive_i;

  in_restart_.open(in_restart_name_);

  if(in_restart_.FieldExist("kbt")) {
    in_restart_.scanField("kbt",kbt_);
  } else { error("No field 'kbt' in restart file"); }
  log.printf("  with kBT = %f\n",kbt_);

  if(in_restart_.FieldExist("update_period")) {
    in_restart_.scanField("update_period",update_period_);
  } else { error("No field 'update_period' in restart file"); }
  log.printf("  Updating every %i steps\n",update_period_);

  if(in_restart_.FieldExist("adaptive")) {
    //note, no version of scanField for boolean
    in_restart_.scanField("adaptive",adaptive_i);
  } else { error("No field 'adaptive' in restart file"); }
  b_adaptive_ = bool(adaptive_i);

  if(in_restart_.FieldExist("seed")) {
    in_restart_.scanField("seed",seed_);
  } else { error("No field 'seed' in restart file"); }
  if(seed_>0) {
    log.printf("  setting random seed = %i",seed_);
    rand_.setSeed(seed_);
  }

  double time, tmp;
  std::vector<double> avg_bias = std::vector<double>(center_.size());
  unsigned int N = 0;
  std::string cv_name;

  while(in_restart_.scanField("time",time)) {

    for(unsigned int i = 0; i<ncvs_; ++i) {
      cv_name = getPntrToArgument(i)->getName();
      in_restart_.scanField(cv_name + "_center", set_coupling_[i]);
      in_restart_.scanField(cv_name + "_set", set_coupling_[i]);
      in_restart_.scanField(cv_name + "_target",target_coupling_[i]);
      in_restart_.scanField(cv_name + "_coupling",current_coupling_[i]);
      in_restart_.scanField(cv_name + "_maxrange",max_coupling_range_[i]);
      in_restart_.scanField(cv_name + "_maxgrad",max_coupling_grad_[i]);
      in_restart_.scanField(cv_name + "_accum",coupling_accum_[i]);
      in_restart_.scanField(cv_name + "_mean",means_[i]);
      //unused due to difference between covar/nocovar
      in_restart_.scanField(cv_name + "_std",tmp);

      avg_bias[i] += current_coupling_[i];
    }
    N++;

    in_restart_.scanField();
  }


  log.printf("  with centers:");
  for(unsigned int i = 0; i<center_.size(); i++) {
    log.printf(" %f",center_[i]);
  }
  log.printf("\n  and scaling:");
  for(unsigned int i = 0; i<scale_.size(); i++) {
    log.printf(" %f",scale_[i]);
  }

  log.printf("\n  with initial ranges / rates:\n");
  for(unsigned int i = 0; i<max_coupling_range_.size(); i++) {
    log.printf("    %f / %f\n",max_coupling_range_[i],max_coupling_grad_[i]);
  }

  if(!b_adaptive_ && update_period_<0) {
    log.printf("  ramping up coupling constants over %i steps\n",-update_period_);
  }

  if(b_mean) {
    log.printf("Loaded in averages for coupling constants...\n");
    for(unsigned int i = 0; i<current_coupling_.size(); i++) current_coupling_[i] = avg_bias[i] / N;
    for(unsigned int i = 0; i<current_coupling_.size(); i++) set_coupling_[i] = avg_bias[i] / N;
  }

  log.printf("  with current coupling constants:\n    ");
  for(unsigned int i = 0; i<current_coupling_.size(); i++) log.printf(" %f",current_coupling_[i]);
  log.printf("\n");
  log.printf("  with initial coupling constants:\n    ");
  for(unsigned int i = 0; i<set_coupling_.size(); i++) log.printf(" %f",set_coupling_[i]);
  log.printf("\n");
  log.printf("  and final coupling constants:\n    ");
  for(unsigned int i = 0; i<target_coupling_.size(); i++) log.printf(" %f",target_coupling_[i]);
  log.printf("\n");

  in_restart_.close();
}

void EDS::setupOutRestart() {
  out_restart_.link(*this);
  out_restart_.fmtField(fmt_);
  out_restart_.open(out_restart_name_);
  out_restart_.setHeavyFlush();

  out_restart_.addConstantField("adaptive").printField("adaptive",b_adaptive_);
  out_restart_.addConstantField("update_period").printField("update_period",update_period_);
  out_restart_.addConstantField("seed").printField("seed",seed_);
  out_restart_.addConstantField("kbt").printField("kbt",kbt_);

}

void EDS::writeOutRestart() {
  std::string cv_name;
  out_restart_.printField("time",getTimeStep()*getStep());

  for(unsigned int i = 0; i<ncvs_; ++i) {
    cv_name = getPntrToArgument(i)->getName();
    out_restart_.printField(cv_name + "_center",center_[i]);
    out_restart_.printField(cv_name + "_set",set_coupling_[i]);
    out_restart_.printField(cv_name + "_target",target_coupling_[i]);
    out_restart_.printField(cv_name + "_coupling",current_coupling_[i]);
    out_restart_.printField(cv_name + "_maxrange",max_coupling_range_[i]);
    out_restart_.printField(cv_name + "_maxgrad",max_coupling_grad_[i]);
    out_restart_.printField(cv_name + "_accum",coupling_accum_[i]);
    out_restart_.printField(cv_name + "_mean",means_[i]);
    if(!b_covar_)
      out_restart_.printField(cv_name + "_std",ssds_[i] / (fmax(1, update_calls_ - 1)));
    else
      out_restart_.printField(cv_name + "_std",covar_(i,i) / (fmax(1, update_calls_ - 1)));

  }
  out_restart_.printField();
}



void EDS::calculate() {

  //get center values from action if necessary
  if(b_c_values_)
    for(unsigned int i = 0; i < ncvs_; ++i)
      center_[i] = center_values_[i]->get();

  apply_bias();

  //adjust parameters according to EDS recipe
  update_calls_++;

  //check if we're ramping or doing normal updates and then restart if needed. The ramping check
  //is complicated because we could be frozen, finished ramping or not ramping.
  //The + 2 is so we have an extra line showing that the bias isn't changing (for my sanity and yours)
  if( b_write_restart_) {
    if(getStep() == 0 ||
        ( (update_period_ < 0 && !b_freeze_ && update_calls_ <= fabs(update_period_) + 2) ||
          (update_period_ > 0 && update_calls_ % update_period_ == 0 ) ) )
      writeOutRestart();
  }

  int b_finished_equil_flag = 1;

  //assume forces already applied and saved


  //are we ramping to a constant value and not done equilibrating?
  if(update_period_ < 0) {
    if(update_calls_ <= fabs(update_period_) && !b_freeze_) {
      for(unsigned int i = 0; i < ncvs_; ++i)
        current_coupling_[i] += (target_coupling_[i]-set_coupling_[i])/fabs(update_period_);
    }
    //make sure we don't reset update calls
    b_finished_equil_flag = 0;
  } else if(update_period_ == 0) { //do we have a no-update case?
    //not updating
    //pass
  } else if(!b_equil_) {
    //if we aren't wating for the bias to equilibrate, collect data
    update_statistics();
  } else {
    // equilibrating
    //check if we've reached the setpoint
    for(unsigned int i = 0; i < ncvs_; ++i) {
      if(coupling_rate_[i] == 0 || pow(current_coupling_[i] - set_coupling_[i],2) < pow(coupling_rate_[i],2)) {
        b_finished_equil_flag &= 1;
      }
      else {
        current_coupling_[i] += coupling_rate_[i];
        b_finished_equil_flag = 0;
      }
    }
  }

  //Update max coupling range if not hard
  if(!b_hard_c_range_) {
    for(unsigned int i = 0; i < ncvs_; ++i) {
      if(fabs(current_coupling_[i])>max_coupling_range_[i]) {
        max_coupling_range_[i]*=c_range_increase_f_;
        max_coupling_grad_[i]*=c_range_increase_f_;
      }
    }
  }

  //reduce all the flags
  if(b_equil_ && b_finished_equil_flag) {
    b_equil_ = false;
    update_calls_ = 0;
  }

  //Now we update coupling constant, if necessary
  if(!b_equil_ && update_period_ > 0 && update_calls_ == update_period_ && !b_freeze_) {
    update_bias();
    update_calls_ = 0;
    avg_coupling_count_++;
    b_equil_ = true; //back to equilibration now
  } //close update if

  //pass couplings out so they are accessible
  for(unsigned int i = 0; i<ncvs_; ++i) {
    out_coupling_[i]->set(current_coupling_[i]);
  }


}

void EDS::apply_bias() {
  //Compute linear force as in "restraint"
  double ene = 0, totf2 = 0, cv, m, f;

  for(unsigned int i = 0; i < ncvs_; ++i) {
    cv = difference(i, center_[i], getArgument(i));
    m = current_coupling_[i];
    f = -m;
    ene += m*cv;
    setOutputForce(i,f);
    totf2 += f*f;
  };

  setBias(ene);
  value_force2_->set(totf2);

}

void EDS::update_statistics()  {
  double s;
  std::vector<double> deltas(ncvs_);
  //Welford, West, and Hanso online variance method
  for(unsigned int i = 0; i < ncvs_; ++i)  {
    deltas[i] = difference(i,means_[i],getArgument(i));
    means_[i] += deltas[i]/fmax(1,update_calls_);
    if(!b_covar_)
      ssds_[i] += deltas[i]*difference(i,means_[i],getArgument(i));
  }
  if(b_covar_) {
    for(unsigned int i = 0; i < ncvs_; ++i) {
      for(unsigned int j = i; j < ncvs_; ++j) {
        s = (update_calls_ - 1) * deltas[i] * deltas[j] / update_calls_ / update_calls_ - covar_(i,j) / update_calls_;
        covar_(i,j) += s;
        //do this so we don't double count
        covar_(j,i) = covar_(i,j);
      }
    }
  }
}

void EDS::reset_statistics() {
  for(unsigned int i = 0; i < ncvs_; ++i)  {
    means_[i] = 0;
    if(!b_covar_)
      ssds_[i] = 0;
  }
  if(b_covar_)
    for(unsigned int i = 0; i < ncvs_; ++i)
      for(unsigned int j = 0; j < ncvs_; ++j)
        covar_(i,j) = 0;
}

void EDS::calc_covar_step_size() {
  //calulcate step size
  //uses scale here, which by default is center
  double tmp;
  for(unsigned int i = 0; i< ncvs_; ++i) {
    tmp = 0;
    for(unsigned int j = 0; j < ncvs_; ++j)
      tmp += difference(i, center_[i], means_[i]) * covar_(i,j);
    step_size_[i] = 2 * tmp / kbt_ / scale_[i] * update_calls_ / (update_calls_ - 1);
  }

}

void EDS::calc_ssd_step_size() {
  double tmp;
  for(unsigned int i = 0; i< ncvs_; ++i) {
    tmp = 2. * difference(i, center_[i], means_[i]) * ssds_[i] / (update_calls_ - 1);
    step_size_[i] = tmp / kbt_/scale_[i];
  }
}

void EDS::update_bias()
{
  if(b_covar_)
    calc_covar_step_size();
  else
    calc_ssd_step_size();

  for(unsigned int i = 0; i< ncvs_; ++i) {

    //check if the step_size exceeds maximum possible gradient
    step_size_[i] = copysign(fmin(fabs(step_size_[i]), max_coupling_grad_[i]), step_size_[i]);

    //reset means/vars
    reset_statistics();

    //multidimesional stochastic step
    if(ncvs_ == 1 || (rand_.RandU01() < (multi_prop_) ) ) {
      coupling_accum_[i] += step_size_[i] * step_size_[i];

      //equation 5 in White and Voth, JCTC 2014
      //no negative sign because it's in step_size
      set_coupling_[i] += max_coupling_range_[i]/sqrt(coupling_accum_[i])*step_size_[i];
      coupling_rate_[i] = (set_coupling_[i]-current_coupling_[i])/update_period_;

    } else {
      //do not change the bias
      coupling_rate_[i] = 0;
    }
  }
}



void EDS::update() {
  //pass
}

EDS::~EDS() {
  out_restart_.close();
}

void EDS::turnOnDerivatives() {
  // do nothing
  // this is to avoid errors triggered when a bias is used as a CV
  // (This is done in ExtendedLagrangian.cpp)
}


}
}//close the 2 namespaces
