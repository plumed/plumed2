/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#include "bias/Bias.h"
#include "core/PlumedMain.h"
#include "core/ActionRegister.h"
#include "core/Atoms.h"
#include "tools/Communicator.h"
#include "tools/File.h"

namespace PLMD {
namespace opes {

//+PLUMEDOC BIAS OPES_UMBRELLAS_TEST
/*
\par Examples

OPES_UMBRELLAS_TEST ...
  LABEL=test
  ARG=cv
  PACE=50
  TEMP=300
  SIGMA=0
  MIN_CV=0
  MAX_CV=1
... OPES_UMBRELLAS_TEST


*/
//+ENDPLUMEDOC

class OPESumbrellas_test : public bias::Bias {

private:
  bool isFirstStep_;
  bool afterCalculate_;
  unsigned NumParallel_;
  unsigned rank_;
  unsigned NumWalkers_;
  unsigned long counter_;

  unsigned stride_;
  unsigned obs_steps_;
  std::vector<double> obs_cv_;
  unsigned tot_umbrellas_;
  unsigned P0_contribution_;

  double beta_;
  double sigma_;
  std::vector<double> center_;
  std::vector<double> deltaF_;
  double rct_;
  double my_rct_;
  double barrier_;
  double current_bias_;

  bool calc_work_;
  double work_;
  std::vector<double> old_deltaF_;

  OFile deltaFsOfile_;
  unsigned print_stride_;

public:
  OPESumbrellas_test(const ActionOptions&);
  void calculate() override;
  void update() override;
  void init_from_obs();
  unsigned estimate_steps(const double,const double,const std::vector<double>&,const std::string) const;
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(OPESumbrellas_test,"OPES_UMBRELLAS_TEST")

void OPESumbrellas_test::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","TEMP","-1","temperature. If not specified tries to get it from MD engine");
  keys.add("compulsory","PACE","10","how often the bias is updated");
  keys.add("compulsory","OBSERVATION_STEPS","1","number of unbiased initial steps to collect statistics for initial deltaFs guess");
  keys.add("optional","BARRIER","a guess of the free energy barrier to be overcome (better to stay higher than lower)");
//umbrella stuff
  keys.add("compulsory","SIGMA","sigma of the umbrella Gaussians");
  keys.add("compulsory","MIN_CV","the minimum of the CV range to be explored");
  keys.add("compulsory","MAX_CV","the maximum of the CV range to be explored");
//deltaFs file
  keys.add("compulsory","FILE","DELTAFS","a file with the estimate of the relative \\f$\\Delta F\\f$ for each component of the target");
  keys.add("optional","PRINT_STRIDE","( default=100 ) stride for printing to DELTAFS file");
  keys.add("optional","FMT","specify format for DELTAFS file");
//miscellaneous
  keys.addFlag("ADD_P0",false,"add also P0 to the target distribution, useful to make sure the target is broader than P0");
  keys.addFlag("CALC_WORK",false,"calculate the work done by the bias between each update");
  keys.addFlag("WALKERS_MPI",false,"switch on MPI version of multiple walkers");
  keys.addFlag("SERIAL",false,"perform calculations in serial");
  keys.use("RESTART");

//output components
  componentsAreNotOptional(keys);
  keys.addOutputComponent("rct","WALKERS_MPI","single walker estimate of \\f$c(t)\\f$: \\f$\\frac{1}{\\beta}\\log \\lange e^{\\beta V} \\rangle\\f$, all walkers should converge to the same value");
  keys.addOutputComponent("work","CALC_WORK","work done by the bias between each update"); //calculating this maybe is only a useless overhead...
}

OPESumbrellas_test::OPESumbrellas_test(const ActionOptions&ao)
  : PLUMED_BIAS_INIT(ao)
  , isFirstStep_(true)
  , afterCalculate_(false)
  , counter_(0)
  , rct_(0)
  , my_rct_(0)
  , work_(0)
  , print_stride_(100)
{
  plumed_massert(getNumberOfArguments()==1,"only one cv is supported");

//set beta_
  const double Kb=plumed.getAtoms().getKBoltzmann();
  double KbT=plumed.getAtoms().getKbT();
  double temp=-1;
  parse("TEMP",temp);
  if(temp>0)
  {
    if(KbT>0 && std::abs(KbT-Kb*temp)>1e-4)
      log.printf(" +++ WARNING +++ using TEMP=%g while MD engine uses %g\n",temp,KbT/Kb);
    KbT=Kb*temp;
  }
  plumed_massert(KbT>0,"your MD engine does not pass the temperature to plumed, you must specify it using TEMP");
  beta_=1./KbT;

//set umbrellas
  parse("SIGMA",sigma_);
  double min_cv;
  double max_cv;
  parse("MIN_CV",min_cv);
  parse("MAX_CV",max_cv);
  plumed_massert(min_cv<max_cv,"MIN_CV should be smaller than MAX_CV");
  tot_umbrellas_=1+std::round((max_cv-min_cv)/sigma_);
  center_.resize(tot_umbrellas_);
  for(unsigned i=0; i<tot_umbrellas_; i++)
    center_[i]=(max_cv-min_cv)/(tot_umbrellas_-1)*i+min_cv;

//set other stuff
  parse("PACE",stride_);
  parse("OBSERVATION_STEPS",obs_steps_);
  plumed_massert(obs_steps_!=0,"minimum is OBSERVATION_STEPS=1");
  obs_cv_.resize(obs_steps_);
  barrier_=std::numeric_limits<double>::infinity();
  parse("BARRIER",barrier_);

//deltaFs file
  std::string deltaFsFileName;
  parse("FILE",deltaFsFileName);
  parse("PRINT_STRIDE",print_stride_);
  std::string fmt;
  parse("FMT",fmt);

  bool add_P0=false;
  parseFlag("ADD_P0",add_P0);
  if(add_P0)
    P0_contribution_=1;
  else
    P0_contribution_=0;

//work flag
  parseFlag("CALC_WORK",calc_work_);

//multiple walkers //TODO implement also external mw for cp2k
  unsigned walker_rank;
  bool walkers_mpi=false;
  parseFlag("WALKERS_MPI",walkers_mpi);
  if(walkers_mpi)
  {
    if(comm.Get_rank()==0)//multi_sim_comm works on first rank only
    {
      NumWalkers_=multi_sim_comm.Get_size();
      walker_rank=multi_sim_comm.Get_rank();
    }
    comm.Bcast(NumWalkers_,0); //if each walker has more than one processor update them all
    comm.Bcast(walker_rank,0);
  }
  else
  {
    NumWalkers_=1;
    walker_rank=0;
  }

//parallelization stuff
  NumParallel_=comm.Get_size();
  rank_=comm.Get_rank();
  bool serial=false;
  parseFlag("SERIAL",serial);
  if(serial)
  {
    log.printf(" -- SERIAL: running without loop parallelization\n");
    NumParallel_=1;
    rank_=0;
  }

  checkRead();

//printing some info
  log.printf("  Beta = %g\n",beta_);
  log.printf("  Running at TEMP = %g\n",1./(Kb*beta_));
  log.printf("  Updating the bias with PACE = %u\n",stride_);
  log.printf("  Total number of umbrellas = %u\n",tot_umbrellas_);
  log.printf("    with SIGMA = %g\n",sigma_);
  log.printf("    in CV range [%g,%g]\n",center_[0],center_[tot_umbrellas_-1]);
  log.printf("  Initial unbiased observation done for OBSERVATION_STEPS = %u\n",obs_steps_);
  if(barrier_!=std::numeric_limits<double>::infinity())
    log.printf("  guess for free energy BARRIER = %g\n",barrier_);
  if(add_P0)
    log.printf(" --- ADD_P0: adding P0 to the target\n");
  if(walkers_mpi)
    log.printf("  WALKERS_MPI: multiple walkers will share the same bias via mpi\n");
  if(NumWalkers_>1)
  {
    log.printf("  Using multiple walkers\n");
    log.printf("    number of walkers: %u\n",NumWalkers_);
    log.printf("    walker rank: %u\n",walker_rank);
  }
  if(NumParallel_>1)
    log.printf("  Using multiple threads per simulation: %u\n",NumParallel_);

//restart if needed
  if(getRestart())
  {
    IFile ifile;
    ifile.link(*this);
    if(ifile.FileExist(deltaFsFileName))
    {
      log.printf("  Restarting from: %s\n",deltaFsFileName.c_str());
      log.printf(" +++ make sure all simulation options are consistent! +++\n");
      ifile.open(deltaFsFileName);
      std::vector<std::string> deltaFsName;
      ifile.scanFieldList(deltaFsName);
      const unsigned pre_f=2; //time and rct
      const unsigned post_f=1; //print_stride
      plumed_massert(deltaFsName.size()-pre_f-post_f>=2,"RESTART - fewer than expected FIELDS found in '"+deltaFsFileName+"' file");
      if(tot_umbrellas_==0) //should never happen
        tot_umbrellas_=deltaFsName.size()-pre_f-post_f;
      plumed_massert(tot_umbrellas_==deltaFsName.size()-pre_f-post_f,"RESTART - something wrong with '"+deltaFsFileName+"' file: steps not matching size");
      const std::string firstW=deltaFsName[pre_f];
      const std::string lastW=deltaFsName[deltaFsName.size()-post_f-1];
      const std::size_t _pos=firstW.find_first_of("_");
      std::size_t test_pos=firstW.find_last_of("_");
      plumed_massert(_pos==test_pos,"RESTART - something wrong with '"+deltaFsFileName+"' file: more than 1 depencency (_) found");
      test_pos=lastW.find_last_of("_");
      plumed_massert(_pos==test_pos,"RESTART - something wrong with '"+deltaFsFileName+"' file: wrong format of FIELDS names");
      const std::string read_min_cv=firstW.substr(_pos+1,firstW.size()-_pos-1);
      const std::string used_min_cv=std::to_string(center_[0]);
      plumed_massert(used_min_cv==read_min_cv,"mismatch between provided MIN_CV and the one in restart");
      const std::string read_max_cv=lastW.substr(_pos+1,lastW.size()-_pos-1);
      const std::string used_max_cv=std::to_string(center_[center_.size()-1]);
      plumed_massert(used_max_cv==read_max_cv,"mismatch between provided MAX_CV and the one in restart");
    //initialize
      deltaF_.resize(tot_umbrellas_);
      obs_steps_=0; //avoid initializing again
    //read steps from file
      int restart_stride;
      ifile.scanField("print_stride",restart_stride);
      plumed_massert(restart_stride==(int)print_stride_,"also PRINT_STRIDE must be consistent to avoid problems with multiple restarts");
      ifile.allowIgnoredFields(); //this allows for multiple restart, but without checking for consistency between them!
      double time;
      while(ifile.scanField("time",time)) //room for improvements: only last line is important
      {
        if(calc_work_)
          old_deltaF_=deltaF_;
        ifile.scanField("rct",rct_);
        for(unsigned i=0; i<tot_umbrellas_; i++)
          ifile.scanField(deltaFsName[pre_f+i],deltaF_[i]);
        ifile.scanField();
        counter_++;
      }
      log.printf("  Successfully read %d lines, up to t=%g\n",counter_,time);
      counter_=(1+(counter_-1)*print_stride_)*NumWalkers_; //adjust counter
      if(NumWalkers_>1)
      {
        my_rct_=rct_; //better than 0
        log.printf(" +++ WARNING +++ the single walker rct estimate is resetted to the global rct\n");
      }
      ifile.reset(false);
      ifile.close();
    }
    else
      log.printf(" +++ WARNING +++ restart requested, but no '%s' file found!\n",deltaFsFileName.c_str());
  }
//sync all walkers to avoid opening files before reding is over (see also METAD)
  comm.Barrier();
  if(comm.Get_rank()==0 && walkers_mpi)
    multi_sim_comm.Barrier();

//setup deltaFs file, without opening it
  deltaFsOfile_.link(*this);
  if(NumWalkers_>1)
  {
    if(walker_rank>0)
      deltaFsFileName="/dev/null"; //only first walker writes on file
    deltaFsOfile_.enforceSuffix("");
  }
  deltaFsOfile_.open(deltaFsFileName);
  if(fmt.length()>0)
    deltaFsOfile_.fmtField(" "+fmt);
  deltaFsOfile_.setHeavyFlush(); //do I need it?
  deltaFsOfile_.addConstantField("print_stride");
  deltaFsOfile_.printField("print_stride",(int)print_stride_);

//add output components
  if(NumWalkers_>1)
  {
    addComponent("rct");
    componentIsNotPeriodic("rct");
    getPntrToComponent("rct")->set(my_rct_);
  }
  if(calc_work_)
  {
    addComponent("work");
    componentIsNotPeriodic("work");
  }

//Bibliography
  log.printf("  Bibliography: ");
  log<<plumed.cite("P. Piaggi and M. Parrinello, Phys. Rev. Lett. 122 (5), 050601 (2019)");
  log<<plumed.cite("M. Invernizzi and M. Parrinello, J. Phys. Chem. Lett. 11, 2731-2736 (2020)");
  log.printf("\n");
}

void OPESumbrellas_test::calculate()
{
  if(obs_steps_>0) //no bias before initialization
    return;

  const double cv=getArgument(0);
  long double sum=0;
  long double der_sum_cv=0;
  for(unsigned i=rank_; i<tot_umbrellas_; i+=NumParallel_)
  {
    const double dist_i=difference(0,center_[i],cv)/sigma_; //this function takes care of PBC if present
    const long double add_i=std::exp(static_cast<long double>(-0.5*std::pow(dist_i,2)+beta_*deltaF_[i]));
    sum+=add_i;
    der_sum_cv-=dist_i/sigma_*add_i;
  }
  if(NumParallel_>1)
  {
    comm.Sum(sum);
    comm.Sum(der_sum_cv);
  }
  sum+=P0_contribution_;
  const unsigned tot_steps=tot_umbrellas_+P0_contribution_;

  current_bias_=-1./beta_*std::log(sum/tot_steps);
  setBias(current_bias_);
  setOutputForce(0,1./beta_*der_sum_cv/sum);

//calculate work
  if(calc_work_)
  {
    long double old_sum=0;
    for(unsigned i=rank_; i<tot_umbrellas_; i+=NumParallel_)
      old_sum+=std::exp(static_cast<long double>(-0.5*std::pow(difference(0,center_[i],cv)/sigma_,2)+beta_*old_deltaF_[i]));
    if(NumParallel_>1)
      comm.Sum(old_sum);
    old_sum+=P0_contribution_;
    work_+=-1./beta_*std::log(sum/old_sum);
  }

  afterCalculate_=true;
}

void OPESumbrellas_test::update()
{
  if(getStep()%stride_!=0)
    return;
  if(isFirstStep_) //skip very first step, as in METAD
  {
    isFirstStep_=false;
    if(obs_steps_!=1) //if obs_steps_==1 go on with initialization
      return;
  }
  if(obs_steps_>0)
  {
    obs_cv_[counter_]=getArgument(0);
    counter_++;
    if(counter_==obs_steps_)
    {
      init_from_obs();
      counter_=NumWalkers_; //all preliminary observations count 1
      obs_steps_=0; //no more observation
    }
    return;
  }
  plumed_massert(afterCalculate_,"OPESumbrellas_test::update() must be called after OPESumbrellas_test::calculate() to work properly");
  afterCalculate_=false;

//work done by the bias in one iteration
  if(calc_work_)
  {
    getPntrToComponent("work")->set(work_);
    work_=0;
    old_deltaF_=deltaF_;
  }

//update averages. this assumes that calculate() always runs before update(), thus uses current_bias_
  const double cv=getArgument(0);
  if(NumWalkers_==1)
  {
    counter_++;
    const double increment=1./beta_*std::log1p(std::exp(static_cast<long double>(beta_*(current_bias_-rct_)))/(counter_-1.));
    for(unsigned i=0; i<tot_umbrellas_; i++)
      deltaF_[i]+=increment-1./beta_*std::log1p(std::exp(static_cast<long double>(-0.5*std::pow(difference(0,center_[i],cv)/sigma_,2)+beta_*(current_bias_-rct_+deltaF_[i])))/(counter_-1.));
    rct_+=increment+1./beta_*std::log1p(-1./counter_);
  }
  else
  {
    std::vector<double> all_bias(NumWalkers_);
    std::vector<double> all_cv(NumWalkers_);
    if(comm.Get_rank()==0)
    {
      multi_sim_comm.Allgather(current_bias_,all_bias);
      multi_sim_comm.Allgather(cv,all_cv);
    }
    comm.Bcast(all_bias,0);
    comm.Bcast(all_cv,0);
    for(unsigned w=0; w<NumWalkers_; w++)
    {
      counter_++;
      const double increment_w=1./beta_*std::log1p(std::exp(static_cast<long double>(beta_*(all_bias[w]-rct_)))/(counter_-1.));
      for(unsigned i=0; i<tot_umbrellas_; i++)
        deltaF_[i]+=increment_w-1./beta_*std::log1p(std::exp(static_cast<long double>(-0.5*std::pow(difference(0,center_[i],all_cv[w])/sigma_,2)+beta_*(all_bias[w]-rct_+deltaF_[i])))/(counter_-1.));
      rct_+=increment_w+1./beta_*std::log1p(-1./counter_);
    }
    //calc single walker rct
    const unsigned single_counter=counter_/NumWalkers_;
    const double increment=1./beta_*std::log1p(std::exp(static_cast<long double>(beta_*(current_bias_-my_rct_)))/(single_counter-1.));
    my_rct_+=increment+1./beta_*std::log1p(-1./single_counter);
    getPntrToComponent("rct")->set(my_rct_);
  }

//write to file
  if((counter_/NumWalkers_-1)%print_stride_==0)
  {
    deltaFsOfile_.printField("time",getTime());
    deltaFsOfile_.printField("rct",rct_);
    for(unsigned i=0; i<tot_umbrellas_; i++)
      deltaFsOfile_.printField("deltaF_"+std::to_string(center_[i]),deltaF_[i]);
    deltaFsOfile_.printField();
  }
}

void OPESumbrellas_test::init_from_obs()
{
//in case of multiple walkers gather all the statistics
  if(NumWalkers_>1)
  {
    obs_steps_*=NumWalkers_;
    std::vector<double> all_obs_cv(obs_steps_);
    if(comm.Get_rank()==0)
      multi_sim_comm.Allgather(obs_cv_,all_obs_cv);
    comm.Bcast(all_obs_cv,0);
    obs_cv_=all_obs_cv;
  }

//initialize deltaF_
  deltaF_.resize(tot_umbrellas_);
  for(unsigned i=0; i<tot_umbrellas_; i++)
  {
    deltaF_[i]=0.5*std::pow(difference(0,center_[i],obs_cv_[0])/sigma_,2)/beta_;
    if(deltaF_[i]>barrier_)
      deltaF_[i]=barrier_;
  }
  for(unsigned i=0; i<tot_umbrellas_; i++)
    for(unsigned t=1; t<obs_steps_; t++) //starts from t=1
      deltaF_[i]+=-1./beta_*std::log1p(std::exp(static_cast<long double>(-0.5*std::pow(difference(0,center_[i],obs_cv_[t])/sigma_,2)+beta_*deltaF_[i]))/t)-1./beta_*std::log1p(-1./(1.+t));

  if(calc_work_)
    old_deltaF_=deltaF_;
  obs_cv_.clear();

//print initialization to file
  deltaFsOfile_.printField("time",getTime());
  deltaFsOfile_.printField("rct",rct_);
  for(unsigned i=0; i<tot_umbrellas_; i++)
    deltaFsOfile_.printField("deltaF_"+std::to_string(center_[i]),deltaF_[i]);
  deltaFsOfile_.printField();
}

}
}
