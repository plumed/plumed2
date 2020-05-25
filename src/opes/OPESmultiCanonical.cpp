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

//+PLUMEDOC BIAS OPES_MULTICANONICAL
/*
\par Examples

OPES_MULTICANONICAL ...
  LABEL=opes
  ARG=ene
  PACE=50
  TEMP=300
  MIN_TEMP=260
  MAX_TEMP=350
... OPES_MULTICANONICAL


*/
//+ENDPLUMEDOC

class OPESmultiCanonical : public bias::Bias {

private:
  bool isFirstStep_;
  bool afterCalculate_;
  unsigned NumParallel_;
  unsigned rank_;
  unsigned NumWalkers_;
  unsigned long counter_;

  unsigned stride_;
  unsigned obs_steps_;
  std::vector<double> obs_ene_;
  unsigned steps_beta_;

  double beta_;
  std::vector<double> beta_p_;
  std::vector<double> deltaF_;
  double rct_;
  double my_rct_;
  double current_bias_;

  double border_weight_;
  double tot_steps_;

  bool calc_work_;
  double work_;
  std::vector<double> old_deltaF_;

  OFile deltaFsOfile_;
  unsigned print_stride_;

public:
  OPESmultiCanonical(const ActionOptions&);
  void calculate() override;
  void update() override;
  void init_integration_grid();
  void init_from_obs();
  unsigned estimate_steps(const double,const double,const std::vector<double>&,const std::string) const;
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(OPESmultiCanonical,"OPES_MULTICANONICAL")

void OPESmultiCanonical::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","PACE","10","how often the bias is updated");
  keys.add("compulsory","OBSERVATION_STEPS","100","number of unbiased initial steps to collect statistics."
                        " When using STEPS_TEMP it is forced to 1");
//tempertature stuff
  keys.add("compulsory","TEMP","-1","temperature. If not specified tries to get it from MD engine");
  keys.add("compulsory","MIN_TEMP","the minimum of the temperature range");
  keys.add("compulsory","MAX_TEMP","the maximum of the temperature range");
  keys.add("optional","STEPS_TEMP","manually set the number of intermediate temperatures");
//deltaFs file
  keys.add("compulsory","FILE","DELTAFS","a file with the estimate of the relative \\f$\\Delta F\\f$ for each component of the target");
  keys.add("optional","PRINT_STRIDE","( default=100 ) stride for printing to DELTAFS file");
  keys.add("optional","FMT","specify format for DELTAFS file");
//miscellaneous
  keys.add("optional","BORDER_WEIGHT","set it greater than 1 to obtain better sampling of the max and min thermodynamics conditions");
  keys.addFlag("CALC_WORK",false,"calculate the work done by the bias between each update");
  keys.addFlag("WALKERS_MPI",false,"switch on MPI version of multiple walkers");
  keys.addFlag("SERIAL",false,"perform calculations in serial");
  keys.use("RESTART");

//output components
  componentsAreNotOptional(keys);
  keys.addOutputComponent("rct","WALKERS_MPI","single walker estimate of \\f$c(t)\\f$: \\f$\\frac{1}{\\beta}\\log \\lange e^{\\beta V} \\rangle\\f$, all walkers should converge to the same value");
  keys.addOutputComponent("work","CALC_WORK","work done by the bias between each update"); //calculating this maybe is only a useless overhead...
}

OPESmultiCanonical::OPESmultiCanonical(const ActionOptions&ao)
  : PLUMED_BIAS_INIT(ao)
  , isFirstStep_(true)
  , afterCalculate_(false)
  , counter_(0)
  , rct_(0)
  , my_rct_(0)
  , work_(0)
  , print_stride_(100)
{
  //TODO is there a way to check that ARG is actually energy?
  plumed_massert(getNumberOfArguments()==1,"olny ENERGY should be given as ARG");

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
  temp=1./(Kb*beta_);

//set temp range
  double min_temp;
  double max_temp;
  parse("MIN_TEMP",min_temp);
  parse("MAX_TEMP",max_temp);
  plumed_massert(max_temp>=min_temp,"MAX_TEMP should be bigger than MIN_TEMP");
  const double tol=1e-3; //if temp is taken from MD engine it might be numerically slightly different
  if(temp<(1-tol)*min_temp || temp>(1+tol)*max_temp)
    log.printf(" +++ WARNING +++ running at TEMP=%g which is outside the chosen temperature range\n",temp);
  beta_p_.resize(2);
  beta_p_[0]=1./(Kb*max_temp); //min_beta
  beta_p_[1]=1./(Kb*min_temp); //max_beta

//set other stuff
  parse("PACE",stride_);
  parse("OBSERVATION_STEPS",obs_steps_);
  plumed_massert(obs_steps_!=0,"minimum is OBSERVATION_STEPS=1");
  steps_beta_=0;
  parse("STEPS_TEMP",steps_beta_);
  if(steps_beta_!=0)
    obs_steps_=1;
  obs_ene_.resize(obs_steps_);

  border_weight_=1;
  parse("BORDER_WEIGHT",border_weight_);

//deltaFs file
  std::string deltaFsFileName;
  parse("FILE",deltaFsFileName);
  parse("PRINT_STRIDE",print_stride_);
  std::string fmt;
  parse("FMT",fmt);

//work flag
  parseFlag("CALC_WORK",calc_work_);

//multiple walkers //external mw for cp2k not supported, but anyway cp2k cannot put bias on energy!
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
  log.printf("  Initial unbiased observation done for OBSERVATION_STEPS = %u\n",obs_steps_);
  if(steps_beta_!=0)
    log.printf(" +++ WARNING +++ STEPS_TEMP is used, thus OBSERVATION_STEP is set to 1\n"
               "                 Custom initial deltaFs can be set with a fake RESTART\n");
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
    //check range consistency
      const unsigned pre_f=2; //time and rct
      const unsigned post_f=1; //print_stride
      plumed_massert(deltaFsName.size()-pre_f-post_f>=2,"RESTART - fewer than expected FIELDS found in '"+deltaFsFileName+"' file");
      if(steps_beta_==0)
        steps_beta_=deltaFsName.size()-pre_f-post_f;
      plumed_massert(steps_beta_==deltaFsName.size()-pre_f-post_f,"RESTART - something wrong with '"+deltaFsFileName+"' file: steps not matching size");
      const std::string firstW=deltaFsName[pre_f];
      const std::string lastW=deltaFsName[deltaFsName.size()-post_f-1];
      const std::size_t _pos=firstW.find_first_of("_");
      std::size_t test_pos=firstW.find_last_of("_");
      plumed_massert(_pos==test_pos,"RESTART - something wrong with '"+deltaFsFileName+"' file: more than 1 depencency (_) found");
      test_pos=lastW.find_last_of("_");
      plumed_massert(_pos==test_pos,"RESTART - something wrong with '"+deltaFsFileName+"' file: wrong format of FIELDS names");
      const std::string read_min_beta=firstW.substr(_pos+1,firstW.size()-_pos-1);
      const std::string used_min_beta=std::to_string(beta_p_[0]);
      plumed_massert(used_min_beta==read_min_beta,"mismatch between provided MAX_TEMP and the one in restart");
      const std::string read_max_beta=lastW.substr(_pos+1,lastW.size()-_pos-1);
      const std::string used_max_beta=std::to_string(beta_p_[1]);
      plumed_massert(used_max_beta==read_max_beta,"mismatch between provided MIN_TEMP and the one in restart");
    //initialize
      init_integration_grid();
      deltaF_.resize(steps_beta_);
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
        for(unsigned i=0; i<steps_beta_; i++)
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

void OPESmultiCanonical::calculate()
{
  if(obs_steps_>0) //no bias before initialization
    return;

  const double ene=getArgument(0);
  long double sum=0;
  long double der_sum_ene=0;
  for(unsigned i=rank_; i<steps_beta_; i+=NumParallel_)
  {
    long double add_i=std::exp(static_cast<long double>((beta_-beta_p_[i])*ene+beta_*deltaF_[i]));
    if(i==0 || i==steps_beta_-1)
      add_i*=border_weight_;
    sum+=add_i;
    der_sum_ene+=(beta_-beta_p_[i])*add_i;
  }
  if(NumParallel_>1)
  {
    comm.Sum(sum);
    comm.Sum(der_sum_ene);
  }

  current_bias_=-1./beta_*std::log(sum/tot_steps_);
  setBias(current_bias_);
  setOutputForce(0,1./beta_*der_sum_ene/sum);

//calculate work
  if(calc_work_)
  {
    long double old_sum=0;
    for(unsigned i=rank_; i<steps_beta_; i+=NumParallel_)
    {
      long double add_i=std::exp(static_cast<long double>((beta_-beta_p_[i])*ene+beta_*old_deltaF_[i]));
      if(i==0 || i==steps_beta_-1)
        add_i*=border_weight_;
      old_sum+=add_i;
    }
    if(NumParallel_>1)
      comm.Sum(old_sum);
    work_+=-1./beta_*std::log(sum/old_sum);
  }

  afterCalculate_=true;
}

void OPESmultiCanonical::update()
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
    obs_ene_[counter_]=getArgument(0);
    counter_++;
    if(counter_==obs_steps_)
    {
      log.printf("\nAction OPES_MULTICANONICAL\n");
      init_from_obs();
      log.printf("Finished initialization\n\n");
      counter_=NumWalkers_; //all preliminary observations count 1
      obs_steps_=0; //no more observation
    }
    return;
  }
  plumed_massert(afterCalculate_,"OPESmultiCanonical::update() must be called after OPESmultiCanonical::calculate() to work properly");
  afterCalculate_=false;

//work done by the bias in one iteration
  if(calc_work_)
  {
    getPntrToComponent("work")->set(work_);
    work_=0;
    old_deltaF_=deltaF_;
  }

//update averages. This assumes that calculate() always runs before update(), thus uses current_bias_
  const double ene=getArgument(0);
  if(NumWalkers_==1)
  {
    counter_++;
    const double increment=1./beta_*std::log1p(std::exp(static_cast<long double>(beta_*(current_bias_-rct_)))/(counter_-1.));
    for(unsigned i=0; i<steps_beta_; i++)
      deltaF_[i]+=increment-1./beta_*std::log1p(std::exp(static_cast<long double>((beta_-beta_p_[i])*ene+beta_*(current_bias_-rct_+deltaF_[i])))/(counter_-1.));
    rct_+=increment+1./beta_*std::log1p(-1./counter_);
  }
  else
  {
    std::vector<double> all_bias(NumWalkers_);
    std::vector<double> all_ene(NumWalkers_);
    if(comm.Get_rank()==0)
    {
      multi_sim_comm.Allgather(current_bias_,all_bias);
      multi_sim_comm.Allgather(ene,all_ene);
    }
    comm.Bcast(all_bias,0);
    comm.Bcast(all_ene,0);
    for(unsigned w=0; w<NumWalkers_; w++)
    {
      counter_++;
      const double increment_w=1./beta_*std::log1p(std::exp(static_cast<long double>(beta_*(all_bias[w]-rct_)))/(counter_-1.));
      for(unsigned i=0; i<steps_beta_; i++)
        deltaF_[i]+=increment_w-1./beta_*std::log1p(std::exp(static_cast<long double>((beta_-beta_p_[i])*all_ene[w]+beta_*(all_bias[w]-rct_+deltaF_[i])))/(counter_-1.));
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
    for(unsigned i=0; i<steps_beta_; i++)
      deltaFsOfile_.printField("deltaF_"+std::to_string(beta_p_[i]),deltaF_[i]);
    deltaFsOfile_.printField();
  }
}

void OPESmultiCanonical::init_integration_grid()
{
  plumed_massert(steps_beta_>0,"zero steps given for integration grid");
//initialize beta grid
  const double min_beta=beta_p_[0];
  const double max_beta=beta_p_[1];
  beta_p_.clear();
  beta_p_.resize(steps_beta_);
  if(steps_beta_==1)
    beta_p_[0]=(min_beta+max_beta)/2.;
  else
    for(unsigned i=0; i<steps_beta_; i++)
      beta_p_[i]=min_beta+i*(max_beta-min_beta)/(steps_beta_-1);

//initialize tot_steps_, depending on total number of steps and border_weight
  tot_steps_=steps_beta_+(border_weight_-1)*2;

//print some info
  log.printf("  Total temp steps (in beta) = %d\n",steps_beta_);
  for(unsigned i=0; i<steps_beta_; i++)
    log.printf("   %2d. beta=% .10f  temp=% g\n",i+1,beta_p_[i],1./(plumed.getAtoms().getKBoltzmann()*beta_p_[i]));
  if(NumParallel_>steps_beta_)
    log.printf(" +++ WARNING +++ number of parallel threads is greater than beta steps. Using SERIAL might be faster\n");
  if(border_weight_!=1)
    log.printf(" --- using a weight different from 1 for the border, BORDER_WEIGHT = %g ---\n",border_weight_);
}

void OPESmultiCanonical::init_from_obs()
{
//in case of multiple walkers gather all the statistics
  if(NumWalkers_>1)
  {
    obs_steps_*=NumWalkers_;
    std::vector<double> all_obs_ene(obs_steps_);
    if(comm.Get_rank()==0)
      multi_sim_comm.Allgather(obs_ene_,all_obs_ene);
    comm.Bcast(all_obs_ene,0);
    obs_ene_=all_obs_ene;
  }

//estimate steps from Neff and initialize integration grid
  if(steps_beta_==0) //beta_p_[1]=max_beta->min_temp and beta_p_[0]=min_beta->max_temp
    steps_beta_=estimate_steps(beta_p_[1]-beta_,beta_p_[0]-beta_,obs_ene_,"TEMP");
  init_integration_grid();

//initialize deltaF_
  deltaF_.resize(steps_beta_);
  for(unsigned i=0; i<steps_beta_; i++)
    deltaF_[i]=(beta_p_[i]/beta_-1.)*obs_ene_[0];
  for(unsigned i=0; i<steps_beta_; i++)
    for(unsigned t=1; t<obs_steps_; t++) //starts from t=1
      deltaF_[i]+=-1./beta_*std::log1p(std::exp(static_cast<long double>((beta_-beta_p_[i])*obs_ene_[t]+beta_*deltaF_[i]))/t)-1./beta_*std::log1p(-1./(1.+t));
  if(calc_work_)
    old_deltaF_=deltaF_;
  obs_ene_.clear();

//print initialization to file
  deltaFsOfile_.printField("time",getTime());
  deltaFsOfile_.printField("rct",rct_);
  for(unsigned i=0; i<steps_beta_; i++)
    deltaFsOfile_.printField("deltaF_"+std::to_string(beta_p_[i]),deltaF_[i]);
  deltaFsOfile_.printField();
}

unsigned OPESmultiCanonical::estimate_steps(const double left_side,const double right_side,const std::vector<double>& obs,const std::string msg) const
{ //uses Neff to estimate the grid spacing
  if(left_side==0 && right_side==0)
  {
    log.printf(" +++ WARNING +++ MIN_%s and MAX_%s are equal to %s, using single step\n",msg.c_str(),msg.c_str(),msg.c_str());
    return 1;
  }
  auto get_neff_HWHM=[](const double side,const std::vector<double>& obs,const double av_obs) //HWHM = half width at half maximum. neff is in general not symmetric
  {
  //func: Neff/N-0.5 is a function between -0.5 and 0.5
    auto func=[](const long double delta,const std::vector<double> obs, const double av_obs)
    {
      long double sum_w=0;
      long double sum_w2=0;
      for(unsigned t=0; t<obs.size(); t++)
      {
        const long double w=std::exp(-delta*(obs[t]-av_obs));
        sum_w+=w;
        sum_w2+=w*w;
      }
      return sum_w*sum_w/sum_w2/obs.size()-0.5;
    };
  //here we find the root of func using the regula falsi (false position) method
  //but any method would be OK, not much precision is needed. src/tools/RootFindingBase.h looked complicated
    const double tolerance=1e-4; //seems to be a good default
    double a=0; //default is right side case
    double func_a=0.5;
    double b=side;
    double func_b=func(side,obs,av_obs);
    if(func_b>=0)
      return 0.0; //no zero is present!
    if(b<0)//left side case
    {
      std::swap(a,b);
      std::swap(func_a,func_b);
    }
    double c=a;
    double func_c=func_a;
    while(std::abs(func_c)>tolerance)
    {
      if(func_a*func_c>0)
      {
        a=c;
        func_a=func_c;
      }
      else
      {
        b=c;
        func_b=func_c;
      }
      c=(a*func_b-b*func_a)/(func_b-func_a);
      func_c=func(c,obs,av_obs);//func is evaluated only here, it might be expensive
    }
    return std::abs(c);
  };

//set average to zero, for numerical stability
  double av_obs=0;
  for(unsigned t=0; t<obs.size(); t++)
    av_obs+=obs[t];
  av_obs/=obs.size();

//estimation
  double left_HWHM=0;
  if(left_side!=0)
    left_HWHM=get_neff_HWHM(left_side,obs,av_obs);
  double right_HWHM=0;
  if(right_side!=0)
    right_HWHM=get_neff_HWHM(right_side,obs,av_obs);
  if(left_HWHM==0)
  {
    right_HWHM*=2;
    if(left_side==0)
      log.printf(" --- MIN_%s is equal to %s\n",msg.c_str(),msg.c_str());
    else
      log.printf(" +++ WARNING +++ MIN_%s is very close to %s\n",msg.c_str(),msg.c_str());
  }
  if(right_HWHM==0)
  {
    left_HWHM*=2;
    if(right_side==0)
      log.printf(" --- MAX_%s is equal to %s\n",msg.c_str(),msg.c_str());
    else
      log.printf(" +++ WARNING +++ MAX_%s is very close to %s\n",msg.c_str(),msg.c_str());
  }
  if(left_HWHM==0 && right_HWHM==0)
  {
    log.printf(" +++ WARNING +++ %s range is very narrow, using MIN_%s and MAX_%s as only steps\n",msg.c_str(),msg.c_str(),msg.c_str());
    return 2;
  }
  const double grid_spacing=left_HWHM+right_HWHM;
  log.printf("  Estimated %s spacing (with beta) = %g\n",msg.c_str(),grid_spacing);
  unsigned steps=std::ceil(std::abs(right_side-left_side)/grid_spacing);
  plumed_massert(steps>1,"something went wrong and estimated grid spacing for "+msg+" gives a step="+std::to_string(steps));
  return steps;
}

}
}
