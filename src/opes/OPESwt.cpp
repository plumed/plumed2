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

//+PLUMEDOC BIAS OPES_WT
/*
On-the-fly probability enhanced sampling (OPES) with well-tempered target distribution \cite Invernizzi2020opeswt .

Contrary to \ref METAD, OPES is not filling the basins, but rather tries wery quickly to get a coarse idea of the full Free Energy Surface (FES), and then slowly refines it.
It is very fast in exploring in the first phase, and then becomes extremely conservative and does not change significantly the shape of the deposited bias any more (quasi-static regime).
For this reason it is possible to use standard umbrella sampling reweighting (see \ref REWEIGHT_BIAS) to analyse the trajectory.
The estimated \f$c(t)\f$ is printed for reference only, since it should converge to a fixed value as the bias converges.
This \f$c(t)\f$ should NOT be used for reweighting.
Similarly, the \f$Z_n\f$ is printed and it should converge when no new region of the CV-space is explored.

Notice that if the employed CV is degenerate and maps different metastable basins onto the same CV-space region, then OPES will practically get stuck, instead of completely reshaping the bias (as \ref METAD would do, especially if BIASFACTOR is high).
This can be useful to diagnostic problems with your collective variable.
If you have no way to improve the set of CVs, you might consider using \ref OPES_WT_EXPLORE instead.

The parameter BARRIER should be set to be at least equal to the highest free energy barrier you wish to overcome.
If it is much lower than that, you will not cross the barrier, if it is much higher, you will be slower in converging.
If you know wich one is the most stable basin of your system you should start your simulation from there.

By default SIGMA is adaptive, estimated from the fluctuations over ADAPTIVE_SIGMA_STRIDE simulation steps (similar to \ref METAD ADAPTIVE=DIFF, but contrary to that, no artifacts will appear and the bias will converge to the correct one).
This might not be an optimal choice.

To use uniform flat target, explicitly set BIASFACTOR=inf (but should be needed only in very specific cases).

Restart can be done from a KERNELS file, but it might not be perfect (due to limited precision when printing numbers to file, or usage of adaptive SIGMA).
For a perfect restart you need to use STATE_RFILE to read a checkpoint with all the needed info.
To save such checkpoints, define a STATE_WFILE and choose how often to print them with STATE_WSTRIDE.
By default this file is overwritten, but you can instead append to it using the flag STORE_STATES.

\par Examples

The following is a minimal working example:

\plumedfile
opes: OPES_WT ARG=cv PACE=500 BARRIER=40
\endplumedfile

Another more articulated one:

\plumedfile
OPES_WT ...
  LABEL=opes
  FILE=Kernels.data
  TEMP=300
  ARG=cv1,cv2
  SIGMA=2.37,5.19
  PACE=500
  BARRIER=60
  BIASFACTOR=inf
  STATE_WFILE=State.data
  STATE_WSTRIDE=50000
  WALKERS_MPI
... OPES_WT
\endplumedfile


*/
//+ENDPLUMEDOC

class OPESwt : public bias::Bias {

private:
  static const bool explore_=false; //FIXME manual toggle between standard and explore mode
  bool isFirstStep_;
  bool afterCalculate_;
  unsigned NumParallel_;
  unsigned rank_;
  unsigned NumWalkers_;
  unsigned walker_rank_;
  unsigned ncv_;
  unsigned long counter_;

  double kbt_;
  double biasfactor_;
  double bias_prefactor_;
  unsigned stride_;
  std::vector<double> sigma0_;
  unsigned adaptive_sigma_stride_;
  unsigned long adaptive_counter_;
  std::vector<double> av_cv_;
  std::vector<double> av_M2_;
  bool fixed_sigma_;
  double epsilon_;
  double sum_weights_;
  double sum_weights2_;
  double current_bias_;

  bool no_Zed_;
  double Zed_;
  double KDEnorm_;

  double threshold2_;
  bool recursive_merge_;
//kernels for now are diagonal truncated Gaussians
  struct kernel
  {
    double height;
    std::vector<double> center;
    std::vector<double> sigma;

    inline void merge_me_with(const kernel & );
    kernel(double h, const std::vector<double> & c,const std::vector<double> & s):
      height(h),center(c),sigma(s) {}
  };
  double cutoff2_;
  double val_at_cutoff_;
  inline double evaluateKernel(const kernel&,const std::vector<double>&) const;
  inline double evaluateKernel(const kernel&,const std::vector<double>&,std::vector<double>&);
  std::vector<kernel> kernels_;
  OFile kernelsOfile_;

  bool calc_work_;
  double work_;
  double old_KDEnorm_;
  double old_Zed_;
  std::vector<kernel> delta_kernels_;

  OFile stateOfile_;
  int wStateStride_;
  bool storeOldStates_;

public:
  explicit OPESwt(const ActionOptions&);
  void calculate() override;
  void update() override;
  double getProbAndDerivatives(const std::vector<double>&,std::vector<double>&);
  void addKernel(const kernel&,const bool);
  void addKernel(const double,const std::vector<double>&,const std::vector<double>&,const bool);
  unsigned getMergeableKernel(const std::vector<double>&,const unsigned);
  void dumpStateToFile();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(OPESwt,"OPES_WT")

void OPESwt::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","TEMP","-1","temperature. If not set, it is taken from MD engine, but not all MD codes provide it");
  keys.add("compulsory","PACE","the frequency for kernel deposition");
  keys.add("compulsory","SIGMA","0","the initial widths of the kernels. If not set, adaptive sigma will be used");
  keys.add("compulsory","BARRIER","the free energy barrier to be overcome. It is used to set BIASFACTOR, EPSILON, and KERNEL_CUTOFF to reasonable values");
  keys.add("compulsory","COMPRESSION_THRESHOLD","1","merge kernels if closer than this threshold, in units of sigma");
//extra options
  keys.add("optional","BIASFACTOR","the \\f$\\gamma\\f$ bias factor used for the well-tempered target \\f$p(\\mathbf{s})\\f$."
           " Set to 'inf' for uniform flat target");
  keys.add("optional","EPSILON","the value of the regularization constant for the probability");
  keys.add("optional","KERNEL_CUTOFF","truncate kernels at this distance, in units of sigma");
  keys.add("optional","ADAPTIVE_SIGMA_STRIDE","number of steps for measuring adaptive sigma. Default is 10xPACE");
  keys.addFlag("NO_ZED",false,"do not normalize over the explored CV space, \\f$Z_n=1\\f$");
  keys.addFlag("FIXED_SIGMA",false,"do not decrease sigma as simulation goes on. Can be added in a RESTART, to keep down the number of compressed kernels");
  keys.addFlag("RECURSIVE_MERGE_OFF",false,"do not recursively attempt kernel merging when a new one is added."
               " Faster, but total number of compressed kernels might grow and slow down things");
//kernels and state files
  keys.add("compulsory","FILE","KERNELS","a file in which the list of deposited kernels is stored");
  keys.add("optional","FMT","specify format for KERNELS file");
  keys.add("optional","STATE_RFILE","read from this file the compressed kernels and all the info needed to RESTART the simulation");
  keys.add("optional","STATE_WFILE","write to this file the compressed kernels and all the info needed to RESTART the simulation");
  keys.add("optional","STATE_WSTRIDE","number of MD steps between writing the STATE_WFILE. Default is only on CPT events (but not all MD codes set them)");
  keys.addFlag("STORE_STATES",false,"append to STATE_WFILE instead of ovewriting it each time");
//miscellaneous
  keys.addFlag("CALC_WORK",false,"calculate the work done by the bias between each update");
  keys.addFlag("WALKERS_MPI",false,"switch on MPI version of multiple walkers");
  keys.addFlag("SERIAL",false,"perform calculations in serial. Might be faster when number of compressed kernels is small, e.g. if only one CV is used");
  keys.use("RESTART");

//output components
  componentsAreNotOptional(keys);
  keys.addOutputComponent("rct","default","estimate of \\f$c(t)\\f$: \\f$\\frac{1}{\\beta}\\log \\lange e^{\\beta V} \\rangle\\f$, should become flat as the simulation converges. Do NOT use for reweighitng");
  keys.addOutputComponent("zed","default","estimate of \\f$Z_n=\\int_\\Omega_n \\Tilde{P}_n(\\mathbf{s})\\, d\\mathbf{s}\\f$, should become flat as no new CV-space region is explored");
  keys.addOutputComponent("neff","default","effective sample size");
  keys.addOutputComponent("nker","default","total number of compressed kernels used to represent the bias");
  keys.addOutputComponent("work","CALC_WORK","work done by the last kernel deposited");
}

OPESwt::OPESwt(const ActionOptions& ao)
  : PLUMED_BIAS_INIT(ao)
  , isFirstStep_(true)
  , afterCalculate_(false)
  , counter_(1)
  , Zed_(1)
  , work_(0)
{
  ncv_=getNumberOfArguments();
//set kbt_
  const double Kb=plumed.getAtoms().getKBoltzmann();
  double temp=-1;
  parse("TEMP",temp);
  kbt_=Kb*temp;
  if(kbt_<0)
  {
    kbt_=plumed.getAtoms().getKbT();
    plumed_massert(kbt_>0,"your MD engine does not pass the temperature to plumed, you must specify it using TEMP");
  }

//other compulsory input
  parse("PACE",stride_);

  double barrier=0;
  parse("BARRIER",barrier);
  plumed_massert(barrier>=0,"the BARRIER should be greater than zero");

  biasfactor_=barrier/kbt_;
  std::string biasfactor_str;
  parse("BIASFACTOR",biasfactor_str);
  if(biasfactor_str=="inf" || biasfactor_str=="INF")
  {
    biasfactor_=std::numeric_limits<double>::infinity();
    bias_prefactor_=1;
  }
  else
  {
    parse("BIASFACTOR",biasfactor_);
    plumed_massert(biasfactor_>1,"BIASFACTOR must be greater than one (use 'inf' for uniform target)");
    bias_prefactor_=1-1./biasfactor_;
  }
  if(explore_)
  {
    plumed_massert(!std::isinf(biasfactor_),"BIASFACTOR=inf is not compatible with EXPLORE mode");
    bias_prefactor_=biasfactor_-1;
  }

  adaptive_sigma_stride_=0;
  parse("ADAPTIVE_SIGMA_STRIDE",adaptive_sigma_stride_);
  parseVector("SIGMA",sigma0_);
  if(sigma0_[0]==0 && sigma0_.size()==1)
  {
    plumed_massert(!std::isinf(biasfactor_),"BIASFACTOR=inf is not compatible with adaptive SIGMA");
    sigma0_.clear();
    adaptive_counter_=0;
    if(adaptive_sigma_stride_==0)
      adaptive_sigma_stride_=10*stride_; //this is arbitrary, chosen from few tests
    av_cv_.resize(ncv_,0);
    av_M2_.resize(ncv_,0);
    plumed_massert(adaptive_sigma_stride_>=stride_,"better to chose ADAPTIVE_SIGMA_STRIDE >= PACE");
  }
  else
  {
    plumed_massert(sigma0_.size()==ncv_,"number of SIGMA parameters does not match number of arguments");
    plumed_massert(adaptive_sigma_stride_==0,"if SIGMA is set then it cannot be adaptive, thus ADAPTIVE_SIGMA_STRIDE should not be set");
    if(explore_)
    {
      for(unsigned i=0; i<ncv_; i++)
        sigma0_[i]*=std::sqrt(biasfactor_); //the sigma of the target is broader F_t(s)=1/gamma*F(s)
    }
  }

  epsilon_=std::exp(-barrier/bias_prefactor_/kbt_);
  parse("EPSILON",epsilon_);
  plumed_massert(epsilon_>0,"you must choose a value for EPSILON greater than zero. Is your BARRIER too high?");
  sum_weights_=std::pow(epsilon_,bias_prefactor_); //to avoid NANs we start with counter_=1 and w0=exp(beta*V0)
  sum_weights2_=sum_weights_*sum_weights_;

  double cutoff=sqrt(2.*barrier/bias_prefactor_/kbt_);
  if(explore_)
    cutoff=sqrt(2.*barrier/kbt_); //otherwise it is too small
  parse("KERNEL_CUTOFF",cutoff);
  plumed_massert(cutoff>0,"you must choose a value for KERNEL_CUTOFF greater than zero");
  cutoff2_=cutoff*cutoff;
  val_at_cutoff_=std::exp(-0.5*cutoff2_);

  threshold2_=1;
  parse("COMPRESSION_THRESHOLD",threshold2_);
  threshold2_*=threshold2_;
  if(threshold2_!=0)
    plumed_massert(threshold2_>0 && threshold2_<cutoff2_,"COMPRESSION_THRESHOLD cannot be bigger than the KERNEL_CUTOFF");

//optional stuff
  no_Zed_=false;
  parseFlag("NO_ZED",no_Zed_);
  if(no_Zed_)
  {//this makes it more gentle in the initial phase
    sum_weights_=1;
    sum_weights2_=1;
  }
  fixed_sigma_=false;
  parseFlag("FIXED_SIGMA",fixed_sigma_);
  bool recursive_merge_off=false;
  parseFlag("RECURSIVE_MERGE_OFF",recursive_merge_off);
  recursive_merge_=!recursive_merge_off;
  parseFlag("CALC_WORK",calc_work_);

//kernels file
  std::string kernelsFileName;
  parse("FILE",kernelsFileName);
  std::string fmt;
  parse("FMT",fmt);

//output checkpoint of current state
  std::string restartFileName;
  parse("STATE_RFILE",restartFileName);
  std::string stateFileName;
  parse("STATE_WFILE",stateFileName);
  wStateStride_=0;
  parse("STATE_WSTRIDE",wStateStride_);
  storeOldStates_=false;
  parseFlag("STORE_STATES",storeOldStates_);
  if(wStateStride_!=0 || storeOldStates_)
    plumed_massert(stateFileName.length()>0,"filename for storing simulation status not specified, use STATE_WFILE");
  if(stateFileName.length()>0 && wStateStride_==0)
    wStateStride_=-1;//will print only on CPT events (checkpoints set by some MD engines, like gromacs)

//multiple walkers //TODO implement also external mw for cp2k
  bool walkers_mpi=false;
  parseFlag("WALKERS_MPI",walkers_mpi);
  if(walkers_mpi)
  {
    if(comm.Get_rank()==0)//multi_sim_comm works on first rank only
    {
      NumWalkers_=multi_sim_comm.Get_size();
      walker_rank_=multi_sim_comm.Get_rank();
    }
    comm.Bcast(NumWalkers_,0); //if each walker has more than one processor update them all
    comm.Bcast(walker_rank_,0);
  }
  else
  {
    NumWalkers_=1;
    walker_rank_=0;
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

//restart if needed
  if(getRestart())
  {
    bool stateRestart=true;
    if(restartFileName.length()==0)
    {
      stateRestart=false;
      restartFileName=kernelsFileName;
    }
    IFile ifile;
    ifile.link(*this);
    if(ifile.FileExist(restartFileName))
    {
      ifile.open(restartFileName);
      log.printf("  RESTART - make sure all used options are compatible\n");
        log.printf("    restarting from: %s\n",restartFileName.c_str());
      if(stateRestart)
        log.printf("    it should be a STATE file (not a KERNELS file)\n");
      else
        log.printf(" +++ WARNING +++ RESTART from KERNELS might be approximate, use STATE_WFILE and STATE_RFILE to restart from the exact state\n");
      std::string old_biasfactor_str;
      ifile.scanField("biasfactor",old_biasfactor_str);
      if(old_biasfactor_str=="inf" || old_biasfactor_str=="INF")
      {
        if(!std::isinf(biasfactor_))
          log.printf(" +++ WARNING +++ previous bias factor was inf while now it is %g\n",biasfactor_);
      }
      else
      {
        double old_biasfactor;
        ifile.scanField("biasfactor",old_biasfactor);
        if(std::abs(biasfactor_-old_biasfactor)>1e-6*biasfactor_)
          log.printf(" +++ WARNING +++ previous bias factor was %g while now it is %g. diff = %g\n",old_biasfactor,biasfactor_,biasfactor_-old_biasfactor);
      }
      double old_epsilon;
      ifile.scanField("epsilon",old_epsilon);
      if(std::abs(epsilon_-old_epsilon)>1e-6*epsilon_)
        log.printf(" +++ WARNING +++ previous epsilon was %g while now it is %g. diff = %g\n",old_epsilon,epsilon_,epsilon_-old_epsilon);
      double old_cutoff;
      ifile.scanField("kernel_cutoff",old_cutoff);
      if(std::abs(cutoff-old_cutoff)>1e-6*cutoff)
        log.printf(" +++ WARNING +++ previous kernel_cutoff was %g while now it is %g. diff = %g\n",old_cutoff,cutoff,cutoff-old_cutoff);
      double old_threshold;
      const double threshold=sqrt(threshold2_);
      ifile.scanField("compression_threshold",old_threshold);
      if(std::abs(threshold-old_threshold)>1e-6*threshold)
        log.printf(" +++ WARNING +++ previous compression_threshold was %g while now it is %g. diff = %g\n",old_threshold,threshold,threshold-old_threshold);
      if(stateRestart)
      {
        ifile.scanField("zed",Zed_);
        ifile.scanField("sum_weights",sum_weights_);
        ifile.scanField("sum_weights2",sum_weights2_);
        std::string str_counter;
        ifile.scanField("counter",str_counter); //scanField does not handle unsigned
        counter_=std::stoul(str_counter);
        if(sigma0_.size()==0)
        {
          ifile.scanField("adaptive_counter",str_counter); //scanField does not handle long unsigned
          adaptive_counter_=std::stoul(str_counter);
          if(NumWalkers_>1)
          {
            for(unsigned w=0; w<NumWalkers_; w++)
            {
              for(unsigned i=0; i<ncv_; i++)
              {
                double tmp1,tmp2;
                ifile.scanField("av_cv_"+getPntrToArgument(i)->getName()+"_"+std::to_string(w),tmp1);
                ifile.scanField("av_M2_"+getPntrToArgument(i)->getName()+"_"+std::to_string(w),tmp2);
                if(w==walker_rank_)
                {
                  av_cv_[i]=tmp1;
                  av_M2_[i]=tmp2;
                }
              }
            }
          }
          else
          {
            for(unsigned i=0; i<ncv_; i++)
            {
              ifile.scanField("av_cv_"+getPntrToArgument(i)->getName(),av_cv_[i]);
              ifile.scanField("av_M2_"+getPntrToArgument(i)->getName(),av_M2_[i]);
            }
          }
        }
      }
      for(unsigned i=0; i<ncv_; i++)
      {
        if(getPntrToArgument(i)->isPeriodic())
        {
          std::string arg_min,arg_max;
          getPntrToArgument(i)->getDomain(arg_min,arg_max);
          std::string file_min,file_max;
          ifile.scanField("min_"+getPntrToArgument(i)->getName(),file_min);
          ifile.scanField("max_"+getPntrToArgument(i)->getName(),file_max);
          plumed_massert(file_min==arg_min,"mismatch between restart and ARG periodicity");
          plumed_massert(file_max==arg_max,"mismatch between restart and ARG periodicity");
        }
      }
      if(stateRestart)
      {
        double time;
        while(ifile.scanField("time",time))
        {
          std::vector<double> center(ncv_);
          std::vector<double> sigma(ncv_);
          double height;
          for(unsigned i=0; i<ncv_; i++)
            ifile.scanField(getPntrToArgument(i)->getName(),center[i]);
          for(unsigned i=0; i<ncv_; i++)
            ifile.scanField("sigma_"+getPntrToArgument(i)->getName(),sigma[i]);
          ifile.scanField("height",height);
          ifile.scanField();
          kernels_.emplace_back(height,center,sigma);
        }
        log.printf("    A total of %d kernels where read\n",kernels_.size());
      }
      else
      {
        ifile.allowIgnoredFields(); //this allows for multiple restart, but without checking for consistency between them!
        double time;
        while(ifile.scanField("time",time))
        {
          std::vector<double> center(ncv_);
          std::vector<double> sigma(ncv_);
          double height;
          double logweight;
          for(unsigned i=0; i<ncv_; i++)
            ifile.scanField(getPntrToArgument(i)->getName(),center[i]);
          for(unsigned i=0; i<ncv_; i++)
            ifile.scanField("sigma_"+getPntrToArgument(i)->getName(),sigma[i]);
          ifile.scanField("height",height);
          ifile.scanField("logweight",logweight);
          ifile.scanField();
          addKernel(height,center,sigma,false);
          const double weight=std::exp(logweight);
          sum_weights_+=weight; //this sum is slightly inaccurate, because when printing some precision is lost
          sum_weights2_+=weight*weight;
          counter_++;
        }
        if(!no_Zed_)
        {
          double sum_uprob=0;
          for(unsigned k=rank_; k<kernels_.size(); k+=NumParallel_)
            for(unsigned kk=0; kk<kernels_.size(); kk++)
              sum_uprob+=evaluateKernel(kernels_[kk],kernels_[k].center);
          if(NumParallel_>1)
            comm.Sum(sum_uprob);
          Zed_=sum_uprob/KDEnorm_/kernels_.size();
        }
        log.printf("    A total of %d kernels where read, and compressed to %d\n",counter_,kernels_.size());
      }
      ifile.reset(false);
      ifile.close();
    }
    else
    { //same behaviour as METAD
      std::string not_found_msg="RESTART requested, but file '"+restartFileName+"' was not found!";
      if(stateRestart)
        error(not_found_msg);
      else
        log.printf(" +++ WARNING +++ %s\n",not_found_msg.c_str());
    }
  }
  else if(restartFileName.length()>0)
    log.printf(" +++ WARNING +++ the provided STATE_RFILE will be ignored, since RESTART was not requested\n");

//sync all walkers to avoid opening files before reading is over (see also METAD)
  comm.Barrier();
  if(comm.Get_rank()==0 && walkers_mpi)
    multi_sim_comm.Barrier();

//setup output kernels file
  kernelsOfile_.link(*this);
  if(NumWalkers_>1)
  {
    if(walker_rank_>0)
      kernelsFileName="/dev/null"; //only first walker writes on file
    kernelsOfile_.enforceSuffix("");
  }
  kernelsOfile_.open(kernelsFileName);
  if(fmt.length()>0)
    kernelsOfile_.fmtField(" "+fmt);
  kernelsOfile_.setHeavyFlush(); //do I need it?
  //define and set const fields
  kernelsOfile_.addConstantField("biasfactor");
  kernelsOfile_.addConstantField("epsilon");
  kernelsOfile_.addConstantField("kernel_cutoff");
  kernelsOfile_.addConstantField("compression_threshold");
  for(unsigned i=0; i<ncv_; i++)
    kernelsOfile_.setupPrintValue(getPntrToArgument(i));
  kernelsOfile_.printField("biasfactor",biasfactor_);
  kernelsOfile_.printField("epsilon",epsilon_);
  kernelsOfile_.printField("kernel_cutoff",sqrt(cutoff2_));
  kernelsOfile_.printField("compression_threshold",sqrt(threshold2_));

//open file for storing state
  if(wStateStride_!=0)
  {
    stateOfile_.link(*this);
    if(NumWalkers_>1)
    {
      if(walker_rank_>0)
        stateFileName="/dev/null"; //only first walker writes on file
      stateOfile_.enforceSuffix("");
    }
    stateOfile_.open(stateFileName);
    if(fmt.length()>0)
      stateOfile_.fmtField(" "+fmt);
  }

//set initial old values for work
  KDEnorm_=explore_?counter_:sum_weights_;
  old_KDEnorm_=KDEnorm_;
  old_Zed_=Zed_;

//add and set output components
  addComponent("rct");
  componentIsNotPeriodic("rct");
  getPntrToComponent("rct")->set(kbt_*std::log(sum_weights_/counter_));
  addComponent("zed");
  componentIsNotPeriodic("zed");
  getPntrToComponent("zed")->set(Zed_);
  addComponent("neff");
  componentIsNotPeriodic("neff");
  getPntrToComponent("neff")->set(std::pow(1+sum_weights_,2)/(1+sum_weights2_));
  addComponent("nker");
  componentIsNotPeriodic("nker");
  getPntrToComponent("nker")->set(kernels_.size());
  if(calc_work_)
  {
    addComponent("work");
    componentIsNotPeriodic("work");
  }

//printing some info
  log.printf("  temperature T = %g\n",kbt_/Kb);
  log.printf("  beta = %g\n",1./kbt_);
  log.printf("  depositing new kernels with PACE = %d\n",stride_);
  log.printf("  expected BARRIER is %g\n",barrier);
  log.printf("  using target distribution with BIASFACTOR gamma = %g\n",biasfactor_);
  if(std::isinf(biasfactor_))
    log.printf("    (thus a uniform flat target distribution, no well-tempering)\n");
  if(sigma0_.size()==0)
  {
    log.printf("  adaptive SIGMA will be used, with ADAPTIVE_SIGMA_STRIDE = %d\n",adaptive_sigma_stride_);
    log.printf("    thus the first n=ADAPTIVE_SIGMA_STRIDE/PACE steps will have no bias, n = %d\n",adaptive_sigma_stride_/stride_);
  }
  else
  {
    log.printf("  kernels have initial SIGMA = ");
    for(unsigned i=0; i<ncv_; i++)
      log.printf(" %g",sigma0_[i]);
    log.printf("\n");
  }
  if(fixed_sigma_)
    log.printf(" -- FIXED_SIGMA: sigma will not decrease as the simulation proceeds\n");
  log.printf("  kernels are truncated with KERNELS_CUTOFF = %g\n",cutoff);
  if(cutoff<3.5)
    log.printf(" +++ WARNING +++ probably kernels are truncated too much\n");
  log.printf("  the value at cutoff is = %g\n",val_at_cutoff_);
  log.printf("  regularization EPSILON = %g\n",epsilon_);
  if(val_at_cutoff_>epsilon_*(1+1e-6))
    log.printf(" +++ WARNING +++ the KERNEL_CUTOFF might be too small for the given EPSILON\n");
  log.printf("  kernels will be compressed when closer than COMPRESSION_THRESHOLD = %g\n",sqrt(threshold2_));
  if(threshold2_==0)
    log.printf(" +++ WARNING +++ kernels will never merge, expect slowdowns\n");
  if(!recursive_merge_)
    log.printf(" -- RECURSIVE_MERGE_OFF: only one merge for each new kernel will be attempted. This is faster only if total number of kernels does not grow too much\n");
  if(no_Zed_)
    log.printf(" -- NO_ZED: using fixed normalization factor = %g\n",Zed_);
  if(wStateStride_!=0 && walker_rank_==0)
    log.printf("  state checkpoints are written on file %s with stride %d\n",stateFileName.c_str(),wStateStride_);
  if(walkers_mpi)
    log.printf(" -- WALKERS_MPI: if present, multiple replicas will communicate\n");
  if(NumWalkers_>1)
  {
    log.printf("  using multiple walkers\n");
    log.printf("    number of walkers: %d\n",NumWalkers_);
    log.printf("    walker rank: %d\n",walker_rank_);
  }
  int mw_warning=0;
  if(!walkers_mpi && comm.Get_rank()==0 && multi_sim_comm.Get_size()>(int)NumWalkers_)
    mw_warning=1;
  comm.Bcast(mw_warning,0);
  if(mw_warning) //log.printf messes up with comm, so never use it without Bcast!
    log.printf(" +++ WARNING +++ multiple replicas will NOT communicate unless the flag WALKERS_MPI is used\n");
  if(NumParallel_>1)
    log.printf("  using multiple threads per simulation: %d\n",NumParallel_);
  log.printf(" Bibliography ");
  log<<plumed.cite("M. Invernizzi and M. Parrinello, J. Phys. Chem. Lett. 11, 2731-2736 (2020)");
  log.printf("\n");
  if(explore_)
    log.printf("We are running in EXPLORE mode!\n");
}

void OPESwt::calculate()
{
  std::vector<double> cv(ncv_);
  for(unsigned i=0; i<ncv_; i++)
    cv[i]=getArgument(i);

  std::vector<double> der_prob(ncv_,0);
  const double prob=getProbAndDerivatives(cv,der_prob);
  current_bias_=kbt_*bias_prefactor_*std::log(prob/Zed_+epsilon_);
  setBias(current_bias_);
  for(unsigned i=0; i<ncv_; i++)
    setOutputForce(i,der_prob[i]==0?0:-kbt_*bias_prefactor_/(prob/Zed_+epsilon_)*der_prob[i]/Zed_);

//calculate work
  if(calc_work_)
  {
    double tot_delta=0;
    for(unsigned d=0; d<delta_kernels_.size(); d++)
      tot_delta+=evaluateKernel(delta_kernels_[d],cv);
    const double old_prob=(prob*KDEnorm_-tot_delta)/old_KDEnorm_;
    work_+=current_bias_-kbt_*bias_prefactor_*std::log(old_prob/old_Zed_+epsilon_);
  }

  afterCalculate_=true;
}

void OPESwt::update()
{
  if(isFirstStep_)//same in MetaD, useful for restarts?
  {
    isFirstStep_=false;
    return;
  }

//dump state if requested
  if( (wStateStride_>0 && getStep()%wStateStride_==0) || (wStateStride_==-1 && getCPT()) )
    dumpStateToFile();

//update variance if adaptive sigma
  if(sigma0_.size()==0)
  {
    adaptive_counter_++;
    unsigned tau=adaptive_sigma_stride_;
    if(adaptive_counter_<adaptive_sigma_stride_)
      tau=adaptive_counter_;
    for(unsigned i=0; i<ncv_; i++)
    { //Welford's online algorithm for standard deviation
      const double cv_i=getArgument(i);
      const double diff=difference(i,av_cv_[i],cv_i);
      av_cv_[i]+=diff/tau; //exponentially decaying average
      av_M2_[i]+=diff*difference(i,av_cv_[i],cv_i);
    }
    if(adaptive_counter_<adaptive_sigma_stride_ && counter_==1) //counter_>1 if restarting
      return;  //do not apply bias before having measured sigma
  }

//do update
  if(getStep()%stride_!=0)
    return;
  plumed_massert(afterCalculate_,"OPESwt::update() must be called after OPESwt::calculate() to work properly");
  afterCalculate_=false; //if needed implementation can be changed to avoid this

//work done by the bias in one iteration uses as zero reference a point at inf, so that the work is always positive
  if(calc_work_)
  {
    const double min_shift=kbt_*bias_prefactor_*std::log(old_Zed_/Zed_*old_KDEnorm_/KDEnorm_);
    getPntrToComponent("work")->set(work_-stride_*min_shift);
    work_=0;
  }
  old_Zed_=Zed_;
  old_KDEnorm_=KDEnorm_;
  delta_kernels_.clear();
  unsigned old_nker=kernels_.size();

//get new kernel height
  double height=std::exp(current_bias_/kbt_); //this assumes that calculate() always runs before update()

//update sum_weights_ and neff
  double sum_heights=height;
  double sum_heights2=height*height;
  if(NumWalkers_>1)
  {
    if(comm.Get_rank()==0)
    {
      multi_sim_comm.Sum(sum_heights);
      multi_sim_comm.Sum(sum_heights2);
    }
    comm.Bcast(sum_heights,0);
    comm.Bcast(sum_heights2,0);
  }
  counter_+=NumWalkers_;
  sum_weights_+=sum_heights;
  sum_weights2_+=sum_heights2;
  const double neff=std::pow(1+sum_weights_,2)/(1+sum_weights2_);
  getPntrToComponent("rct")->set(kbt_*std::log(sum_weights_/counter_));
  getPntrToComponent("neff")->set(neff);
  if(explore_)
  {
    KDEnorm_=counter_;
  //in opes explore the kernel height=1, because it is not multiplied by the weight
    height=1;
  }
  else
    KDEnorm_=sum_weights_;

//if needed, rescale sigma and height
  std::vector<double> sigma=sigma0_;
  if(sigma0_.size()==0)
  {
    sigma.resize(ncv_);
    if(explore_)
    {
      for(unsigned i=0; i<ncv_; i++)
        sigma[i]=std::sqrt(av_M2_[i]/adaptive_counter_);
    }
    else
    {
      if(counter_==1+NumWalkers_)
      { //very first estimate is from unbiased, thus must be adjusted
        for(unsigned i=0; i<ncv_; i++)
          av_M2_[i]/=(1-bias_prefactor_);
      }
      for(unsigned i=0; i<ncv_; i++)
        sigma[i]=std::sqrt(av_M2_[i]/adaptive_counter_*(1-bias_prefactor_));
    }
  }
  if(!fixed_sigma_)
  {
    const double size=explore_?counter_:neff; //for EXPLORE neff is not relevant
    const double s_rescaling=std::pow(size*(ncv_+2.)/4.,-1./(4+ncv_));
    for(unsigned i=0; i<ncv_; i++)
      sigma[i]*=s_rescaling;
  //the height should be divided by sqrt(2*pi)*sigma,
  //but this overall factor would be canceled when dividing by Zed
  //thus we skip it altogether, but keep the s_rescaling
    height/=std::pow(s_rescaling,ncv_);
  }

//get new kernel center
  std::vector<double> center(ncv_);
  for(unsigned i=0; i<ncv_; i++)
    center[i]=getArgument(i);

//add new kernel(s)
  if(NumWalkers_>1)
  {
    std::vector<double> all_height(NumWalkers_,0.0);
    std::vector<double> all_center(NumWalkers_*ncv_,0.0);
    std::vector<double> all_sigma(NumWalkers_*ncv_,0.0);
    if(comm.Get_rank()==0)
    {
      multi_sim_comm.Allgather(height,all_height); //TODO heights should be communicated only once
      multi_sim_comm.Allgather(center,all_center);
      multi_sim_comm.Allgather(sigma,all_sigma);
    }
    comm.Bcast(all_height,0);
    comm.Bcast(all_center,0);
    comm.Bcast(all_sigma,0);
    for(unsigned w=0; w<NumWalkers_; w++)
    {
      std::vector<double> center_w(all_center.begin()+ncv_*w,all_center.begin()+ncv_*(w+1));
      std::vector<double> sigma_w(all_sigma.begin()+ncv_*w,all_sigma.begin()+ncv_*(w+1));
      addKernel(all_height[w],center_w,sigma_w,true);
    }
  }
  else
    addKernel(height,center,sigma,true);
  getPntrToComponent("nker")->set(kernels_.size());

  //update Zed_
  if(!no_Zed_)
  {
    double sum_uprob=0;
    const unsigned ks=kernels_.size();
    const unsigned ds=delta_kernels_.size();
    const bool few_kernels=(ks*ks<(3*ks*ds+2*ds*ds*NumParallel_+100)); //this seems reasonable, but is not rigorous...
    if(few_kernels) //really needed? Probably is almost always false
    {
      for(unsigned k=rank_; k<kernels_.size(); k+=NumParallel_)
        for(unsigned kk=0; kk<kernels_.size(); kk++)
          sum_uprob+=evaluateKernel(kernels_[kk],kernels_[k].center);
      if(NumParallel_>1)
        comm.Sum(sum_uprob);
    }
    else
    {
    // Here instead of redoing the full summation, we add only the changes, knowing that
    // uprob = old_uprob + delta_uprob
    // and we also need to consider that in the new sum there are some novel centers and some disappeared ones
      double delta_sum_uprob=0;
      for(unsigned k=rank_; k<kernels_.size(); k+=NumParallel_)
      {
        for(unsigned d=0; d<delta_kernels_.size(); d++)
        {
          const double sign=delta_kernels_[d].height<0?-1:1; //take away contribution from kernels that are gone, and add the one from new ones
          delta_sum_uprob+=evaluateKernel(delta_kernels_[d],kernels_[k].center)+sign*evaluateKernel(kernels_[k],delta_kernels_[d].center);
        }
      }
      if(NumParallel_>1)
        comm.Sum(delta_sum_uprob);
      for(unsigned d=0; d<delta_kernels_.size(); d++)
      {
        for(unsigned dd=0; dd<delta_kernels_.size(); dd++)
        { //now subtract the delta_uprob added before, but not needed
          const double sign=delta_kernels_[d].height<0?-1:1;
          delta_sum_uprob-=sign*evaluateKernel(delta_kernels_[dd],delta_kernels_[d].center);
        }
      }
      sum_uprob=Zed_*old_KDEnorm_*old_nker+delta_sum_uprob;
    }
    Zed_=sum_uprob/KDEnorm_/kernels_.size();
    getPntrToComponent("zed")->set(Zed_);
  }
}

double OPESwt::getProbAndDerivatives(const std::vector<double> &cv,std::vector<double> &der_prob)
{
  double prob=0.0;
  for(unsigned k=rank_; k<kernels_.size(); k+=NumParallel_) //TODO add neighbor list
    prob+=evaluateKernel(kernels_[k],cv,der_prob);
  if(NumParallel_>1)
  {
    comm.Sum(prob);
    comm.Sum(der_prob);
  }
  //normalize the estimate
  prob/=KDEnorm_;
  for(unsigned i=0; i<ncv_; i++)
    der_prob[i]/=KDEnorm_;

  return prob;
}

void OPESwt::addKernel(const kernel &new_kernel,const bool write_to_file)
{
  addKernel(new_kernel.height,new_kernel.center,new_kernel.sigma,write_to_file);
}

void OPESwt::addKernel(const double height,const std::vector<double>& center,const std::vector<double>& sigma,const bool write_to_file)
{
  bool no_match=true;
  if(threshold2_!=0)
  {
    unsigned taker_k=getMergeableKernel(center,kernels_.size());
    if(taker_k<kernels_.size())
    {
      no_match=false;
      delta_kernels_.emplace_back(-1*kernels_[taker_k].height,kernels_[taker_k].center,kernels_[taker_k].sigma);
      kernels_[taker_k].merge_me_with(kernel(height,center,sigma));
      delta_kernels_.push_back(kernels_[taker_k]);
      if(recursive_merge_) //the overhead is worth it if it keeps low the total number of kernels
      {
      //TODO this second check could run only through the kernels closer than, say, 2*threshold
      //     the function getMergeableKernel could return a list of such neighbors
        unsigned giver_k=taker_k;
        taker_k=getMergeableKernel(kernels_[giver_k].center,giver_k);
        while(taker_k<kernels_.size())
        {
          delta_kernels_.pop_back();
          delta_kernels_.emplace_back(-1*kernels_[taker_k].height,kernels_[taker_k].center,kernels_[taker_k].sigma);
          if(taker_k>giver_k) //saves time when erasing
            std::swap(taker_k,giver_k);
          kernels_[taker_k].merge_me_with(kernels_[giver_k]);
          delta_kernels_.push_back(kernels_[taker_k]);
          kernels_.erase(kernels_.begin()+giver_k);
          giver_k=taker_k;
          taker_k=getMergeableKernel(kernels_[giver_k].center,giver_k);
        }
      }
    }
  }
  if(no_match)
  {
    kernels_.emplace_back(height,center,sigma);
    delta_kernels_.emplace_back(height,center,sigma);
  }

//write to file
  if(write_to_file)
  {
    kernelsOfile_.printField("time",getTime());
    for(unsigned i=0; i<ncv_; i++)
      kernelsOfile_.printField(getPntrToArgument(i),center[i]);
    for(unsigned i=0; i<ncv_; i++)
      kernelsOfile_.printField("sigma_"+getPntrToArgument(i)->getName(),sigma[i]);
    kernelsOfile_.printField("height",height);
    kernelsOfile_.printField("logweight",current_bias_/kbt_);
    kernelsOfile_.printField();
  }
}

unsigned OPESwt::getMergeableKernel(const std::vector<double> &giver_center,const unsigned giver_k)
{ //returns kernels_.size() if no match is found
  unsigned min_k=kernels_.size();
  double min_dist2=threshold2_;
  for(unsigned k=rank_; k<kernels_.size(); k+=NumParallel_) //TODO add neighbor list!
  {
    if(k==giver_k) //a kernel should not be merged with itself
      continue;
    double dist2=0;
    for(unsigned i=0; i<ncv_; i++)
    { //TODO implement merging through the border for periodic CVs
      const double d=(kernels_[k].center[i]-giver_center[i])/kernels_[k].sigma[i];
      dist2+=d*d;
      if(dist2>=min_dist2)
        break;
    }
    if(dist2<min_dist2)
    {
      min_dist2=dist2;
      min_k=k;
    }
  }
  if(NumParallel_>1)
  {
    std::vector<double> all_min_dist2(NumParallel_);
    std::vector<unsigned> all_min_k(NumParallel_);
    comm.Allgather(min_dist2,all_min_dist2);
    comm.Allgather(min_k,all_min_k);
    const unsigned best=std::distance(std::begin(all_min_dist2),std::min_element(std::begin(all_min_dist2),std::end(all_min_dist2)));
    if(all_min_dist2[best]<threshold2_)
      min_k=all_min_k[best];
  }
  return min_k;
}

void OPESwt::dumpStateToFile()
{
  if(storeOldStates_)
    stateOfile_.clearFields();
  else if(walker_rank_==0)
    stateOfile_.rewind();

  stateOfile_.addConstantField("biasfactor");
  stateOfile_.addConstantField("epsilon");
  stateOfile_.addConstantField("kernel_cutoff");
  stateOfile_.addConstantField("compression_threshold");
  stateOfile_.addConstantField("zed");
  stateOfile_.addConstantField("sum_weights");
  stateOfile_.addConstantField("sum_weights2");
  stateOfile_.addConstantField("counter");
  if(sigma0_.size()==0)
  {
    stateOfile_.addConstantField("adaptive_counter");
    if(NumWalkers_>1)
    {
      for(unsigned w=0; w<NumWalkers_; w++)
      {
        for(unsigned i=0; i<ncv_; i++)
        {
          stateOfile_.addConstantField("av_cv_"+getPntrToArgument(i)->getName()+"_"+std::to_string(w));
          stateOfile_.addConstantField("av_M2_"+getPntrToArgument(i)->getName()+"_"+std::to_string(w));
        }
      }
    }
    else
    {
      for(unsigned i=0; i<ncv_; i++)
      {
        stateOfile_.addConstantField("av_cv_"+getPntrToArgument(i)->getName());
        stateOfile_.addConstantField("av_M2_"+getPntrToArgument(i)->getName());
      }
    }
  }
  for(unsigned i=0; i<ncv_; i++) //print periodicity of CVs
    stateOfile_.setupPrintValue(getPntrToArgument(i));
  stateOfile_.printField("biasfactor",biasfactor_);
  stateOfile_.printField("epsilon",epsilon_);
  stateOfile_.printField("kernel_cutoff",sqrt(cutoff2_));
  stateOfile_.printField("compression_threshold",sqrt(threshold2_));
  stateOfile_.printField("zed",Zed_);
  stateOfile_.printField("sum_weights",sum_weights_);
  stateOfile_.printField("sum_weights2",sum_weights2_);
  stateOfile_.printField("counter",std::to_string(counter_)); //printField does not handle unsigned
  if(sigma0_.size()==0)
  {
    stateOfile_.printField("adaptive_counter",std::to_string(adaptive_counter_)); //printField does not handle long unsigned
    if(NumWalkers_>1)
    {
      std::vector<double> all_av_cv(NumWalkers_*ncv_,0.0);
      std::vector<double> all_av_M2(NumWalkers_*ncv_,0.0);
      if(comm.Get_rank()==0)
      {
        multi_sim_comm.Allgather(av_cv_,all_av_cv);
        multi_sim_comm.Allgather(av_M2_,all_av_M2);
      }
      comm.Bcast(all_av_cv,0);
      comm.Bcast(all_av_M2,0);
      for(unsigned w=0; w<NumWalkers_; w++)
      {
        for(unsigned i=0; i<ncv_; i++)
        {
          stateOfile_.printField("av_cv_"+getPntrToArgument(i)->getName()+"_"+std::to_string(w),all_av_cv[w*ncv_+i]);
          stateOfile_.printField("av_M2_"+getPntrToArgument(i)->getName()+"_"+std::to_string(w),all_av_M2[w*ncv_+i]);
        }
      }
    }
    else
    {
      for(unsigned i=0; i<ncv_; i++)
      {
        stateOfile_.printField("av_cv_"+getPntrToArgument(i)->getName(),av_cv_[i]);
        stateOfile_.printField("av_M2_"+getPntrToArgument(i)->getName(),av_M2_[i]);
      }
    }
  }
  for(unsigned k=0; k<kernels_.size(); k++)
  {
    stateOfile_.printField("time",getTime()); //this is not very usefull
    for(unsigned i=0; i<ncv_; i++)
      stateOfile_.printField(getPntrToArgument(i),kernels_[k].center[i]);
    for(unsigned i=0; i<ncv_; i++)
      stateOfile_.printField("sigma_"+getPntrToArgument(i)->getName(),kernels_[k].sigma[i]);
    stateOfile_.printField("height",kernels_[k].height);
    stateOfile_.printField();
  }
  if(!storeOldStates_)
    stateOfile_.flush();
}

inline double OPESwt::evaluateKernel(const kernel& G,const std::vector<double>& x) const
{ //NB: cannot be a method of kernel class, because uses external variables (for cutoff)
  double norm2=0;
  for(unsigned i=0; i<ncv_; i++)
  {
    const double diff_i=difference(i,G.center[i],x[i])/G.sigma[i];
    norm2+=diff_i*diff_i;
    if(norm2>=cutoff2_)
      return 0;
  }
  return G.height*(std::exp(-0.5*norm2)-val_at_cutoff_);
}

inline double OPESwt::evaluateKernel(const kernel& G,const std::vector<double>& x, std::vector<double> & acc_der)
{ //NB: cannot be a method of kernel class, because uses external variables (for cutoff)
  double norm2=0;
  std::vector<double> diff(ncv_);
  for(unsigned i=0; i<ncv_; i++)
  {
    diff[i]=difference(i,G.center[i],x[i])/G.sigma[i];
    norm2+=diff[i]*diff[i];
    if(norm2>=cutoff2_)
      return 0;
  }
  const double val=G.height*(std::exp(-0.5*norm2)-val_at_cutoff_);
  for(unsigned i=0; i<ncv_; i++)
    acc_der[i]-=diff[i]/G.sigma[i]*val; //NB: we accumulate the derivative into der
  return val;
}

inline void OPESwt::kernel::merge_me_with(const kernel & other)
{
  const double h=height+other.height;
  for(unsigned i=0; i<center.size(); i++)
  {
    const double c_i=(height*center[i]+other.height*other.center[i])/h;
    const double s_my_part=height*(sigma[i]*sigma[i]+center[i]*center[i]);
    const double s_other_part=other.height*(other.sigma[i]*other.sigma[i]+other.center[i]*other.center[i]);
    const double s2_i=(s_my_part+s_other_part)/h-c_i*c_i;
    center[i]=c_i;
    sigma[i]=sqrt(s2_i);
  }
  height=h;
}

}
}
