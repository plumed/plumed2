/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2020 of Michele Invernizzi.

The opes module is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The opes module is distributed in the hope that it will be useful,
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
#include "tools/OpenMP.h"

namespace PLMD {
namespace opes {

//+PLUMEDOC OPES_BIAS OPES_METAD
/*
On-the-fly probability enhanced sampling (\ref OPES "OPES") with metadynamics-like target distribution \cite Invernizzi2020rethinking.

This OPES_METAD action samples target distributions defined via their marginal \f$p^{\text{tg}}(\mathbf{s})\f$ over some collective variables (CVs), \f$\mathbf{s}=\mathbf{s}(\mathbf{x})\f$.
By default OPES_METAD targets the well-tempered distribution, \f$p^{\text{WT}}(\mathbf{s})\propto [P(\mathbf{s})]^{1/\gamma}\f$, where \f$\gamma\f$ is known as BIASFACTOR.
Similarly to \ref METAD, OPES_METAD optimizes the bias on-the-fly, with a given PACE.
It does so by reweighting via kernel density estimation the unbiased distribution in the CV space, \f$P(\mathbf{s})\f$.
A compression algorithm is used to prevent the number of kernels from growing linearly with the simulation time.
The bias at step \f$n\f$ is
\f[
V_n(\mathbf{s}) = (1-1/\gamma)\frac{1}{\beta}\log\left(\frac{\tilde{P}_n(\mathbf{s})}{Z_n}+\epsilon\right)\, .
\f]
See Ref.\cite Invernizzi2020rethinking for a complete description of the method.

As an intuitive picture, rather than gradually filling the metastable basins, OPES_METAD quickly tries to get a coarse idea of the full free energy surface (FES), and then slowly refines its details.
It has a fast initial exploration phase, and then becomes extremely conservative and does not significantly change the shape of the deposited bias any more, reaching a regime of quasi-static bias.
For this reason, it is possible to use standard umbrella sampling reweighting (see \ref REWEIGHT_BIAS) to analyse the trajectory.
At <a href="https://github.com/invemichele/opes/tree/master/postprocessing">this link</a> you can find some python scripts that work in a similar way to \ref sum_hills, but the preferred way to obtain a FES with OPES is via reweighting.
The estimated \f$c(t)\f$ is printed for reference only, since it should converge to a fixed value as the bias converges.
This \f$c(t)\f$ should NOT be used for reweighting.
Similarly, the \f$Z_n\f$ factor is printed only for reference, and it should converge when no new region of the CV-space is explored.

Notice that OPES_METAD is more sensitive to degenerate CVs than \ref METAD.
If the employed CVs map different metastable basins onto the same CV-space region, then OPES_METAD will remain stuck rather than completely reshaping the bias.
This can be useful to diagnose problems with your collective variable.
If it is not possible to improve the set of CVs and remove this degeneracy, then you might instead consider to use \ref OPES_METAD_EXPLORE or \ref METAD.
In this way you will be able to obtain an estimate of the FES, but be aware that you most likely will not reach convergence and thus this estimate will be subjected to systematic errors (see e.g. Fig.3 in \cite Pietrucci2017review).
On the contrary, if your CVs are not degenerate but only suboptimal, you should converge faster by using OPES_METAD instead of \ref METAD \cite Invernizzi2020rethinking.

The parameter BARRIER should be set to be at least equal to the highest free energy barrier you wish to overcome.
If it is much lower than that, you will not cross the barrier, if it is much higher, convergence might take a little longer.
If the system has a basin that is clearly more stable than the others, it is better to start the simulation from there.

By default, the kernels SIGMA is adaptive, estimated from the fluctuations over ADAPTIVE_SIGMA_STRIDE simulation steps (similar to \ref METAD ADAPTIVE=DIFF, but contrary to that, no artifacts are introduced and the bias will converge to the correct one).
However, notice that depending on the system this might not be the optimal choice for SIGMA.

You can target a uniform flat distribution by explicitly setting BIASFACTOR=inf.
However, this should be useful only in very specific cases.

Restart can be done from a KERNELS file, but it might not be perfect (due to limited precision when printing kernels to file, or usage of adaptive SIGMA).
For an exact restart you must use STATE_RFILE to read a checkpoint with all the needed info.
To save such checkpoints, define a STATE_WFILE and choose how often to print them with STATE_WSTRIDE.
By default this file is overwritten, but you can instead append to it using the flag STORE_STATES.

Multiple walkers are supported only with MPI communication, via the keyword WALKERS_MPI.

\par Examples

The following is a minimal working example:

\plumedfile
cv: DISTANCE ATOMS=1,2
opes: OPES_METAD ARG=cv PACE=100 BARRIER=40
PRINT STRIDE=100 FILE=COLVAR ARG=*
\endplumedfile

Another more articulated one:

\plumedfile
phi: TORSION ATOMS=5,7,9,15
psi: TORSION ATOMS=7,9,15,17

opes: OPES_METAD ...
  FILE=Kernels.data
  TEMP=300
  ARG=phi,psi
  SIGMA=0.15,0.15
  PACE=500
  BARRIER=50
  BIASFACTOR=inf
  STATE_RFILE=Restart.data
  STATE_WFILE=State.data
  STATE_WSTRIDE=50000
  STORE_STATES
  WALKERS_MPI
...

PRINT FMT=%g STRIDE=500 FILE=Colvar.data ARG=phi,psi,opes.*
\endplumedfile


*/
//+ENDPLUMEDOC

template <class mode>
class OPESmetad : public bias::Bias {

private:
  bool isFirstStep_;
  bool afterCalculate_;
  unsigned NumOMP_;
  unsigned NumParallel_;
  unsigned rank_;
  unsigned NumWalkers_;
  unsigned walker_rank_;
  unsigned long counter_;
  std::size_t ncv_;

  double kbt_;
  double biasfactor_;
  double bias_prefactor_;
  unsigned stride_;
  std::vector<double> sigma0_;
  std::vector<double> sigma_min_;
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
//kernels are truncated diagonal Gaussians
  struct kernel
  {
    double height;
    std::vector<double> center;
    std::vector<double> sigma;
    kernel(double h, const std::vector<double>& c,const std::vector<double>& s):
      height(h),center(c),sigma(s) {}
  };
  double cutoff2_;
  double val_at_cutoff_;
  inline void mergeKernels(kernel&,const kernel&); //merge the second one into the first one
  inline double evaluateKernel(const kernel&,const std::vector<double>&) const;
  inline double evaluateKernel(const kernel&,const std::vector<double>&,std::vector<double>&,std::vector<double>&);
  std::vector<kernel> kernels_; //all compressed kernels
  OFile kernelsOfile_;
//neighbour list stuff
  bool nlist_;
  double nlist_param_[2];
  std::vector<unsigned> nlist_index_;
  std::vector<double> nlist_center_;
  std::vector<double> nlist_dev2_;
  unsigned nlist_steps_;
  bool nlist_update_;
  bool nlist_pace_reset_;

  bool calc_work_;
  double work_;
  double old_KDEnorm_;
  double old_Zed_;
  std::vector<kernel> delta_kernels_;

  OFile stateOfile_;
  int wStateStride_;
  bool storeOldStates_;

public:
  explicit OPESmetad(const ActionOptions&);
  void calculate() override;
  void update() override;
  static void registerKeywords(Keywords& keys);

  double getProbAndDerivatives(const std::vector<double>&,std::vector<double>&);
  void addKernel(const kernel&,const bool);
  void addKernel(const double,const std::vector<double>&,const std::vector<double>&,const bool);
  unsigned getMergeableKernel(const std::vector<double>&,const unsigned);
  void updateNlist(const std::vector<double>&);
  void dumpStateToFile();
};

struct convergence { static const bool explore=false; };
typedef OPESmetad<convergence> OPESmetad_c;
PLUMED_REGISTER_ACTION(OPESmetad_c,"OPES_METAD")

//OPES_METAD_EXPLORE is very similar from the point of view of the code,
//but conceptually it is better to make it a separate BIAS action

//+PLUMEDOC OPES_BIAS OPES_METAD_EXPLORE
/*
On-the-fly probability enhanced sampling (\ref OPES "OPES") with well-tempered target distribution, exploration mode \cite future_paper .

This OPES_METAD_EXPLORE action samples the well-tempered target distribution, that is defined via its marginal \f$p^{\text{WT}}(\mathbf{s})\propto [P(\mathbf{s})]^{1/\gamma}\f$ over some collective variables (CVs), \f$\mathbf{s}=\mathbf{s}(\mathbf{x})\f$.
While \ref OPES_METAD does so by estimating the unbiased distribution \f$P(\mathbf{s})\f$, OPES_METAD_EXPLORE instead estimates on-the-fly the target \f$p^{\text{WT}}(\mathbf{s})\f$ and uses it to define the bias.
The bias at step \f$n\f$ is
\f[
V_n(\mathbf{s}) = (\gamma-1)\frac{1}{\beta}\log\left(\frac{\tilde{P}^{\text{WT}}_n(\mathbf{s})}{Z_n}+\epsilon\right)\, .
\f]
See Ref.\cite future_paper for a complete description of the method.

Compared to \ref OPES_METAD, OPES_METAD_EXPLORE is more similar to \ref METAD, because it allows the bias to vary significantly, thus enhancing exploration.
This goes at the expenses of a possibly slower convergence of the reweight estimate.
It is useful to look around when you have no idea of the BARRIER, or if you want to quickly test the effectiveness of a new CV, and see if it is degenerate or not.

Like \ref OPES_METAD, also OPES_METAD_EXPLORE uses a kernel density estimation with an on-the-fly compression algorithm.
The only difference is that it does not perfom reweight, since it estimates the sampled distribution and not the unbiased one.

\par Examples

The following is a minimal working example:

\plumedfile
cv: DISTANCE ATOMS=1,2
opes: OPES_METAD_EXPLORE ARG=cv PACE=500 BARRIER=40
PRINT STRIDE=100 FILE=COLVAR ARG=cv,opes.*
\endplumedfile
*/
//+ENDPLUMEDOC

struct exploration { static const bool explore=true; };
typedef OPESmetad<exploration> OPESmetad_e;
// For some reason, this is not seen correctly by cppcheck
// cppcheck-suppress unknownMacro
PLUMED_REGISTER_ACTION(OPESmetad_e,"OPES_METAD_EXPLORE")

template <class mode>
void OPESmetad<mode>::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","TEMP","-1","temperature. If not set, it is taken from MD engine, but not all MD codes provide it");
  keys.add("compulsory","PACE","the frequency for kernel deposition");
  keys.add("compulsory","SIGMA","ADAPTIVE","the initial widths of the kernels. If not set, adaptive sigma will be used and the ADAPTIVE_SIGMA_STRIDE should be set");
  keys.add("compulsory","BARRIER","the free energy barrier to be overcome. It is used to set BIASFACTOR, EPSILON, and KERNEL_CUTOFF to reasonable values");
  keys.add("compulsory","COMPRESSION_THRESHOLD","1","merge kernels if closer than this threshold, in units of sigma");
//extra options
  keys.add("optional","ADAPTIVE_SIGMA_STRIDE","number of steps for measuring adaptive sigma. Default is 10xPACE");
  keys.add("optional","SIGMA_MIN","never reduce SIGMA below this value");
  std::string info_biasfactor("the \\f$\\gamma\\f$ bias factor used for the well-tempered target \\f$p(\\mathbf{s})\\f$. ");
  if(mode::explore)
    info_biasfactor+="Cannot be 'inf'";
  else
    info_biasfactor+="Set to 'inf' for uniform flat target";
  keys.add("optional","BIASFACTOR",info_biasfactor);
  keys.add("optional","EPSILON","the value of the regularization constant for the probability");
  keys.add("optional","KERNEL_CUTOFF","truncate kernels at this distance, in units of sigma");
  keys.add("optional","NLIST_PARAMETERS","( default=3.,0.5 ) the two cutoff parameters for the kernels neighbor list");
  keys.addFlag("NLIST",false,"use neighbor list for kernels summation, faster but experimental");
  keys.addFlag("NLIST_PACE_RESET",false,"force the reset of the neighbor list at each PACE. Can be useful with WALKERS_MPI");
  keys.addFlag("FIXED_SIGMA",false,"do not decrease sigma as simulation goes on. Can be added in a RESTART, to keep in check the number of compressed kernels");
  keys.addFlag("RECURSIVE_MERGE_OFF",false,"do not recursively attempt kernel merging when a new one is added");
  keys.addFlag("NO_ZED",false,"do not normalize over the explored CV space, \\f$Z_n=1\\f$");
//kernels and state files
  keys.add("compulsory","FILE","KERNELS","a file in which the list of all deposited kernels is stored");
  keys.add("optional","FMT","specify format for KERNELS file");
  keys.add("optional","STATE_RFILE","read from this file the compressed kernels and all the info needed to RESTART the simulation");
  keys.add("optional","STATE_WFILE","write to this file the compressed kernels and all the info needed to RESTART the simulation");
  keys.add("optional","STATE_WSTRIDE","number of MD steps between writing the STATE_WFILE. Default is only on CPT events (but not all MD codes set them)");
  keys.addFlag("STORE_STATES",false,"append to STATE_WFILE instead of ovewriting it each time");
//miscellaneous
  keys.addFlag("CALC_WORK",false,"calculate the work done by the bias between each update");
  keys.addFlag("WALKERS_MPI",false,"switch on MPI version of multiple walkers");
  keys.addFlag("SERIAL",false,"perform calculations in serial");
  keys.use("RESTART");
  keys.use("UPDATE_FROM");
  keys.use("UPDATE_UNTIL");

//output components
  keys.addOutputComponent("rct","default","estimate of \\f$c(t)\\f$: \\f$\\frac{1}{\\beta}\\log \\langle e^{\\beta V} \\rangle\\f$, should become flat as the simulation converges. Do NOT use for reweighting");
  keys.addOutputComponent("zed","default","estimate of \\f$Z_n\\f$, should become flat as no new CV-space region is explored");
  keys.addOutputComponent("neff","default","effective sample size");
  keys.addOutputComponent("nker","default","total number of compressed kernels used to represent the bias");
  keys.addOutputComponent("work","CALC_WORK","work done by the last kernel deposited");
  keys.addOutputComponent("nlker","NLIST","number of kernels in the neighbor list");
  keys.addOutputComponent("nlsteps","NLIST","number of steps from last neighbor list update");
}

template <class mode>
OPESmetad<mode>::OPESmetad(const ActionOptions& ao)
  : PLUMED_BIAS_INIT(ao)
  , isFirstStep_(true)
  , afterCalculate_(false)
  , counter_(1)
  , ncv_(getNumberOfArguments())
  , Zed_(1)
  , work_(0)
{
  std::string error_in_input1("Error in input in action "+getName()+" with label "+getLabel()+": the keyword ");
  std::string error_in_input2(" could not be read correctly");

//set kbt_
  const double Kb=plumed.getAtoms().getKBoltzmann();
  kbt_=plumed.getAtoms().getKbT();
  double temp=-1;
  parse("TEMP",temp);
  if(temp>0)
  {
    if(kbt_>0 && std::abs(kbt_-Kb*temp)>1e-4)
      log.printf(" +++ WARNING +++ using TEMP=%g while MD engine uses %g\n",temp,kbt_/Kb);
    kbt_=Kb*temp;
  }
  plumed_massert(kbt_>0,"your MD engine does not pass the temperature to plumed, you must specify it using TEMP");

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
    if(biasfactor_str.length()>0)
      plumed_massert(Tools::convert(biasfactor_str,biasfactor_),error_in_input1+"BIASFACTOR"+error_in_input2);
    plumed_massert(biasfactor_>1,"BIASFACTOR must be greater than one (use 'inf' for uniform target)");
    bias_prefactor_=1-1./biasfactor_;
  }
  if(mode::explore)
  {
    plumed_massert(!std::isinf(biasfactor_),"BIASFACTOR=inf is not compatible with EXPLORE mode");
    bias_prefactor_=biasfactor_-1;
  }

  adaptive_sigma_stride_=0;
  parse("ADAPTIVE_SIGMA_STRIDE",adaptive_sigma_stride_);
  std::vector<std::string> sigma_str;
  parseVector("SIGMA",sigma_str);
  double dummy;
  if(sigma_str.size()==1 && !Tools::convert(sigma_str[0],dummy))
  {
    plumed_massert(sigma_str[0]=="ADAPTIVE" || sigma_str[0]=="adaptive",error_in_input1+"SIGMA"+error_in_input2);
    plumed_massert(!std::isinf(biasfactor_),"BIASFACTOR=inf is not compatible with adaptive SIGMA");
    adaptive_counter_=0;
    if(adaptive_sigma_stride_==0)
      adaptive_sigma_stride_=10*stride_; //NB: this is arbitrary, chosen from few tests
    av_cv_.resize(ncv_,0);
    av_M2_.resize(ncv_,0);
    plumed_massert(adaptive_sigma_stride_>=stride_,"better to chose ADAPTIVE_SIGMA_STRIDE >= PACE");
  }
  else
  {
    plumed_massert(sigma_str.size()==ncv_,"number of SIGMA parameters does not match number of arguments");
    plumed_massert(adaptive_sigma_stride_==0,"if SIGMA is not ADAPTIVE you cannot set an ADAPTIVE_SIGMA_STRIDE");
    sigma0_.resize(ncv_);
    for(unsigned i=0; i<ncv_; i++)
    {
      plumed_massert(Tools::convert(sigma_str[i],sigma0_[i]),error_in_input1+"SIGMA"+error_in_input2);
      if(mode::explore)
        sigma0_[i]*=std::sqrt(biasfactor_); //the sigma of the target is broader F_t(s)=1/gamma*F(s)
    }
  }
  parseVector("SIGMA_MIN",sigma_min_);
  plumed_massert(sigma_min_.size()==0 || sigma_min_.size()==ncv_,"number of SIGMA_MIN does not match number of arguments");

  epsilon_=std::exp(-barrier/bias_prefactor_/kbt_);
  parse("EPSILON",epsilon_);
  plumed_massert(epsilon_>0,"you must choose a value for EPSILON greater than zero. Is your BARRIER too high?");
  sum_weights_=std::pow(epsilon_,bias_prefactor_); //to avoid NANs we start with counter_=1 and w0=exp(beta*V0)
  sum_weights2_=sum_weights_*sum_weights_;

  double cutoff=sqrt(2.*barrier/bias_prefactor_/kbt_);
  if(mode::explore)
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

//setup neighbor list
  nlist_=false;
  parseFlag("NLIST",nlist_);
  nlist_pace_reset_=false;
  parseFlag("NLIST_PACE_RESET",nlist_pace_reset_);
  if(nlist_pace_reset_)
    nlist_=true;
  std::vector<double> nlist_param;
  parseVector("NLIST_PARAMETERS",nlist_param);
  if(nlist_param.size()==0)
  {
    nlist_param_[0]=3.0;//*cutoff2_ -> max distance of neighbors
    nlist_param_[1]=0.5;//*nlist_dev2_[i] -> condition for rebuilding
  }
  else
  {
    nlist_=true;
    plumed_massert(nlist_param.size()==2,"two cutoff parameters are needed for the neighbor list");
    plumed_massert(nlist_param[0]>1.0,"the first of NLIST_PARAMETERS must be greater than 1. The smaller the first, the smaller should be the second as well");
    const double min_PARAM_1=(1.-1./std::sqrt(nlist_param[0]))+0.16;
    plumed_massert(nlist_param[1]>0,"the second of NLIST_PARAMETERS must be greater than 0");
    plumed_massert(nlist_param[1]<=min_PARAM_1,"the second of NLIST_PARAMETERS must be smaller to avoid systematic errors. Largest suggested value is: 1.16-1/sqrt(PARAM_0) = "+std::to_string(min_PARAM_1));
    nlist_param_[0]=nlist_param[0];
    nlist_param_[1]=nlist_param[1];
  }
  nlist_center_.resize(ncv_);
  nlist_dev2_.resize(ncv_,0.);
  nlist_steps_=0;
  nlist_update_=true;

//optional stuff
  no_Zed_=false;
  parseFlag("NO_ZED",no_Zed_);
  if(no_Zed_)
  { //this makes it more gentle in the initial phase
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
  NumOMP_=OpenMP::getNumThreads();
  NumParallel_=comm.Get_size();
  rank_=comm.Get_rank();
  bool serial=false;
  parseFlag("SERIAL",serial);
  if(serial)
  {
    NumOMP_=1;
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
      bool tmp_nlist=nlist_;
      nlist_=false; // NLIST is not needed while restarting
      ifile.open(restartFileName);
      log.printf("  RESTART - make sure all used options are compatible\n");
      log.printf("    restarting from: %s\n",restartFileName.c_str());
      std::string action_name=getName();
      if(stateRestart)
      {
        log.printf("    it should be a STATE file (not a KERNELS file)\n");
        action_name+="_state";
      }
      else
      {
        log.printf(" +++ WARNING +++ RESTART from KERNELS might be approximate, use STATE_WFILE and STATE_RFILE to restart from the exact state\n");
        action_name+="_kernels";
      }
      std::string old_action_name;
      ifile.scanField("action",old_action_name);
      plumed_massert(action_name==old_action_name,"RESTART - mismatch between old and new action name. Expected '"+action_name+"', but found '"+old_action_name+"'");
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
        ifile.scanField("counter",counter_);
        if(sigma0_.size()==0)
        {
          ifile.scanField("adaptive_counter",adaptive_counter_);
          if(NumWalkers_==1)
          {
            for(unsigned i=0; i<ncv_; i++)
            {
              ifile.scanField("av_cv_"+getPntrToArgument(i)->getName(),av_cv_[i]);
              ifile.scanField("av_M2_"+getPntrToArgument(i)->getName(),av_M2_[i]);
            }
          }
          else
          {
            for(unsigned w=0; w<NumWalkers_; w++)
              for(unsigned i=0; i<ncv_; i++)
              {
                double tmp1,tmp2;
                const std::string arg_iw=getPntrToArgument(i)->getName()+"_"+std::to_string(w);
                ifile.scanField("av_cv_"+arg_iw,tmp1);
                ifile.scanField("av_M2_"+arg_iw,tmp2);
                if(w==walker_rank_)
                {
                  av_cv_[i]=tmp1;
                  av_M2_[i]=tmp2;
                }
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
          plumed_massert(file_min==arg_min,"RESTART - mismatch between old and new ARG periodicity");
          plumed_massert(file_max==arg_max,"RESTART - mismatch between old and new ARG periodicity");
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
        log.printf("    a total of %lu kernels where read\n",kernels_.size());
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
        KDEnorm_=mode::explore?counter_:sum_weights_;
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
        log.printf("    a total of %lu kernels where read, and compressed to %lu\n",counter_,kernels_.size());
      }
      ifile.reset(false);
      ifile.close();
      nlist_=tmp_nlist;
    }
    else //same behaviour as METAD
      plumed_merror("RESTART requested, but file '"+restartFileName+"' was not found!\n  Set RESTART=NO or provide a restart file");
    if(NumWalkers_>1) //make sure that all walkers are doing the same thing
    {
      const unsigned kernels_size=kernels_.size();
      std::vector<unsigned> all_kernels_size(NumWalkers_);
      if(comm.Get_rank()==0)
        multi_sim_comm.Allgather(kernels_size,all_kernels_size);
      comm.Bcast(all_kernels_size,0);
      bool same_number_of_kernels=true;
      for(unsigned w=1; w<NumWalkers_; w++)
        if(all_kernels_size[0]!=all_kernels_size[w])
          same_number_of_kernels=false;
      plumed_massert(same_number_of_kernels,"RESTART - not all walkers are reading the same file!");
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
  kernelsOfile_.addConstantField("action");
  kernelsOfile_.addConstantField("biasfactor");
  kernelsOfile_.addConstantField("epsilon");
  kernelsOfile_.addConstantField("kernel_cutoff");
  kernelsOfile_.addConstantField("compression_threshold");
  for(unsigned i=0; i<ncv_; i++)
    kernelsOfile_.setupPrintValue(getPntrToArgument(i));
  kernelsOfile_.printField("action",getName()+"_kernels");
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

//set initial old values
  KDEnorm_=mode::explore?counter_:sum_weights_;
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
  if(nlist_)
  {
    addComponent("nlker");
    componentIsNotPeriodic("nlker");
    addComponent("nlsteps");
    componentIsNotPeriodic("nlsteps");
  }

//printing some info
  log.printf("  temperature = %g\n",kbt_/Kb);
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
  if(sigma_min_.size()>0)
  {
    log.printf("  kernels have a SIGMA_MIN = ");
    for(unsigned i=0; i<ncv_; i++)
      log.printf(" %g",sigma_min_[i]);
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
  if(nlist_)
    log.printf(" -- NLIST: using neighbor list for kernels, with parameters: %g,%g\n",nlist_param_[0],nlist_param_[1]);
  if(nlist_pace_reset_)
    log.printf(" -- NLIST_PACE_RESET: forcing the neighbor list to update every PACE\n");
  if(no_Zed_)
    log.printf(" -- NO_ZED: using fixed normalization factor = %g\n",Zed_);
  if(wStateStride_!=0 && walker_rank_==0)
    log.printf("  state checkpoints are written on file %s with stride %d\n",stateFileName.c_str(),wStateStride_);
  if(walkers_mpi)
    log.printf(" -- WALKERS_MPI: if multiple replicas are present, they will share the same bias via MPI\n");
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
  if(serial)
    log.printf(" -- SERIAL: running without loop parallelization\n");
  log.printf(" Bibliography ");
  log<<plumed.cite("M. Invernizzi and M. Parrinello, J. Phys. Chem. Lett. 11, 2731-2736 (2020)");
  log.printf("\n");
}

template <class mode>
void OPESmetad<mode>::calculate()
{
//get cv
  std::vector<double> cv(ncv_);
  for(unsigned i=0; i<ncv_; i++)
    cv[i]=getArgument(i);

//check neighbor list
  if(nlist_)
  {
    nlist_steps_++;
    if(getExchangeStep())
      nlist_update_=true;
    else
    {
      for(unsigned i=0; i<ncv_; i++)
      {
        const double diff_i=difference(i,cv[i],nlist_center_[i]);
        if(diff_i*diff_i>nlist_param_[1]*nlist_dev2_[i])
        {
          nlist_update_=true;
          break;
        }
      }
    }
    if(nlist_update_)
      updateNlist(cv);
  }

//set bias and forces
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

template <class mode>
void OPESmetad<mode>::update()
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
      const double diff_i=difference(i,av_cv_[i],cv_i);
      av_cv_[i]+=diff_i/tau; //exponentially decaying average
      av_M2_[i]+=diff_i*difference(i,av_cv_[i],cv_i);
    }
    if(adaptive_counter_<adaptive_sigma_stride_ && counter_==1) //counter_>1 if restarting
      return; //do not apply bias before having measured sigma
  }

//do update
  if(getStep()%stride_!=0)
    return;
  plumed_massert(afterCalculate_,"OPESmetad::update() must be called after OPESmetad::calculate() to work properly");
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
  if(mode::explore)
  {
    KDEnorm_=counter_;
    height=1; //plain KDE, bias reweight does not enter here
  }
  else
    KDEnorm_=sum_weights_;

//if needed, rescale sigma and height
  std::vector<double> sigma=sigma0_;
  if(sigma0_.size()==0)
  {
    sigma.resize(ncv_);
    if(mode::explore)
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
    const double size=mode::explore?counter_:neff; //for EXPLORE neff is not relevant
    const double s_rescaling=std::pow(size*(ncv_+2.)/4.,-1./(4+ncv_));
    if(sigma_min_.size()==0)
    {
      for(unsigned i=0; i<ncv_; i++)
        sigma[i]*=s_rescaling;
      //the height should be divided by sqrt(2*pi)*sigma,
      //but this overall factor would be canceled when dividing by Zed
      //thus we skip it altogether, but keep the s_rescaling
      height/=std::pow(s_rescaling,ncv_);
    }
    else
    {
      for(unsigned i=0; i<ncv_; i++)
      {
        const double s_rescaling_i=std::max(s_rescaling,sigma_min_[i]/sigma[i]);
        sigma[i]*=s_rescaling_i;
        height/=s_rescaling_i;
      }
    }
  }

//get new kernel center
  std::vector<double> center(ncv_);
  for(unsigned i=0; i<ncv_; i++)
    center[i]=getArgument(i);

//add new kernel(s)
  if(NumWalkers_==1)
    addKernel(height,center,sigma,true);
  else
  {
    std::vector<double> all_height(NumWalkers_,0.0);
    std::vector<double> all_center(NumWalkers_*ncv_,0.0);
    std::vector<double> all_sigma(NumWalkers_*ncv_,0.0);
    if(comm.Get_rank()==0)
    {
      multi_sim_comm.Allgather(height,all_height); //heights were communicated also before...
      multi_sim_comm.Allgather(center,all_center);
      multi_sim_comm.Allgather(sigma,all_sigma);
    }
    comm.Bcast(all_height,0);
    comm.Bcast(all_center,0);
    comm.Bcast(all_sigma,0);
    if(nlist_)
    { //gather all the nlist_index_, so merging can be done using it
      std::vector<int> all_nlist_size(NumWalkers_);
      if(comm.Get_rank()==0)
      {
        all_nlist_size[walker_rank_]=nlist_index_.size();
        multi_sim_comm.Sum(all_nlist_size);
      }
      comm.Bcast(all_nlist_size,0);
      unsigned tot_size=0;
      for(unsigned w=0; w<NumWalkers_; w++)
        tot_size+=all_nlist_size[w];
      if(tot_size>0)
      {
        std::vector<int> disp(NumWalkers_);
        for(unsigned w=0; w<NumWalkers_-1; w++)
          disp[w+1]=disp[w]+all_nlist_size[w];
        std::vector<unsigned> all_nlist_index(tot_size);
        if(comm.Get_rank()==0)
          multi_sim_comm.Allgatherv(nlist_index_,all_nlist_index,&all_nlist_size[0],&disp[0]);
        comm.Bcast(all_nlist_index,0);
        std::set<unsigned> nlist_index_set(all_nlist_index.begin(),all_nlist_index.end()); //remove duplicates and sort
        nlist_index_.assign(nlist_index_set.begin(),nlist_index_set.end());
      }
    }
    for(unsigned w=0; w<NumWalkers_; w++)
    {
      std::vector<double> center_w(all_center.begin()+ncv_*w,all_center.begin()+ncv_*(w+1));
      std::vector<double> sigma_w(all_sigma.begin()+ncv_*w,all_sigma.begin()+ncv_*(w+1));
      addKernel(all_height[w],center_w,sigma_w,true);
    }
  }
  getPntrToComponent("nker")->set(kernels_.size());
  if(nlist_)
  {
    getPntrToComponent("nlker")->set(nlist_index_.size());
    if(nlist_pace_reset_)
      nlist_update_=true;
  }

  //update Zed_
  if(!no_Zed_)
  {
    double sum_uprob=0;
    const unsigned ks=kernels_.size();
    const unsigned ds=delta_kernels_.size();
    const bool few_kernels=(ks*ks<(3*ks*ds+2*ds*ds*NumParallel_+100)); //this seems reasonable, but is not rigorous...
    if(few_kernels) //really needed? Probably is almost always false
    {
      #pragma omp parallel num_threads(NumOMP_)
      {
        #pragma omp for reduction(+:sum_uprob) nowait
        for(unsigned k=rank_; k<kernels_.size(); k+=NumParallel_)
          for(unsigned kk=0; kk<kernels_.size(); kk++)
            sum_uprob+=evaluateKernel(kernels_[kk],kernels_[k].center);
      }
      if(NumParallel_>1)
        comm.Sum(sum_uprob);
    }
    else
    {
      // Here instead of redoing the full summation, we add only the changes, knowing that
      // uprob = old_uprob + delta_uprob
      // and we also need to consider that in the new sum there are some novel centers and some disappeared ones
      double delta_sum_uprob=0;
      if(!nlist_)
      {
        #pragma omp parallel num_threads(NumOMP_)
        {
          #pragma omp for reduction(+:delta_sum_uprob) nowait
          for(unsigned k=rank_; k<kernels_.size(); k+=NumParallel_)
          {
            for(unsigned d=0; d<delta_kernels_.size(); d++)
            {
              const double sign=delta_kernels_[d].height<0?-1:1; //take away contribution from kernels that are gone, and add the one from new ones
              delta_sum_uprob+=evaluateKernel(delta_kernels_[d],kernels_[k].center)+sign*evaluateKernel(kernels_[k],delta_kernels_[d].center);
            }
          }
        }
      }
      else
      {
        #pragma omp parallel num_threads(NumOMP_)
        {
          #pragma omp for reduction(+:delta_sum_uprob) nowait
          for(unsigned nk=rank_; nk<nlist_index_.size(); nk+=NumParallel_)
          {
            const unsigned k=nlist_index_[nk];
            for(unsigned d=0; d<delta_kernels_.size(); d++)
            {
              const double sign=delta_kernels_[d].height<0?-1:1; //take away contribution from kernels that are gone, and add the one from new ones
              delta_sum_uprob+=evaluateKernel(delta_kernels_[d],kernels_[k].center)+sign*evaluateKernel(kernels_[k],delta_kernels_[d].center);
            }
          }
        }
      }
      if(NumParallel_>1)
        comm.Sum(delta_sum_uprob);
      #pragma omp parallel num_threads(NumOMP_)
      {
        #pragma omp for reduction(+:delta_sum_uprob)
        for(unsigned d=0; d<delta_kernels_.size(); d++)
        {
          for(unsigned dd=0; dd<delta_kernels_.size(); dd++)
          { //now subtract the delta_uprob added before, but not needed
            const double sign=delta_kernels_[d].height<0?-1:1;
            delta_sum_uprob-=sign*evaluateKernel(delta_kernels_[dd],delta_kernels_[d].center);
          }
        }
      }
      sum_uprob=Zed_*old_KDEnorm_*old_nker+delta_sum_uprob;
    }
    Zed_=sum_uprob/KDEnorm_/kernels_.size();
    getPntrToComponent("zed")->set(Zed_);
  }
}

template <class mode>
double OPESmetad<mode>::getProbAndDerivatives(const std::vector<double>& cv,std::vector<double>& der_prob)
{
  double prob=0.0;
  if(!nlist_)
  {
    if(NumOMP_==1 || (unsigned)kernels_.size()<2*NumOMP_*NumParallel_)
    {
      // for performances and thread safety
      std::vector<double> dist(ncv_);
      for(unsigned k=rank_; k<kernels_.size(); k+=NumParallel_)
        prob+=evaluateKernel(kernels_[k],cv,der_prob,dist);
    }
    else
    {
      #pragma omp parallel num_threads(NumOMP_)
      {
        std::vector<double> omp_deriv(der_prob.size(),0.);
        // for performances and thread safety
        std::vector<double> dist(ncv_);
        #pragma omp for reduction(+:prob) nowait
        for(unsigned k=rank_; k<kernels_.size(); k+=NumParallel_)
          prob+=evaluateKernel(kernels_[k],cv,omp_deriv,dist);
        #pragma omp critical
        for(unsigned i=0; i<ncv_; i++)
          der_prob[i]+=omp_deriv[i];
      }
    }
  }
  else
  {
    if(NumOMP_==1 || (unsigned)nlist_index_.size()<2*NumOMP_*NumParallel_)
    {
      // for performances and thread safety
      std::vector<double> dist(ncv_);
      for(unsigned nk=rank_; nk<nlist_index_.size(); nk+=NumParallel_)
        prob+=evaluateKernel(kernels_[nlist_index_[nk]],cv,der_prob,dist);
    }
    else
    {
      #pragma omp parallel num_threads(NumOMP_)
      {
        std::vector<double> omp_deriv(der_prob.size(),0.);
        // for performances and thread safety
        std::vector<double> dist(ncv_);
        #pragma omp for reduction(+:prob) nowait
        for(unsigned nk=rank_; nk<nlist_index_.size(); nk+=NumParallel_)
          prob+=evaluateKernel(kernels_[nlist_index_[nk]],cv,omp_deriv,dist);
        #pragma omp critical
        for(unsigned i=0; i<ncv_; i++)
          der_prob[i]+=omp_deriv[i];
      }
    }
  }
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

template <class mode>
void OPESmetad<mode>::addKernel(const kernel &new_kernel,const bool write_to_file)
{
  addKernel(new_kernel.height,new_kernel.center,new_kernel.sigma,write_to_file);
}

template <class mode>
void OPESmetad<mode>::addKernel(const double height,const std::vector<double>& center,const std::vector<double>& sigma,const bool write_to_file)
{
  bool no_match=true;
  if(threshold2_!=0)
  {
    unsigned taker_k=getMergeableKernel(center,kernels_.size());
    if(taker_k<kernels_.size())
    {
      no_match=false;
      delta_kernels_.emplace_back(-1*kernels_[taker_k].height,kernels_[taker_k].center,kernels_[taker_k].sigma);
      mergeKernels(kernels_[taker_k],kernel(height,center,sigma));
      delta_kernels_.push_back(kernels_[taker_k]);
      if(recursive_merge_) //the overhead is worth it if it keeps low the total number of kernels
      {
        unsigned giver_k=taker_k;
        taker_k=getMergeableKernel(kernels_[giver_k].center,giver_k);
        while(taker_k<kernels_.size())
        {
          delta_kernels_.pop_back();
          delta_kernels_.emplace_back(-1*kernels_[taker_k].height,kernels_[taker_k].center,kernels_[taker_k].sigma);
          if(taker_k>giver_k) //saves time when erasing
            std::swap(taker_k,giver_k);
          mergeKernels(kernels_[taker_k],kernels_[giver_k]);
          delta_kernels_.push_back(kernels_[taker_k]);
          kernels_.erase(kernels_.begin()+giver_k);
          if(nlist_)
          {
            unsigned giver_nk=0;
            bool found_giver=false;
            for(unsigned nk=0; nk<nlist_index_.size(); nk++)
            {
              if(found_giver)
                nlist_index_[nk]--; //all the indexes shift due to erase
              if(nlist_index_[nk]==giver_k)
              {
                giver_nk=nk;
                found_giver=true;
              }
            }
            plumed_dbg_massert(found_giver,"problem with merging and NLIST");
            nlist_index_.erase(nlist_index_.begin()+giver_nk);
          }
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
    if(nlist_)
      nlist_index_.push_back(kernels_.size()-1);
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

template <class mode>
unsigned OPESmetad<mode>::getMergeableKernel(const std::vector<double>& giver_center,const unsigned giver_k)
{ //returns kernels_.size() if no match is found
  unsigned min_k=kernels_.size();
  double min_norm2=threshold2_;
  if(!nlist_)
  {
    #pragma omp parallel num_threads(NumOMP_)
    {
      unsigned min_k_omp = min_k;
      double min_norm2_omp = threshold2_;
      #pragma omp for nowait
      for(unsigned k=rank_; k<kernels_.size(); k+=NumParallel_)
      {
        if(k==giver_k) //a kernel should not be merged with itself
          continue;
        double norm2=0;
        for(unsigned i=0; i<ncv_; i++)
        {
          const double dist_i=difference(i,giver_center[i],kernels_[k].center[i])/kernels_[k].sigma[i];
          norm2+=dist_i*dist_i;
          if(norm2>=min_norm2_omp)
            break;
        }
        if(norm2<min_norm2_omp)
        {
          min_norm2_omp=norm2;
          min_k_omp=k;
        }
      }
      #pragma omp critical
      {
        if(min_norm2_omp < min_norm2)
        {
          min_norm2 = min_norm2_omp;
          min_k = min_k_omp;
        }
      }
    }
  }
  else
  {
    #pragma omp parallel num_threads(NumOMP_)
    {
      unsigned min_k_omp = min_k;
      double min_norm2_omp = threshold2_;
      #pragma omp for nowait
      for(unsigned nk=rank_; nk<nlist_index_.size(); nk+=NumParallel_)
      {
        const unsigned k=nlist_index_[nk];
        if(k==giver_k) //a kernel should not be merged with itself
          continue;
        double norm2=0;
        for(unsigned i=0; i<ncv_; i++)
        {
          const double dist_i=difference(i,giver_center[i],kernels_[k].center[i])/kernels_[k].sigma[i];
          norm2+=dist_i*dist_i;
          if(norm2>=min_norm2)
            break;
        }
        if(norm2<min_norm2_omp)
        {
          min_norm2_omp=norm2;
          min_k_omp=k;
        }
      }
      #pragma omp critical
      {
        if(min_norm2_omp < min_norm2)
        {
          min_norm2 = min_norm2_omp;
          min_k = min_k_omp;
        }
      }
    }
  }
  if(NumParallel_>1)
  {
    std::vector<double> all_min_norm2(NumParallel_);
    std::vector<unsigned> all_min_k(NumParallel_);
    comm.Allgather(min_norm2,all_min_norm2);
    comm.Allgather(min_k,all_min_k);
    const unsigned best=std::distance(std::begin(all_min_norm2),std::min_element(std::begin(all_min_norm2),std::end(all_min_norm2)));
    min_k=all_min_k[best];
  }
  return min_k;
}

template <class mode>
void OPESmetad<mode>::updateNlist(const std::vector<double>& new_center)
{
  if(kernels_.size()==0) //no need to check for neighbors
    return;

  nlist_center_=new_center;
  nlist_index_.clear();
  //first we gather all the nlist_index
  if(NumOMP_==1 || (unsigned)kernels_.size()<2*NumOMP_*NumParallel_)
  {
    for(unsigned k=rank_; k<kernels_.size(); k+=NumParallel_)
    {
      double norm2_k=0;
      for(unsigned i=0; i<ncv_; i++)
      {
        const double dist_ik=difference(i,nlist_center_[i],kernels_[k].center[i])/kernels_[k].sigma[i];
        norm2_k+=dist_ik*dist_ik;
      }
      if(norm2_k<=nlist_param_[0]*cutoff2_)
        nlist_index_.push_back(k);
    }
  }
  else
  {
    #pragma omp parallel num_threads(NumOMP_)
    {
      std::vector<unsigned> private_nlist_index;
      #pragma omp for nowait
      for(unsigned k=rank_; k<kernels_.size(); k+=NumParallel_)
      {
        double norm2_k=0;
        for(unsigned i=0; i<ncv_; i++)
        {
          const double dist_ik=difference(i,nlist_center_[i],kernels_[k].center[i])/kernels_[k].sigma[i];
          norm2_k+=dist_ik*dist_ik;
        }
        if(norm2_k<=nlist_param_[0]*cutoff2_)
          private_nlist_index.push_back(k);
      }
      #pragma omp critical
      nlist_index_.insert(nlist_index_.end(),private_nlist_index.begin(),private_nlist_index.end());
    }
    if(recursive_merge_)
      std::sort(nlist_index_.begin(),nlist_index_.end());
  }
  if(NumParallel_>1)
  {
    std::vector<int> all_nlist_size(NumParallel_);
    all_nlist_size[rank_]=nlist_index_.size();
    comm.Sum(all_nlist_size);
    unsigned tot_size=0;
    for(unsigned r=0; r<NumParallel_; r++)
      tot_size+=all_nlist_size[r];
    if(tot_size>0)
    {
      std::vector<int> disp(NumParallel_);
      for(unsigned r=0; r<NumParallel_-1; r++)
        disp[r+1]=disp[r]+all_nlist_size[r];
      std::vector<unsigned> local_nlist_index=nlist_index_;
      nlist_index_.resize(tot_size);
      comm.Allgatherv(local_nlist_index,nlist_index_,&all_nlist_size[0],&disp[0]);
      if(recursive_merge_)
        std::sort(nlist_index_.begin(),nlist_index_.end());
    }
  }
  //calculate the square deviation
  std::vector<double> dev2(ncv_,0.);
  for(unsigned k=rank_; k<nlist_index_.size(); k+=NumParallel_)
  {
    for(unsigned i=0; i<ncv_; i++)
    {
      const double diff_ik=difference(i,nlist_center_[i],kernels_[nlist_index_[k]].center[i]);
      dev2[i]+=diff_ik*diff_ik;
    }
  }
  if(NumParallel_>1)
    comm.Sum(dev2);
  for(unsigned i=0; i<ncv_; i++)
  {
    if(dev2[i]==0) //e.g. if nlist_index_.size()==0
      nlist_dev2_[i]=std::pow(kernels_.back().sigma[i],2);
    else
      nlist_dev2_[i]=dev2[i]/nlist_index_.size();
  }
  getPntrToComponent("nlker")->set(nlist_index_.size());
  getPntrToComponent("nlsteps")->set(nlist_steps_);
  nlist_steps_=0;
  nlist_update_=false;
}

template <class mode>
void OPESmetad<mode>::dumpStateToFile()
{
//gather adaptive sigma info if needed
//doing this while writing to file can lead to misterious slowdowns
  std::vector<double> all_av_cv;
  std::vector<double> all_av_M2;
  if(sigma0_.size()==0 && NumWalkers_>1)
  {
    all_av_cv.resize(NumWalkers_*ncv_);
    all_av_M2.resize(NumWalkers_*ncv_);
    if(comm.Get_rank()==0)
    {
      multi_sim_comm.Allgather(av_cv_,all_av_cv);
      multi_sim_comm.Allgather(av_M2_,all_av_M2);
    }
    comm.Bcast(all_av_cv,0);
    comm.Bcast(all_av_M2,0);
  }

//rewrite header or rewind file
  if(storeOldStates_)
    stateOfile_.clearFields();
  else if(walker_rank_==0)
    stateOfile_.rewind();
//define fields
  stateOfile_.addConstantField("action");
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
    if(NumWalkers_==1)
    {
      for(unsigned i=0; i<ncv_; i++)
      {
        stateOfile_.addConstantField("av_cv_"+getPntrToArgument(i)->getName());
        stateOfile_.addConstantField("av_M2_"+getPntrToArgument(i)->getName());
      }
    }
    else
    {
      for(unsigned w=0; w<NumWalkers_; w++)
        for(unsigned i=0; i<ncv_; i++)
        {
          const std::string arg_iw=getPntrToArgument(i)->getName()+"_"+std::to_string(w);
          stateOfile_.addConstantField("av_cv_"+arg_iw);
          stateOfile_.addConstantField("av_M2_"+arg_iw);
        }
    }
  }
//print fields
  for(unsigned i=0; i<ncv_; i++) //periodicity of CVs
    stateOfile_.setupPrintValue(getPntrToArgument(i));
  stateOfile_.printField("action",getName()+"_state");
  stateOfile_.printField("biasfactor",biasfactor_);
  stateOfile_.printField("epsilon",epsilon_);
  stateOfile_.printField("kernel_cutoff",sqrt(cutoff2_));
  stateOfile_.printField("compression_threshold",sqrt(threshold2_));
  stateOfile_.printField("zed",Zed_);
  stateOfile_.printField("sum_weights",sum_weights_);
  stateOfile_.printField("sum_weights2",sum_weights2_);
  stateOfile_.printField("counter",counter_);
  if(sigma0_.size()==0)
  {
    stateOfile_.printField("adaptive_counter",adaptive_counter_);
    if(NumWalkers_==1)
    {
      for(unsigned i=0; i<ncv_; i++)
      {
        stateOfile_.printField("av_cv_"+getPntrToArgument(i)->getName(),av_cv_[i]);
        stateOfile_.printField("av_M2_"+getPntrToArgument(i)->getName(),av_M2_[i]);
      }
    }
    else
    {
      for(unsigned w=0; w<NumWalkers_; w++)
        for(unsigned i=0; i<ncv_; i++)
        {
          const std::string arg_iw=getPntrToArgument(i)->getName()+"_"+std::to_string(w);
          stateOfile_.printField("av_cv_"+arg_iw,all_av_cv[w*ncv_+i]);
          stateOfile_.printField("av_M2_"+arg_iw,all_av_M2[w*ncv_+i]);
        }
    }
  }
//print kernels
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
//make sure file is written even if small
  if(!storeOldStates_)
    stateOfile_.flush();
}

template <class mode>
inline double OPESmetad<mode>::evaluateKernel(const kernel& G,const std::vector<double>& x) const
{ //NB: cannot be a method of kernel class, because uses external variables (for cutoff)
  double norm2=0;
  for(unsigned i=0; i<ncv_; i++)
  {
    const double dist_i=difference(i,G.center[i],x[i])/G.sigma[i];
    norm2+=dist_i*dist_i;
    if(norm2>=cutoff2_)
      return 0;
  }
  return G.height*(std::exp(-0.5*norm2)-val_at_cutoff_);
}

template <class mode>
inline double OPESmetad<mode>::evaluateKernel(const kernel& G,const std::vector<double>& x, std::vector<double>& acc_der, std::vector<double>& dist)
{ //NB: cannot be a method of kernel class, because uses external variables (for cutoff)
  double norm2=0;
  for(unsigned i=0; i<ncv_; i++)
  {
    dist[i]=difference(i,G.center[i],x[i])/G.sigma[i];
    norm2+=dist[i]*dist[i];
    if(norm2>=cutoff2_)
      return 0;
  }
  const double val=G.height*(std::exp(-0.5*norm2)-val_at_cutoff_);
  for(unsigned i=0; i<ncv_; i++)
    acc_der[i]-=dist[i]/G.sigma[i]*val; //NB: we accumulate the derivative into der
  return val;
}

template <class mode>
inline void OPESmetad<mode>::mergeKernels(kernel& k1,const kernel& k2)
{
  const double h=k1.height+k2.height;
  for(unsigned i=0; i<k1.center.size(); i++)
  {
    const bool isPeriodic_i=getPntrToArgument(i)->isPeriodic();
    if(isPeriodic_i)
      k1.center[i]=k2.center[i]+difference(i,k2.center[i],k1.center[i]); //fix PBC
    const double c_i=(k1.height*k1.center[i]+k2.height*k2.center[i])/h;
    const double ss_k1_part=k1.height*(k1.sigma[i]*k1.sigma[i]+k1.center[i]*k1.center[i]);
    const double ss_k2_part=k2.height*(k2.sigma[i]*k2.sigma[i]+k2.center[i]*k2.center[i]);
    const double ss_i=(ss_k1_part+ss_k2_part)/h-c_i*c_i;
    if(isPeriodic_i)
      k1.center[i]=bringBackInPbc(i,c_i);
    else
      k1.center[i]=c_i;
    k1.sigma[i]=sqrt(ss_i);
  }
  k1.height=h;
}

}
}
