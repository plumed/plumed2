/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The VES code team
   (see the PEOPLE-VES file at the root of this folder for a list of names)

   See http://www.ves-code.org for more information.

   This file is part of VES code module.

   The VES code module is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The VES code module is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the VES code module.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include "bias/Bias.h"
#include "core/PlumedMain.h"
#include "core/ActionRegister.h"
#include "core/Atoms.h"
#include "tools/Communicator.h"
#include "tools/Grid.h"
#include "tools/File.h"
//#include <algorithm> //std::fill

namespace PLMD {
namespace ves {

//+PLUMEDOC VES_BIAS VES_DELTA_F
/*
Implementation of VES\f$\Delta F\f$ method \cite Invernizzi2019vesdeltaf (step two only).

\warning
  Notice that this is a stand-alone bias Action, it does not need any of the other VES module components

First you should create some estimate of the local free energy basins of your system,
using e.g. multiple \ref METAD short runs, and combining them with the \ref sum_hills utility.
Once you have them, you can use this bias Action to perform the VES optimization part of the method.

These \f$N+1\f$ local basins are used to model the global free energy.
In particular, given the conditional probabilities \f$P(\mathbf{s}|i)\propto e^{-\beta F_i(\mathbf{s})}\f$
and the probabilities of being in a given basin \f$P_i\f$, we can write:
\f[
  e^{-\beta F(\mathbf{s})}\propto P(\mathbf{s})=\sum_{i=0}^N P(\mathbf{s}|i)P_i \, .
\f]
We use this free energy model and the chosen bias factor \f$\gamma\f$ to build the bias potential:
\f$V(\mathbf{s})=-(1-1/\gamma)F(\mathbf{s})\f$.
Or, more explicitly:
\f[
  V(\mathbf{s})=(1-1/\gamma)\frac{1}{\beta}\log\left[e^{-\beta F_0(\mathbf{s})}
  +\sum_{i=1}^{N} e^{-\beta F_i(\mathbf{s})} e^{-\beta \alpha_i}\right] \, ,
\f]
where the parameters \f$\boldsymbol{\alpha}\f$ are the \f$N\f$ free energy differences (see below) from the \f$F_0\f$ basin.

By default the \f$F_i(\mathbf{s})\f$ are shifted so that \f$\min[F_i(\mathbf{s})]=0\f$ for all \f$i=\{0,...,N\}\f$.
In this case the optimization parameters \f$\alpha_i\f$ are the difference in height between the minima of the basins.
Using the keyword `NORMALIZE`, you can also decide to normalize the local free energies so that
\f$\int d\mathbf{s}\, e^{-\beta F_i(\mathbf{s})}=1\f$.
In this case the parameters will represent not the difference in height (which depends on the chosen CVs),
but the actual free energy difference, \f$\alpha_i=\Delta F_i\f$.

However, as discussed in Ref. \cite Invernizzi2019vesdeltaf, a better estimate of \f$\Delta F_i\f$ should be obtained through the reweighting procedure.

\par Examples

The following performs the optimization of the free energy difference between two metastable basins:

\plumedfile
VES_DELTA_F ...
  LABEL=ves
  ARG=cv
  TEMP=300
  FILE_F0=../fesA.data
  FILE_F1=../fesB.data
  BIASFACTOR=10.0
  M_STEP=0.1
  AV_STRIDE=500
  PRINT_STRIDE=100
... VES_DELTA_F

PRINT FMT=%g STRIDE=500 FILE=Colvar.data ARG=cv,ves.bias,ves.rct
\endplumedfile

The local FES files can be obtained as described in Sec. 4.2 of Ref. \cite Invernizzi2019vesdeltaf, i.e. for example:
- run 4 indipendent MetaD runs, all starting from basin A, and kill them as soon as they make the first transition (see e.g. \ref COMMITTOR)
- \verbatim cat HILLS* > all_HILLS \endverbatim
- \verbatim plumed sum_hills --hills all_HILLS --oufile all_fesA.dat --mintozero --min -1 --max 1 --bin 100 \endverbatim
- \verbatim awk -v n_rep=4 '{if($1!="#!" && $1!="") {for(i=1+(NF-1)/2; i<=NF; i++) $i/=n_rep;} print $0}' all_fesA.dat > fesA.data \endverbatim

The header of the file should be similar to the following:

\verbatim
#! FIELDS cv file.free der_cv
#! SET min_cv -1
#! SET max_cv 1
#! SET nbins_cv  100
#! SET periodic_cv false
\endverbatim

*/
//+ENDPLUMEDOC

class VesDeltaF : public bias::Bias {

private:
  double beta_;
  unsigned NumParallel_;
  unsigned rank_;
  unsigned NumWalkers_;
  bool isFirstStep_;
  bool afterCalculate_;

//prob
  double tot_prob_;
  std::vector<double> prob_;
  std::vector< std::vector<double> > der_prob_;

//local basins
  std::vector< std::unique_ptr<Grid> > grid_p_; //pointers because of GridBase::create
  std::vector<double> norm_;

//optimizer-related stuff
  long unsigned mean_counter_;
  unsigned mean_weight_tau_;
  unsigned alpha_size_;
  unsigned sym_alpha_size_;
  std::vector<double> mean_alpha_;
  std::vector<double> inst_alpha_;
  std::vector<double> past_increment2_;
  double minimization_step_;
  bool damping_off_;
//'tg' -> 'target distribution'
  double inv_gamma_;
  unsigned tg_counter_;
  unsigned tg_stride_;
  std::vector<double> tg_dV_dAlpha_;
  std::vector<double> tg_d2V_dAlpha2_;
//'av' -> 'ensemble average'
  unsigned av_counter_;
  unsigned av_stride_;
  std::vector<double> av_dV_dAlpha_;
  std::vector<double> av_dV_dAlpha_prod_;
  std::vector<double> av_d2V_dAlpha2_;
//printing
  unsigned print_stride_;
  OFile alphaOfile_;
//other
  std::vector<double> exp_alpha_;
  std::vector<double> prev_exp_alpha_;
  double work_;

//functions
  void update_alpha();
  void update_tg_and_rct();
  inline unsigned get_index(const unsigned, const unsigned) const;

public:
  explicit VesDeltaF(const ActionOptions&);
  void calculate() override;
  void update() override;
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(VesDeltaF,"VES_DELTA_F")

void VesDeltaF::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("optional","TEMP","temperature is compulsory, but it can be sometimes fetched from the MD engine");
//local free energies
  keys.add("numbered","FILE_F","names of files containing local free energies and derivatives. "
           "The first one, FILE_F0, is used as reference for all the free energy differences.");
  keys.reset_style("FILE_F","compulsory");
  keys.addFlag("NORMALIZE",false,"normalize all local free energies so that alpha will be (approx) \\f$\\Delta F\\f$");
  keys.addFlag("NO_MINTOZERO",false,"leave local free energies as provided, without shifting them to zero min");
//target distribution
  keys.add("compulsory","BIASFACTOR","0","the \\f$\\gamma\\f$ bias factor used for well-tempered target \\f$p(\\mathbf{s})\\f$."
           " Set to 0 for non-tempered flat target");
  keys.add("optional","TG_STRIDE","( default=1 ) number of AV_STRIDEs between updates"
           " of target \\f$p(\\mathbf{s})\\f$ and reweighing factor \\f$c(t)\\f$");
//optimization
  keys.add("compulsory","M_STEP","1.0","the \\f$\\mu\\f$ step used for the \\f$\\Omega\\f$ functional minimization");
  keys.add("compulsory","AV_STRIDE","500","number of simulation steps between alpha updates");
  keys.add("optional","TAU_MEAN","exponentially decaying average for alpha (rescaled using AV_STRIDE)."
           " Should be used only in very specific cases");
  keys.add("optional","INITIAL_ALPHA","( default=0 ) an initial guess for the bias potential parameter alpha");
  keys.addFlag("DAMPING_OFF",false,"do not use an AdaGrad-like term to rescale M_STEP");
//output parameters file
  keys.add("compulsory","ALPHA_FILE","ALPHA","file name for output minimization parameters");
  keys.add("optional","PRINT_STRIDE","( default=10 ) stride for printing to ALPHA_FILE");
  keys.add("optional","FMT","specify format for ALPHA_FILE");
//debug flags
  keys.addFlag("SERIAL",false,"perform the calculation in serial even if multiple tasks are available");
  keys.addFlag("MULTIPLE_WALKERS",false,"use multiple walkers connected via MPI for the optimization");
  keys.use("RESTART");

//output components
  componentsAreNotOptional(keys);
  keys.addOutputComponent("rct","default","the reweighting factor \\f$c(t)\\f$");
  keys.addOutputComponent("work","default","the work done by the bias in one AV_STRIDE");
}

VesDeltaF::VesDeltaF(const ActionOptions&ao)
  : PLUMED_BIAS_INIT(ao)
  , isFirstStep_(true)
  , afterCalculate_(false)
  , mean_counter_(0)
  , av_counter_(0)
  , work_(0)
{
//set beta
  const double Kb=plumed.getAtoms().getKBoltzmann();
  double temp=0;
  parse("TEMP",temp);
  double KbT=Kb*temp;
  if(KbT==0)
  {
    KbT=plumed.getAtoms().getKbT();
    plumed_massert(KbT>0,"your MD engine does not pass the temperature to plumed, you must specify it using TEMP");
  }
  beta_=1.0/KbT;

//initialize probability grids using local free energies
  bool spline=true;
  bool sparsegrid=false;
  std::string funcl="file.free"; //typical name given by sum_hills

  std::vector<std::string> fes_names;
  for(unsigned n=0;; n++)//NB: here we start from FILE_F0 not from FILE_F1
  {
    std::string filename;
    if(!parseNumbered("FILE_F",n,filename))
      break;
    fes_names.push_back(filename);
    IFile gridfile;
    gridfile.open(filename);
    auto g=GridBase::create(funcl,getArguments(),gridfile,sparsegrid,spline,true);
// we assume this cannot be sparse. in case we want it to be sparse, some of the methods
// that are available only in Grid should be ported to GridBase
    auto gg=dynamic_cast<Grid*>(g.get());
// if this throws, g is deleted
    plumed_assert(gg);
// release ownership in order to transfer it to emplaced pointer
    g.release();
    grid_p_.emplace_back(gg);
  }
  plumed_massert(grid_p_.size()>1,"at least 2 basins must be defined, starting from FILE_F0");
  alpha_size_=grid_p_.size()-1;
  sym_alpha_size_=alpha_size_*(alpha_size_+1)/2; //useful for symmetric matrix [alpha_size_]x[alpha_size_]
  //check for consistency with first local free energy
  for(unsigned n=1; n<grid_p_.size(); n++)
  {
    std::string error_tag="FILE_F"+std::to_string(n)+" '"+fes_names[n]+"' not compatible with reference one, FILE_F0";
    plumed_massert(grid_p_[n]->getSize()==grid_p_[0]->getSize(),error_tag);
    plumed_massert(grid_p_[n]->getMin()==grid_p_[0]->getMin(),error_tag);
    plumed_massert(grid_p_[n]->getMax()==grid_p_[0]->getMax(),error_tag);
    plumed_massert(grid_p_[n]->getBinVolume()==grid_p_[0]->getBinVolume(),error_tag);
  }

  bool no_mintozero=false;
  parseFlag("NO_MINTOZERO",no_mintozero);
  if(!no_mintozero)
  {
    for(unsigned n=0; n<grid_p_.size(); n++)
      grid_p_[n]->setMinToZero();
  }
  bool normalize=false;
  parseFlag("NORMALIZE",normalize);
  norm_.resize(grid_p_.size(),0);
  std::vector<double> c_norm(grid_p_.size());
  //convert the FESs to probability distributions
  //NB: the spline interpolation will be done on the probability distributions, not on the given FESs
  const unsigned ncv=getNumberOfArguments(); //just for ease
  for(unsigned n=0; n<grid_p_.size(); n++)
  {
    for(Grid::index_t t=0; t<grid_p_[n]->getSize(); t++)
    {
      std::vector<double> der(ncv);
      const double val=std::exp(-beta_*grid_p_[n]->getValueAndDerivatives(t,der));
      for(unsigned s=0; s<ncv; s++)
        der[s]*=-beta_*val;
      grid_p_[n]->setValueAndDerivatives(t,val,der);
      norm_[n]+=val;
    }
    c_norm[n]=1./beta_*std::log(norm_[n]);
    if(normalize)
    {
      grid_p_[n]->scaleAllValuesAndDerivatives(1./norm_[n]);
      norm_[n]=1;
    }
  }

//get target
  double biasfactor=0;
  parse("BIASFACTOR",biasfactor);
  plumed_massert(biasfactor==0 || biasfactor>1,"BIASFACTOR must be zero (for uniform target) or greater than one");
  if(biasfactor==0)
    inv_gamma_=0;
  else
    inv_gamma_=1./biasfactor;
  tg_counter_=0;
  tg_stride_=1;
  parse("TG_STRIDE",tg_stride_);
  tg_dV_dAlpha_.resize(alpha_size_,0);
  tg_d2V_dAlpha2_.resize(sym_alpha_size_,0);

//setup optimization stuff
  minimization_step_=1;
  parse("M_STEP",minimization_step_);

  av_stride_=500;
  parse("AV_STRIDE",av_stride_);
  av_dV_dAlpha_.resize(alpha_size_,0);
  av_dV_dAlpha_prod_.resize(sym_alpha_size_,0);
  av_d2V_dAlpha2_.resize(sym_alpha_size_,0);

  mean_weight_tau_=0;
  parse("TAU_MEAN",mean_weight_tau_);
  if(mean_weight_tau_!=1) //set it to 1 for basic SGD
  {
    plumed_massert((mean_weight_tau_==0 || mean_weight_tau_>av_stride_),"TAU_MEAN is rescaled with AV_STRIDE, so it has to be greater");
    mean_weight_tau_/=av_stride_; //this way you can look at the number of simulation steps to choose TAU_MEAN
  }

  parseVector("INITIAL_ALPHA",mean_alpha_);
  if(mean_alpha_.size()>0)
  {
    plumed_massert(mean_alpha_.size()==alpha_size_,"provide one INITIAL_ALPHA for each basin beyond the first one");
  }
  else
    mean_alpha_.resize(alpha_size_,0);
  inst_alpha_=mean_alpha_;
  exp_alpha_.resize(alpha_size_);
  for(unsigned i=0; i<alpha_size_; i++)
    exp_alpha_[i]=std::exp(-beta_*mean_alpha_[i]);
  prev_exp_alpha_=exp_alpha_;

  damping_off_=false;
  parseFlag("DAMPING_OFF",damping_off_);
  if(damping_off_)
    past_increment2_.resize(alpha_size_,1);
  else
    past_increment2_.resize(alpha_size_,0);

//file printing options
  std::string alphaFileName("ALPHA");
  parse("ALPHA_FILE",alphaFileName);
  print_stride_=10;
  parse("PRINT_STRIDE",print_stride_);
  std::string fmt;
  parse("FMT",fmt);

//other flags, mainly for debugging
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

  bool multiple_walkers=false;
  parseFlag("MULTIPLE_WALKERS",multiple_walkers);
  if(!multiple_walkers)
    NumWalkers_=1;
  else
  {
    if(comm.Get_rank()==0)//multi_sim_comm works well on first rank only
      NumWalkers_=multi_sim_comm.Get_size();
    if(comm.Get_size()>1) //if each walker has more than one processor update them all
      comm.Bcast(NumWalkers_,0);
  }

  checkRead();

//restart if needed
  if(getRestart())
  {
    IFile ifile;
    ifile.link(*this);
    if(NumWalkers_>1)
      ifile.enforceSuffix("");
    if(ifile.FileExist(alphaFileName))
    {
      log.printf("  Restarting from: %s\n",alphaFileName.c_str());
      log.printf("    all options (also PRINT_STRIDE) must be consistent!\n");
      log.printf("    any INITIAL_ALPHA will be overwritten\n");
      ifile.open(alphaFileName);
      double time;
      std::vector<double> damping(alpha_size_);
      while(ifile.scanField("time",time)) //room for improvements: only last line is important
      {
        for(unsigned i=0; i<alpha_size_; i++)
        {
          const std::string index(std::to_string(i+1));
          prev_exp_alpha_[i]=std::exp(-beta_*mean_alpha_[i]);
          ifile.scanField("alpha_"+index,mean_alpha_[i]);
          ifile.scanField("auxiliary_"+index,inst_alpha_[i]);
          ifile.scanField("damping_"+index,damping[i]);
        }
        ifile.scanField();
        mean_counter_+=print_stride_;
      }
      for(unsigned i=0; i<alpha_size_; i++)
      {
        exp_alpha_[i]=std::exp(-beta_*mean_alpha_[i]);
        past_increment2_[i]=damping[i]*damping[i];
      }
      //sync all walkers and treads. Not sure is mandatory but is no harm
      comm.Barrier();
      if(comm.Get_rank()==0)
        multi_sim_comm.Barrier();
    }
    else
      log.printf("  -- WARNING: restart requested, but no '%s' file found!\n",alphaFileName.c_str());
  }

//setup output file with Alpha values
  alphaOfile_.link(*this);
  if(NumWalkers_>1)
  {
    if(comm.Get_rank()==0 && multi_sim_comm.Get_rank()>0)
      alphaFileName="/dev/null"; //only first walker writes on file
    alphaOfile_.enforceSuffix("");
  }
  alphaOfile_.open(alphaFileName);
  if(fmt.length()>0)
    alphaOfile_.fmtField(" "+fmt);

//add other output components
  addComponent("rct"); componentIsNotPeriodic("rct");
  addComponent("work"); componentIsNotPeriodic("work");

//print some info
  log.printf("  Temperature T: %g\n",1./(Kb*beta_));
  log.printf("  Beta (1/Kb*T): %g\n",beta_);
  log.printf("  Local free energy basins files and normalization constants:\n");
  for(unsigned n=0; n<grid_p_.size(); n++)
    log.printf("    F_%d filename: %s  c_%d=%g\n",n,fes_names[n].c_str(),n,c_norm[n]);
  if(no_mintozero)
    log.printf(" -- NO_MINTOZERO: local free energies are not shifted to be zero at minimum\n");
  if(normalize)
    log.printf(" -- NORMALIZE: F_n+=c_n, alpha=DeltaF\n");
  log.printf("  Using target distribution with 1/gamma = %g\n",inv_gamma_);
  log.printf("    and updated with stride %d\n",tg_stride_);
  log.printf("  Step for the minimization algorithm: %g\n",minimization_step_);
  log.printf("  Stride for the ensemble average: %d\n",av_stride_);
  if(mean_weight_tau_>1)
    log.printf("  Exponentially decaying average with weight=tau/av_stride=%d\n",mean_weight_tau_);
  if(mean_weight_tau_==1)
    log.printf(" +++ WARNING +++ setting TAU_MEAN=1 is equivalent to use simple SGD, without mean alpha nor hessian contribution\n");
  log.printf("  Initial guess for alpha:\n");
  for(unsigned i=0; i<alpha_size_; i++)
    log.printf("    alpha_%d = %g\n",i+1,mean_alpha_[i]);
  if(damping_off_)
    log.printf(" -- DAMPING_OFF: the minimization step will NOT become smaller as the simulation goes on\n");
  log.printf("  Printing on file %s with stride %d\n",alphaFileName.c_str(),print_stride_);
  if(serial)
    log.printf(" -- SERIAL: running without loop parallelization\n");
  if(NumParallel_>1)
    log.printf("  Using multiple threads per simulation: %d\n",NumParallel_);
  if(multiple_walkers)
  {
    log.printf(" -- MULTIPLE_WALKERS: multiple simulations will combine statistics for the optimization\n");
    if(NumWalkers_>1)
    {
      log.printf("    number of walkers: %d\n",NumWalkers_);
      log.printf("    walker rank: %d\n",multi_sim_comm.Get_rank()); //only comm.Get_rank()=0 prints, so this is fine
    }
    else
      log.printf(" +++ WARNING +++ only one replica found: are you sure you are running MPI-connected simulations?\n");
  }
  log.printf(" Bibliography ");
  log<<plumed.cite("Invernizzi and Parrinello, J. Chem. Theory Comput. 15, 2187-2194 (2019)");
  log<<plumed.cite("Valsson and Parrinello, Phys. Rev. Lett. 113, 090601 (2014)");
  if(inv_gamma_>0)
    log<<plumed.cite("Valsson and Parrinello, J. Chem. Theory Comput. 11, 1996-2002 (2015)");

//last initializations
  prob_.resize(grid_p_.size());
  der_prob_.resize(grid_p_.size(),std::vector<double>(getNumberOfArguments()));
  update_tg_and_rct();
}

void VesDeltaF::calculate()
{
//get CVs
  const unsigned ncv=getNumberOfArguments(); //just for ease
  std::vector<double> cv(ncv);
  for(unsigned s=0; s<ncv; s++)
    cv[s]=getArgument(s);
//get probabilities for each basin, and total one
  for(unsigned n=0; n<grid_p_.size(); n++)
    prob_[n]=grid_p_[n]->getValueAndDerivatives(cv,der_prob_[n]);
  tot_prob_=prob_[0];
  for(unsigned i=0; i<alpha_size_; i++)
    tot_prob_+=prob_[i+1]*exp_alpha_[i];

//update bias and forces: V=-(1-inv_gamma_)*fes
  setBias((1-inv_gamma_)/beta_*std::log(tot_prob_));
  for(unsigned s=0; s<ncv; s++)
  {
    double dProb_dCV_s=der_prob_[0][s];
    for(unsigned i=0; i<alpha_size_; i++)
      dProb_dCV_s+=der_prob_[i+1][s]*exp_alpha_[i];
    setOutputForce(s,-(1-inv_gamma_)/beta_/tot_prob_*dProb_dCV_s);
  }
  afterCalculate_=true;
}

void VesDeltaF::update()
{
//skip first step to sync getTime() and av_counter_, as in METAD
  if(isFirstStep_)
  {
    isFirstStep_=false;
    return;
  }
  plumed_massert(afterCalculate_,"VesDeltaF::update() must be called after VesDeltaF::calculate() to work properly");
  afterCalculate_=false;

//calculate derivatives for ensemble averages
  std::vector<double> dV_dAlpha(alpha_size_);
  std::vector<double> d2V_dAlpha2(sym_alpha_size_);
  for(unsigned i=0; i<alpha_size_; i++)
    dV_dAlpha[i]=-(1-inv_gamma_)/tot_prob_*prob_[i+1]*exp_alpha_[i];
  for(unsigned i=0; i<alpha_size_; i++)
  {
    d2V_dAlpha2[get_index(i,i)]=-beta_*dV_dAlpha[i];
    for(unsigned j=i; j<alpha_size_; j++)
      d2V_dAlpha2[get_index(i,j)]-=beta_/(1-inv_gamma_)*dV_dAlpha[i]*dV_dAlpha[j];
  }
//update ensemble averages
  av_counter_++;
  for(unsigned i=0; i<alpha_size_; i++)
  {
    av_dV_dAlpha_[i]+=(dV_dAlpha[i]-av_dV_dAlpha_[i])/av_counter_;
    for(unsigned j=i; j<alpha_size_; j++)
    {
      const unsigned ij=get_index(i,j);
      av_dV_dAlpha_prod_[ij]+=(dV_dAlpha[i]*dV_dAlpha[j]-av_dV_dAlpha_prod_[ij])/av_counter_;
      av_d2V_dAlpha2_[ij]+=(d2V_dAlpha2[ij]-av_d2V_dAlpha2_[ij])/av_counter_;
    }
  }
//update work
  double prev_tot_prob=prob_[0];
  for(unsigned i=0; i<alpha_size_; i++)
    prev_tot_prob+=prob_[i+1]*prev_exp_alpha_[i];
  work_+=(1-inv_gamma_)/beta_*std::log(tot_prob_/prev_tot_prob);

//update coefficients
  if(av_counter_==av_stride_)
  {
    update_alpha();
    tg_counter_++;
    if(tg_counter_==tg_stride_)
    {
      update_tg_and_rct();
      tg_counter_=0;
    }
    //reset the ensemble averages
    av_counter_=0;
    std::fill(av_dV_dAlpha_.begin(),av_dV_dAlpha_.end(),0);
    std::fill(av_dV_dAlpha_prod_.begin(),av_dV_dAlpha_prod_.end(),0);
    std::fill(av_d2V_dAlpha2_.begin(),av_d2V_dAlpha2_.end(),0);
  }
}

void VesDeltaF::update_tg_and_rct()
{
//calculate target averages
  double Z_0=norm_[0];
  for(unsigned i=0; i<alpha_size_; i++)
    Z_0+=norm_[i+1]*exp_alpha_[i];
  double Z_tg=0;
  std::fill(tg_dV_dAlpha_.begin(),tg_dV_dAlpha_.end(),0);
  std::fill(tg_d2V_dAlpha2_.begin(),tg_d2V_dAlpha2_.end(),0);
  for(Grid::index_t t=rank_; t<grid_p_[0]->getSize(); t+=NumParallel_)
  { //TODO can we recycle some code?
    std::vector<double> prob(grid_p_.size());
    for(unsigned n=0; n<grid_p_.size(); n++)
      prob[n]=grid_p_[n]->getValue(t);
    double tot_prob=prob[0];
    for(unsigned i=0; i<alpha_size_; i++)
      tot_prob+=prob[i+1]*exp_alpha_[i];
    std::vector<double> dV_dAlpha(alpha_size_);
    std::vector<double> d2V_dAlpha2(sym_alpha_size_);
    for(unsigned i=0; i<alpha_size_; i++)
      dV_dAlpha[i]=-(1-inv_gamma_)/tot_prob*prob[i+1]*exp_alpha_[i];
    for(unsigned i=0; i<alpha_size_; i++)
    {
      d2V_dAlpha2[get_index(i,i)]=-beta_*dV_dAlpha[i];
      for(unsigned j=i; j<alpha_size_; j++)
        d2V_dAlpha2[get_index(i,j)]-=beta_/(1-inv_gamma_)*dV_dAlpha[i]*dV_dAlpha[j];
    }
    const double unnorm_tg_p=std::pow(tot_prob,inv_gamma_);
    Z_tg+=unnorm_tg_p;
    for(unsigned i=0; i<alpha_size_; i++)
      tg_dV_dAlpha_[i]+=unnorm_tg_p*dV_dAlpha[i];
    for(unsigned ij=0; ij<sym_alpha_size_; ij++)
      tg_d2V_dAlpha2_[ij]+=unnorm_tg_p*d2V_dAlpha2[ij];
  }
  if(NumParallel_>1)
  {
    comm.Sum(Z_tg);
    comm.Sum(tg_dV_dAlpha_);
    comm.Sum(tg_d2V_dAlpha2_);
  }
  for(unsigned i=0; i<alpha_size_; i++)
    tg_dV_dAlpha_[i]/=Z_tg;
  for(unsigned ij=0; ij<sym_alpha_size_; ij++)
    tg_d2V_dAlpha2_[ij]/=Z_tg;
  getPntrToComponent("rct")->set(-1./beta_*std::log(Z_tg/Z_0)); //Z_tg is the best available estimate of Z_V
}

void VesDeltaF::update_alpha()
{
//combining the averages of multiple walkers
  if(NumWalkers_>1)
  {
    if(comm.Get_rank()==0) //sum only once: in the first rank of each walker
    {
      multi_sim_comm.Sum(av_dV_dAlpha_);
      multi_sim_comm.Sum(av_dV_dAlpha_prod_);
      multi_sim_comm.Sum(av_d2V_dAlpha2_);
      for(unsigned i=0; i<alpha_size_; i++)
        av_dV_dAlpha_[i]/=NumWalkers_;
      for(unsigned ij=0; ij<sym_alpha_size_; ij++)
      {
        av_dV_dAlpha_prod_[ij]/=NumWalkers_;
        av_d2V_dAlpha2_[ij]/=NumWalkers_;
      }
    }
    if(comm.Get_size()>1)//if there are more ranks for each walker, everybody has to know
    {
      comm.Bcast(av_dV_dAlpha_,0);
      comm.Bcast(av_dV_dAlpha_prod_,0);
      comm.Bcast(av_d2V_dAlpha2_,0);
    }
  }
  //set work and reset it
  getPntrToComponent("work")->set(work_);
  work_=0;

//build the gradient and the Hessian of the functional
  std::vector<double> grad_omega(alpha_size_);
  std::vector<double> hess_omega(sym_alpha_size_);
  for(unsigned i=0; i<alpha_size_; i++)
  {
    grad_omega[i]=tg_dV_dAlpha_[i]-av_dV_dAlpha_[i];
    for(unsigned j=i; j<alpha_size_; j++)
    {
      const unsigned ij=get_index(i,j);
      hess_omega[ij]=beta_*(av_dV_dAlpha_prod_[ij]-av_dV_dAlpha_[i]*av_dV_dAlpha_[j])+tg_d2V_dAlpha2_[ij]-av_d2V_dAlpha2_[ij];
    }
  }
//calculate the increment and update alpha
  mean_counter_++;
  long unsigned mean_weight=mean_counter_;
  if(mean_weight_tau_>0 && mean_weight_tau_<mean_counter_)
    mean_weight=mean_weight_tau_;
  std::vector<double> damping(alpha_size_);
  for(unsigned i=0; i<alpha_size_; i++)
  {
    double increment_i=grad_omega[i];
    for(unsigned j=0; j<alpha_size_; j++)
      increment_i+=hess_omega[get_index(i,j)]*(inst_alpha_[j]-mean_alpha_[j]);
    if(!damping_off_)
      past_increment2_[i]+=increment_i*increment_i;
    damping[i]=std::sqrt(past_increment2_[i]);
    prev_exp_alpha_[i]=std::exp(-beta_*mean_alpha_[i]);
    inst_alpha_[i]-=minimization_step_/damping[i]*increment_i;
    mean_alpha_[i]+=(inst_alpha_[i]-mean_alpha_[i])/mean_weight;
    exp_alpha_[i]=std::exp(-beta_*mean_alpha_[i]);
  }

//update the Alpha file
  if(mean_counter_%print_stride_==0)
  {
    alphaOfile_.printField("time",getTime());
    for(unsigned i=0; i<alpha_size_; i++)
    {
      const std::string index(std::to_string(i+1));
      alphaOfile_.printField("alpha_"+index,mean_alpha_[i]);
      alphaOfile_.printField("auxiliary_"+index,inst_alpha_[i]);
      alphaOfile_.printField("damping_"+index,damping[i]);
    }
    alphaOfile_.printField();
  }
}

//mapping of a [alpha_size_]x[alpha_size_] symmetric matrix into a vector of size sym_alpha_size_, useful for the communicator
inline unsigned VesDeltaF::get_index(const unsigned i, const unsigned j) const
{
  if(i<=j)
    return j+i*(alpha_size_-1)-i*(i-1)/2;
  else
    return get_index(j,i);
}

}
}
