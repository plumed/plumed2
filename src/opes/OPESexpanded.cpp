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
#include "core/ActionSet.h"
#include "tools/Communicator.h"
#include "tools/File.h"

#include "ExpansionCVs.h"

namespace PLMD {
namespace opes {

//+PLUMEDOC OPES_BIAS OPES_EXPANDED
/*
On-the-fly probability enhanced sampling (\ref OPES "OPES") with expanded ensembles target distribution (replica-exchange-like) \cite Invernizzi2020unified.

An expanded ensemble is obtained by summing a set of ensembles at slightly different termodynamic conditions, or with slightly different Hamiltonians.
Such ensembles can be sampled via methods like replica exchange, or this OPES_EXPANDED bias action.
A typical example is a mutlticanonical simulation, in which a whole range of temperatures is sampled instead of a single one.

In oreder to define an expanded target ensemble we use \ref EXPANSION_CV "expansion collective variables" (ECVs), \f$\Delta u_\lambda(\mathbf{x})\f$.
The bias at step \f$n\f$ is
\f[
  V_n(\mathbf{x})=-\frac{1}{\beta}\log \left(\frac{1}{N_{\{\lambda\}}}\sum_\lambda e^{-\Delta u_\lambda(\mathbf{x})+\beta\Delta F_n(\lambda)}\right)\, .
\f]
See Ref.\cite Invernizzi2020unified for more details on the method.

Notice that the estimates in the DELTAFS file are expressed in energy units, and should be multiplied by \f$\beta\f$ to be dimensionless as in Ref.\cite Invernizzi2020unified.
The DELTAFS file also contains an estimate of \f$c(t)=\frac{1}{\beta} \log \langle e^{\beta V}\rangle\f$.
Similarly to \ref OPES_METAD, it is printed only for reference, since \f$c(t)\f$ reaches a fixed value as the bias converges, and should NOT be used for reweighting.
Its value is also needed for restarting a simulation.

Contrary to \ref OPES_METAD, OPES_EXPANDED does not use kernel density estimation.

\par Examples

\plumedfile
ene: ENERGY
ecv: ECV_MULTICANONICAL ARG=ene MAX_TEMP=1000
opes: OPES_EXPANDED ARG=ecv.* PACE=500
PRINT FILE=COLVAR STRIDE=500 ARG=ene,opes.bias
\endplumedfile

You can easily combine multiple ECVs.
The OPES_EXPANDED bias will create a multidimensional target grid to sample all the combinations.

\plumedfile
ene: ENERGY
dst: DISTANCE ATOMS=1,2

ecv1: ECV_MULTICANONICAL ARG=ene SET_ALL_TEMPS=200,300,500,1000
ecv2: ECV_UMBRELLAS_LINE ARG=dst MIN_CV=1.2 MAX_CV=4.3 SIGMA=0.5
opes: OPES_EXPANDED ARG=ecv1.*,ecv2.* PACE=500 OBSERVATION_STEPS=1

PRINT FILE=COLVAR STRIDE=500 ARG=ene,dst,opes.bias
\endplumedfile

If an ECV is based on more than one CV you must provide all the output component, in the proper order.
You can use regex for that, like in the following example.

\plumedfile
ene: ENERGY
vol: VOLUME
mtp: ECV_MULTITHERMAL_MULTIBARIC ...
  ARG=ene,vol
  TEMP=300
  MIN_TEMP=200
  MAX_TEMP=800
  PRESSURE=0.06022140857*1000 #1 kbar
  MIN_PRESSURE=0
  MAX_PRESSURE=0.06022140857*2000 #2 kbar
...

cv1: DISTANCE ATOMS=1,2
cv2: DISTANCE ATOMS=3,4
umb: ECV_UMBRELLAS_LINE ARG=cv1,cv2 TEMP=300 MIN_CV=0.1,0.1 MAX_CV=1.5,1.5 SIGMA=0.2 BARRIER=70

opes: OPES_EXPANDED ARG=mtp.*,umb.* PACE=500 WALKERS_MPI PRINT_STRIDE=1000

PRINT FILE=COLVAR STRIDE=500 ARG=ene,vol,cv1,cv2,opes.bias
\endplumedfile


*/
//+ENDPLUMEDOC

class OPESexpanded : public bias::Bias {

private:
  bool isFirstStep_;
  bool afterCalculate_;
  unsigned NumParallel_;
  unsigned rank_;
  unsigned NumWalkers_;
  unsigned long counter_;
  std::size_t ncv_;

  std::vector<const double *> ECVs_;
  std::vector<const double *> derECVs_;
  std::vector<opes::ExpansionCVs*> pntrToECVsClass_;
  std::vector< std::vector<unsigned> > index_k_;
// A note on indexes usage:
//  j -> underlying CVs
//  i -> DeltaFs
//  k -> single ECVs, which might not be trivially numbered
//  l -> groups of ECVs, pntrToECVsClass
//  h -> subgroups of ECVs, arguments in ECVsClass
//  w -> walkers

  double kbt_;
  unsigned stride_;
  std::vector<double> deltaF_;
  double rct_;
  double current_bias_;

  unsigned obs_steps_;
  std::vector<double> obs_cvs_;

  bool calc_work_;
  double work_;
  std::vector<double> old_deltaF_;

  unsigned print_stride_;
  OFile deltaFsOfile_;
  std::vector<std::string> deltaF_name_;

public:
  OPESexpanded(const ActionOptions&);
  void calculate() override;
  void update() override;
  static void registerKeywords(Keywords& keys);

  void init_pntrToECVsClass();
  void init_linkECVs();
  void init_from_obs();
  inline void update_deltaF(double);
  inline double getExpansion(const unsigned) const;
};

PLUMED_REGISTER_ACTION(OPESexpanded,"OPES_EXPANDED")

void OPESexpanded::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.remove("ARG");
  keys.add("compulsory","ARG","the label of the ECVs that define the expansion. You can use an * to make sure all the output components of the ECVs are used, as in the examples above");
  keys.add("compulsory","PACE","how often the bias is updated");
  keys.add("compulsory","OBSERVATION_STEPS","100","number of unbiased initial PACE steps to collect statistics for initialization");
//deltaFs file
  keys.add("compulsory","FILE","DELTAFS","a file with the estimate of the relative \\f$\\Delta F\\f$ for each component of the target and of the global \\f$c(t)\\f$");
  keys.add("compulsory","PRINT_STRIDE","100","stride for printing to DELTAFS file, in units of PACE");
  keys.add("optional","FMT","specify format for DELTAFS file");
//miscellaneous
  keys.addFlag("CALC_WORK",false,"calculate the work done by the bias between each update");
  keys.addFlag("WALKERS_MPI",false,"switch on MPI version of multiple walkers");
  keys.addFlag("SERIAL",false,"perform calculations in serial");
  keys.use("RESTART");

//output components
  componentsAreNotOptional(keys);
  keys.addOutputComponent("work","CALC_WORK","work done by the bias between each update");
}

OPESexpanded::OPESexpanded(const ActionOptions&ao)
  : PLUMED_BIAS_INIT(ao)
  , isFirstStep_(true)
  , afterCalculate_(false)
  , counter_(0)
  , ncv_(getNumberOfArguments())
  , rct_(0)
  , work_(0)
{
//set pace
  parse("PACE",stride_);
  parse("OBSERVATION_STEPS",obs_steps_);
  plumed_massert(obs_steps_!=0,"minimum is OBSERVATION_STEPS=1");
  obs_cvs_.resize(obs_steps_*ncv_);

//deltaFs file
  std::string deltaFsFileName;
  parse("FILE",deltaFsFileName);
  parse("PRINT_STRIDE",print_stride_);
  std::string fmt;
  parse("FMT",fmt);

//work flag
  parseFlag("CALC_WORK",calc_work_);

//multiple walkers //external MW for cp2k not supported, but anyway cp2k cannot put bias on energy!
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
    NumParallel_=1;
    rank_=0;
  }

  checkRead();

//check ECVs and link them
  init_pntrToECVsClass();
//set kbt_
  kbt_=pntrToECVsClass_[0]->getKbT();
  for(unsigned l=0; l<pntrToECVsClass_.size(); l++)
    plumed_massert(std::abs(kbt_-pntrToECVsClass_[l]->getKbT())<1e-4,"must set same TEMP for each ECV");

//restart if needed
  if(getRestart())
  {
    IFile ifile;
    ifile.link(*this);
    if(ifile.FileExist(deltaFsFileName))
    {
      log.printf("  RESTART from: %s\n",deltaFsFileName.c_str());
      log.printf(" +++ make sure all simulation options are consistent! +++\n");
      ifile.open(deltaFsFileName);
      //get deltaFs names
      ifile.scanFieldList(deltaF_name_);
      plumed_massert(deltaF_name_.size()>=4,"RESTART - fewer than expected FIELDS found in '"+deltaFsFileName+"' file");
      plumed_massert(deltaF_name_[deltaF_name_.size()-1]=="print_stride","RESTART - coult not find expected FIELDS in '"+deltaFsFileName+"' file");
      plumed_massert(deltaF_name_[0]=="time","RESTART - coult not find expected FIELDS in '"+deltaFsFileName+"' file");
      plumed_massert(deltaF_name_[1]=="rct","RESTART - coult not find expected FIELDS in '"+deltaFsFileName+"' file");
      deltaF_name_.pop_back();
      deltaF_name_.erase(deltaF_name_.begin(),deltaF_name_.begin()+2);
      std::size_t pos=5;//each name starts with "DeltaF"
      for(unsigned j=0; j<ncv_; j++)
        pos = deltaF_name_[0].find("_", pos+1); //checking only first one, hopefully is enough
      plumed_massert(pos<deltaF_name_[0].length(),"RESTART - fewer '_' than expected in DeltaF fields: did you remove any CV?");
      pos = deltaF_name_[0].find("_", pos+1);
      plumed_massert(pos>deltaF_name_[0].length(),"RESTART - more '_' than expected in DeltaF fields: did you add new CV?");
      //get lambdas, init ECVs and Link them
      auto getLambdaName=[](const std::string& name, const unsigned start, const unsigned dim)
      {
        std::size_t pos_start=5;//each name starts with "DeltaF"
        for(unsigned j=0; j<=start; j++)
          pos_start=name.find("_",pos_start+1);
        std::size_t pos_end=pos_start;
        for(unsigned j=0; j<dim; j++)
          pos_end=name.find("_",pos_end+1);
        pos_start++;//do not include heading "_"
        return name.substr(pos_start,pos_end-pos_start);
      };
      unsigned index_j=ncv_;
      unsigned sizeSkip=1;
      for(int l=pntrToECVsClass_.size()-1; l>=0; l--)
      {
        const unsigned dim_l=pntrToECVsClass_[l]->getNumberOfArguments();
        index_j-=dim_l;
        std::vector<std::string> lambdas_l(1);
        lambdas_l[0]=getLambdaName(deltaF_name_[0],index_j,dim_l);
        for(unsigned i=sizeSkip; i<deltaF_name_.size(); i+=sizeSkip)
        {
          std::string tmp_lambda=getLambdaName(deltaF_name_[i],index_j,dim_l);
          if(tmp_lambda==lambdas_l[0])
            break;
          lambdas_l.push_back(tmp_lambda);
        }
        pntrToECVsClass_[l]->initECVs_restart(lambdas_l);
        sizeSkip*=lambdas_l.size();
      }
      plumed_massert(sizeSkip==deltaF_name_.size(),"RESTART - this should not happen");
      deltaF_.resize(deltaF_name_.size());
      init_linkECVs(); //link ECVs and initializes index_k_
      log.printf(" ->%4lu DeltaFs in total\n",deltaF_.size());
      obs_steps_=0; //avoid initializing again
      //read steps from file
      unsigned restart_stride;
      ifile.scanField("print_stride",restart_stride);
      plumed_massert(restart_stride==print_stride_,"you can change PACE, but not PRINT_STRIDE. It would cause problems with multiple restarts");
      ifile.allowIgnoredFields(); //this allows for multiple restart, but without checking for consistency between them!
      double time;
      while(ifile.scanField("time",time)) //room for improvements: only number of lines and last line is important
      {
        if(calc_work_)
          old_deltaF_=deltaF_;
        ifile.scanField("rct",rct_);
        for(unsigned i=0; i<deltaF_.size(); i++)
          ifile.scanField(deltaF_name_[i],deltaF_[i]);
        ifile.scanField();
        counter_++;
      }
      log.printf("  successfully read %lu lines, up to t=%g\n",counter_,time);
      counter_=(1+(counter_-1)*print_stride_)*NumWalkers_; //adjust counter
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
  deltaFsOfile_.printField("print_stride",print_stride_);

//add output components
  if(calc_work_)
  {
    addComponent("work");
    componentIsNotPeriodic("work");
  }

//printing some info
  log.printf("  updating the bias with PACE = %u\n",stride_);
  log.printf("  initial unbiased OBSERVATION_STEPS = %u (in units of PACE)\n",obs_steps_);
  if(walkers_mpi)
    log.printf(" -- WALKERS_MPI: if multiple replicas are present, they will share the same bias via MPI\n");
  if(NumWalkers_>1)
  {
    log.printf("  using multiple walkers\n");
    log.printf("    number of walkers: %u\n",NumWalkers_);
    log.printf("    walker rank: %u\n",walker_rank);
  }
  int mw_warning=0;
  if(!walkers_mpi && comm.Get_rank()==0 && multi_sim_comm.Get_size()>(int)NumWalkers_)
    mw_warning=1;
  comm.Bcast(mw_warning,0);
  if(mw_warning) //log.printf messes up with comm, so never use it without Bcast!
    log.printf(" +++ WARNING +++ multiple replicas will NOT communicate unless the flag WALKERS_MPI is used\n");
  if(NumParallel_>1)
    log.printf("  using multiple threads per simulation: %u\n",NumParallel_);
  if(serial)
    log.printf(" -- SERIAL: running without loop parallelization\n");
//Bibliography
  log.printf("  Bibliography: ");
  log<<plumed.cite("M. Invernizzi, P.M. Piaggi, and M. Parrinello, arXiv:2007.03055 (2020)");
  log.printf("\n");
}

void OPESexpanded::calculate()
{
  if(deltaF_.size()==0) //no bias before initialization
    return;

  long double sum=0;
  std::vector<long double> der_sum_cv(ncv_,0);
  for(unsigned i=rank_; i<deltaF_.size(); i+=NumParallel_)
  {
    long double add_i=std::exp(static_cast<long double>(-getExpansion(i)+deltaF_[i]/kbt_));
    sum+=add_i;
    //set derivatives
    for(unsigned j=0; j<ncv_; j++)
      der_sum_cv[j]-=derECVs_[j][index_k_[i][j]]*add_i;
  }
  if(NumParallel_>1)
  {
    comm.Sum(sum);
    comm.Sum(der_sum_cv);
  }

  current_bias_=-kbt_*std::log(sum/deltaF_.size());
  setBias(current_bias_);
  for(unsigned j=0; j<ncv_; j++)
    setOutputForce(j,kbt_*der_sum_cv[j]/sum);

//calculate work
  if(calc_work_)
  {
    long double old_sum=0;
    for(unsigned i=rank_; i<deltaF_.size(); i+=NumParallel_)
      old_sum+=std::exp(static_cast<long double>(-getExpansion(i)+old_deltaF_[i]/kbt_));
    if(NumParallel_>1)
      comm.Sum(old_sum);
    work_+=-kbt_*std::log(sum/old_sum);
  }

  afterCalculate_=true;
}

void OPESexpanded::update()
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
    for(unsigned j=0; j<ncv_; j++)
      obs_cvs_[counter_*ncv_+j]=getArgument(j);
    counter_++;
    if(counter_==obs_steps_)
    {
      log.printf("\nAction %s\n",getName().c_str());
      init_from_obs();
      log.printf("Finished initialization\n\n");
      counter_=NumWalkers_; //all preliminary observations count 1
      obs_steps_=0; //no more observation
    }
    return;
  }
  plumed_massert(afterCalculate_,"OPESexpanded::update() must be called after OPESexpanded::calculate() to work properly");
  afterCalculate_=false;

//work done by the bias in one iteration
  if(calc_work_)
  {
    getPntrToComponent("work")->set(work_);
    work_=0;
    old_deltaF_=deltaF_;
  }

//update averages. This assumes that calculate() always runs before update(), thus uses current_bias_
  if(NumWalkers_==1)
    update_deltaF(current_bias_);
  else
  {
    std::vector<double> cvs(ncv_);
    for(unsigned j=0; j<ncv_; j++)
      cvs[j]=getArgument(j);
    std::vector<double> all_bias(NumWalkers_);
    std::vector<double> all_cvs(NumWalkers_*ncv_);
    if(comm.Get_rank()==0)
    {
      multi_sim_comm.Allgather(current_bias_,all_bias);
      multi_sim_comm.Allgather(cvs,all_cvs);
    }
    comm.Bcast(all_bias,0);
    comm.Bcast(all_cvs,0);
    for(unsigned w=0; w<NumWalkers_; w++)
    {
      //calculate ECVs
      unsigned index_wj=w*ncv_;
      for(unsigned k=0; k<pntrToECVsClass_.size(); k++)
      {
        pntrToECVsClass_[k]->calculateECVs(&all_cvs[index_wj]);
        index_wj+=pntrToECVsClass_[k]->getNumberOfArguments();
      }
      update_deltaF(all_bias[w]);
    }
  }

//write to file
  if((counter_/NumWalkers_-1)%print_stride_==0)
  {
    deltaFsOfile_.printField("time",getTime());
    deltaFsOfile_.printField("rct",rct_);
    for(unsigned i=0; i<deltaF_.size(); i++)
      deltaFsOfile_.printField(deltaF_name_[i],deltaF_[i]);
    deltaFsOfile_.printField();
  }
}

void OPESexpanded::init_pntrToECVsClass()
{
  std::vector<opes::ExpansionCVs*> all_pntrToECVsClass=plumed.getActionSet().select<opes::ExpansionCVs*>();
  plumed_massert(all_pntrToECVsClass.size()>0,"no Expansion CVs found");
  for(unsigned j=0; j<ncv_; j++)
  {
    std::string error_notECV("all the ARGs of "+getName()+" must be Expansion Collective Variables (ECV)");
    const unsigned dot_pos=getPntrToArgument(j)->getName().find(".");
    plumed_massert(dot_pos<getPntrToArgument(j)->getName().size(),error_notECV+", thus contain a dot in the name");
    unsigned foundECV_l=all_pntrToECVsClass.size();
    for(unsigned l=0; l<all_pntrToECVsClass.size(); l++)
    {
      if(getPntrToArgument(j)->getName().substr(0,dot_pos)==all_pntrToECVsClass[l]->getLabel())
      {
        foundECV_l=l;
        pntrToECVsClass_.push_back(all_pntrToECVsClass[l]);
        std::string missing_arg="some ECV component is missing from ARG";
        plumed_massert(j+all_pntrToECVsClass[l]->getNumberOfArguments()<=getNumberOfArguments(),missing_arg);
        for(unsigned h=0; h<all_pntrToECVsClass[l]->getNumberOfArguments(); h++)
        {
          std::string argName=getPntrToArgument(j+h)->getName();
          std::string expectedECVname=all_pntrToECVsClass[l]->getComponentsVector()[h];
          plumed_massert(argName==expectedECVname,missing_arg+", or is in wrong order: given ARG="+argName+" expected ARG="+expectedECVname);
        }
        j+=all_pntrToECVsClass[l]->getNumberOfArguments()-1;
        break;
      }
    }
    plumed_massert(foundECV_l<all_pntrToECVsClass.size(),error_notECV);
  }
  for(unsigned l=0; l<pntrToECVsClass_.size(); l++)
    for(unsigned ll=l+1; ll<pntrToECVsClass_.size(); ll++)
      plumed_massert(pntrToECVsClass_[l]->getLabel()!=pntrToECVsClass_[ll]->getLabel(),"cannot use same ECV twice");
}

void OPESexpanded::init_linkECVs()
{
  plumed_massert(deltaF_.size()>0,"init_linkECVs must run after deltaF_ has been resized");
  ECVs_.resize(ncv_);
  derECVs_.resize(ncv_);
  index_k_.resize(deltaF_.size(),std::vector<unsigned>(ncv_));
  unsigned index_j=0;
  unsigned sizeSkip=deltaF_.size();
  for(unsigned l=0; l<pntrToECVsClass_.size(); l++)
  {
    std::vector< std::vector<unsigned> > l_index_k(pntrToECVsClass_[l]->getIndex_k());
    plumed_massert(deltaF_.size()%l_index_k.size()==0,"buggy ECV: mismatch between getTotNumECVs() and getIndex_k().size()");
    plumed_massert(l_index_k[0].size()==pntrToECVsClass_[l]->getNumberOfArguments(),"buggy ECV: mismatch between number of ARG and underlying CVs");
    sizeSkip/=l_index_k.size();
    for(unsigned h=0; h<pntrToECVsClass_[l]->getNumberOfArguments(); h++)
    {
      ECVs_[index_j+h]=pntrToECVsClass_[l]->getPntrToECVs(h);
      derECVs_[index_j+h]=pntrToECVsClass_[l]->getPntrToDerECVs(h);
      for(unsigned i=0; i<deltaF_.size(); i++)
        index_k_[i][index_j+h]=l_index_k[(i/sizeSkip)%l_index_k.size()][h];
    }
    index_j+=pntrToECVsClass_[l]->getNumberOfArguments();
  }
  plumed_massert(sizeSkip==1,"this should not happen!");
}

void OPESexpanded::init_from_obs() //TODO improve speed?
{
//in case of multiple walkers gather all the statistics
  if(NumWalkers_>1)
  {
    std::vector<double> all_obs_cv(ncv_*obs_steps_*NumWalkers_);
    if(comm.Get_rank()==0)
      multi_sim_comm.Allgather(obs_cvs_,all_obs_cv);
    comm.Bcast(all_obs_cv,0);
    obs_cvs_=all_obs_cv;//could this lead to memory issues?
    obs_steps_*=NumWalkers_;
  }

//initialize ECVs from observations
  unsigned index_j=0;
  unsigned totNumECVs=1;
  for(unsigned l=0; l<pntrToECVsClass_.size(); l++)
  {
    pntrToECVsClass_[l]->initECVs_observ(obs_cvs_,ncv_,index_j);
    totNumECVs*=pntrToECVsClass_[l]->getTotNumECVs(); //ECVs from different exansions will be combined
    index_j+=pntrToECVsClass_[l]->getNumberOfArguments();
  }
  plumed_massert(index_j==getNumberOfArguments(),"mismatch between number of linked CVs and number of ARG");
//link ECVs and initialize index_k_, mapping each deltaF to a single ECVs set
  deltaF_.resize(totNumECVs);
  init_linkECVs();

//initialize deltaF_ from obs
//for the first point, t=0, the ECVs are calculated by initECVs_observ, setting also any initial guess
  index_j=0;
  for(unsigned i=0; i<deltaF_.size(); i++)
    for(unsigned j=0; j<ncv_; j++)
      deltaF_[i]+=kbt_*ECVs_[j][index_k_[i][j]];
  for(unsigned t=1; t<obs_steps_; t++) //starts from t=1
  {
    unsigned index_j=0;
    for(unsigned l=0; l<pntrToECVsClass_.size(); l++)
    {
      pntrToECVsClass_[l]->calculateECVs(&obs_cvs_[t*ncv_+index_j]);
      index_j+=pntrToECVsClass_[l]->getNumberOfArguments();
    }
    for(unsigned i=0; i<deltaF_.size(); i++)
    {
      const long double diff_i=static_cast<long double>(-getExpansion(i)+deltaF_[i]/kbt_);
      deltaF_[i]-=kbt_*(std::log1p(std::exp(diff_i)/t)+std::log1p(-1./(1.+t)));
    }
  }
  obs_cvs_.clear();
  if(calc_work_)
    old_deltaF_=deltaF_;

//set deltaF_name_
  deltaF_name_.resize(deltaF_.size(),"DeltaF");
  unsigned sizeSkip=deltaF_.size();
  for(unsigned l=0; l<pntrToECVsClass_.size(); l++)
  {
    std::vector<std::string> lambdas_l=pntrToECVsClass_[l]->getLambdas();
    plumed_massert(lambdas_l.size()==pntrToECVsClass_[l]->getTotNumECVs(),"buggy ECV: mismatch between getTotNumECVs() and getLambdas().size()");
    sizeSkip/=lambdas_l.size();
    for(unsigned i=0; i<deltaF_.size(); i++)
      deltaF_name_[i]+="_"+lambdas_l[(i/sizeSkip)%lambdas_l.size()];
  }

//print initialization to file
  log.printf(" ->%4lu DeltaFs in total\n",deltaF_.size());
  deltaFsOfile_.printField("time",getTime());
  deltaFsOfile_.printField("rct",rct_);
  for(unsigned i=0; i<deltaF_.size(); i++)
    deltaFsOfile_.printField(deltaF_name_[i],deltaF_[i]);
  deltaFsOfile_.printField();
}

inline void OPESexpanded::update_deltaF(double bias)
{
  //plumed_massert(counter_>1,"something went wrong");
  counter_++;
  const double increment=kbt_*std::log1p(std::exp(static_cast<long double>((bias-rct_)/kbt_))/(counter_-1.));
  for(unsigned i=0; i<deltaF_.size(); i++)
  {
    const long double diff_i=static_cast<long double>(-getExpansion(i)+(bias-rct_+deltaF_[i])/kbt_);
    deltaF_[i]+=increment-kbt_*std::log1p(std::exp(diff_i)/(counter_-1.));
  }
  rct_+=increment+kbt_*std::log1p(-1./counter_);
}

inline double OPESexpanded::getExpansion(unsigned i) const
{
  double expansion=0;
  for(unsigned j=0; j<ncv_; j++)
    expansion+=ECVs_[j][index_k_[i][j]];//the index_k could be trivially guessed for most ECVs, but unfourtunately not all
  return expansion;
}

}
}
