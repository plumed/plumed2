/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015,2016 The plumed team
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
/* 
 The sigma_mean optimisation method was developed by Thomas Loehr and
 Carlo Camilloni 
*/


#include "Bias.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "core/Value.h"
#include "tools/File.h"
#include "tools/OpenMP.h"
#include "tools/Random.h"
#include <cmath>
#include <ctime>

using namespace std;

namespace PLMD{
namespace bias{

//+PLUMEDOC BIAS METAINFERENCE
/*
Calculate the Metainference Score for a set of back calculated experimental data.


The back calculated data, that are expected to be averages over replicas (NR=1,2,..,N)
The functional form of this bias can be chosen between three variants selected
with NOISE=GAUSS,MGAUSS,OUTLIERS,MOUTLIERS which correspond to modelling the noise for
the arguments as a single gaussian common to all the data points, a gaussian per data
point, a lognormal per data point, or a single long-tailed gaussian common to all the data points.

As from Metainference theory there are two sigma values: SIGMA_MEAN represent the
error of calculating an average quanity using a finite set of replica and should
be set as small as possible following the guidelines for replica-averaged simulations
in the framework of the Maximum Entropy Principle. 
SIGMA_BIAS is an uncertainty parameter, sampled by a MC algorithm in the bounded interval 
defined by SIGMA_MIN and SIGMA_MAX. The initial value is set at SIGMA0. The MC move is a 
random displacement of maximum value equal to DSIGMA.

\par Examples

In the following example we calculate a set of \ref RDC, take the replica-average of
them and comparing them with a set of experimental values. RDCs are compared with
the experimental data but for a multiplication factor SCALE that is also sampled by
MC on-the-fly

\verbatim
RDC ...
LABEL=rdc
SCALE=0.0001
GYROM=-72.5388
ATOMS1=22,23
ATOMS2=25,27
ATOMS3=29,31
ATOMS4=33,34
... RDC

ardc: ENSEMBLE ARG=rdc.*

METAINFERENCE ...
ARG=ardc.*
NOISETYPE=MGAUSS
PARAMETERS=1.9190,2.9190,3.9190,4.9190
SCALEDATA SCALE0=1 SCALE_MIN=0.00001 SCALE_MAX=3 DSCALE=0.00 
SIGMA0=0.01 SIGMA_MIN=0.00001 SIGMA_MAX=3 DSIGMA=0.00 
SIGMA_MEAN=0.001
TEMP=300
LABEL=spe
... METAINFERENCE 

PRINT ARG=spe.bias FILE=BIAS STRIDE=1 
\endverbatim

in the following example instead of using one uncertainty parameter per data point we use
a single uncertainty value in a long-tailed gaussian to take into account for outliers.

\verbatim
METAINFERENCE ...
ARG=ardc.*
NOISETYPE=OUTLIERS
PARAMETERS=1.9190,2.9190,3.9190,4.9190
SCALEDATA SCALE0=1 SCALE_MIN=0.00001 SCALE_MAX=3 DSCALE=0.00 
SIGMA0=0.01 SIGMA_MIN=0.00001 SIGMA_MAX=3 DSIGMA=0.00 
SIGMA_MEAN=0.001
TEMP=300
LABEL=spe
... METAINFERENCE 
\endverbatim

(See also \ref RDC and \ref ENSEMBLE).

*/
//+ENDPLUMEDOC

class Metainference : public Bias
{
  // experimental values
  vector<double> parameters;
  // noise type
  unsigned noise_type_;
  enum { GAUSS, MGAUSS, OUTLIERS, MOUTLIERS };
  // scale is data scaling factor
  // noise type
  unsigned scale_prior_;
  enum { SC_GAUSS, SC_FLAT };
  bool   doscale_;
  double scale_;
  double scale_mu_;
  double scale_sigma_;
  double scale_min_;
  double scale_max_;
  double Dscale_;
  // sigma is data uncertainty
  vector<double> sigma_;
  vector<double> sigma_min_;
  vector<double> sigma_max_;
  vector<double> Dsigma_;
  // sigma_mean is uncertainty in the mean estimate
  vector<double> sigma_mean_;
  vector<double> variance_;

  // sigma_mean rescue params
  double sm_mod_;
  double sm_mod_min_;
  double sm_mod_max_;
  double Dsm_mod_;
  double max_force_;

  // temperature in kbt
  double   kbt_;
  // number of data points
  unsigned ndata_;

  // Monte Carlo stuff
  vector<Random> random;
  unsigned MCsteps_;
  unsigned MCstride_;
  long unsigned MCaccept_;
  long unsigned MCtrial_;
  unsigned MCchunksize_;

  // output
  Value*   valueScale;
  Value*   valueAccept;
  Value*   valueRSigmaMean;
  vector<Value*> valueSigma;
  vector<Value*> valueSigmaMean;
  Value*   valueSMmod;
  Value*   valueMaxForceMD;

  // restart
  unsigned write_stride_;
  OFile    sfile_;

  // others
  bool     master;
  bool     do_reweight;
  unsigned do_optsigmamean_;
  unsigned nrep_;
  unsigned replica_;
  unsigned narg;

  // we need this for the forces
  Atoms& atoms;

  double getEnergySP(const vector<double> &mean, const vector<double> &sigma, const double scale, const double modifier);
  double getEnergySPE(const vector<double> &mean, const vector<double> &sigma, const double scale, const double modifier);
  double getEnergyGJ(const vector<double> &mean, const vector<double> &sigma, const double scale, const double modifier);
  double getEnergyGJE(const vector<double> &mean, const vector<double> &sigma, const double scale, const double modifier);
  void   doMonteCarlo(const vector<double> &mean, const double modifier);
  double getEnergyForceSP(const vector<double> &mean, const double fact, const double modifier);
  double getEnergyForceSPE(const vector<double> &mean, const double fact, const double modifier);
  double getEnergyForceGJ(const vector<double> &mean, const double fact, const double modifier);
  double getEnergyForceGJE(const vector<double> &mean, const double fact, const double modifier);
  void   writeStatus();
  
public:
  explicit Metainference(const ActionOptions&);
  ~Metainference();
  void calculate();
  void update();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(Metainference,"METAINFERENCE")

void Metainference::registerKeywords(Keywords& keys){
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("optional","PARARG","reference values for the experimental data, these can be provided as arguments without derivatives"); 
  keys.add("optional","PARAMETERS","reference values for the experimental data");
  keys.add("compulsory","NOISETYPE","functional form of the noise (GAUSS,MGAUSS,OUTLIERS,MOUTLIERS)");
  keys.addFlag("REWEIGHT",false,"simple REWEIGHT using the latest ARG as energy"); 
  keys.addFlag("SCALEDATA",false,"Set to TRUE if you want to sample a scaling factor common to all values and replicas");  
  keys.add("compulsory","SCALE0","initial value of the uncertainty parameter");
  keys.add("compulsory","SCALE_PRIOR","FLAT","either FLAT or GAUSSIAN");
  keys.add("optional","SCALE_MIN","minimum value of the uncertainty parameter");
  keys.add("optional","SCALE_MAX","maximum value of the uncertainty parameter");
  keys.add("optional","SCALE_SIGMA","maximum value of the uncertainty parameter");
  keys.add("optional","DSCALE","maximum MC move of the uncertainty parameter");
  keys.add("compulsory","SIGMA0","initial value of the uncertainty parameter");
  keys.add("compulsory","SIGMA_MIN","minimum value of the uncertainty parameter");
  keys.add("compulsory","SIGMA_MAX","maximum value of the uncertainty parameter");
  keys.add("optional","DSIGMA","maximum MC move of the uncertainty parameter");
  keys.add("compulsory","OPTSIGMAMEAN","NONE","Set to NONE/SEM/FULL to manually set sigma mean, or to estimate it on the fly and use the safety check on forces");  
  keys.add("optional","SIGMA_MEAN0","starting value for the uncertainty in the mean estimate");
  keys.add("optional","SIGMA_MEAN_MOD0","starting value for sm modifier");
  keys.add("optional","SIGMA_MEAN_MOD_MIN","starting value for sm modifier");
  keys.add("optional","SIGMA_MEAN_MOD_MAX","starting value for sm modifier");
  keys.add("optional","DSIGMA_MEAN_MOD","step value for sm modifier");
  keys.add("optional","MAX_FORCE","maximum allowable force");
  keys.add("optional","TEMP","the system temperature - this is only needed if code doesnt' pass the temperature to plumed");
  keys.add("optional","MC_STEPS","number of MC steps");
  keys.add("optional","MC_STRIDE","MC stride");
  keys.add("optional","MC_CHUNKSIZE","MC chunksize");
  keys.add("optional","STATUS_FILE","write a file with all the data usefull for restart/continuation of Metainference");
  keys.add("compulsory","WRITE_STRIDE","write the status to a file every N steps, this can be used for restart/continuation");
  keys.use("RESTART");
  useCustomisableComponents(keys);
  keys.addOutputComponent("sigma",        "default",      "uncertainty parameter");
  keys.addOutputComponent("sigmaMean",    "default",      "uncertainty in the mean estimate");
  keys.addOutputComponent("rewSigmaMean", "default",      "sigma mean multiplier (give by scale, reweight and optimisation)");
  keys.addOutputComponent("accept",       "default",      "MC acceptance");
  keys.addOutputComponent("weight",       "REWEIGHT",     "weights of the weighted average");
  keys.addOutputComponent("MetaDf",       "REWEIGHT",     "force on metadynamics");
  keys.addOutputComponent("scale",        "SCALEDATA",    "scale parameter");
  keys.addOutputComponent("maxForceMD",   "OPTSIGMAMEAN", "max force on atoms");
  keys.addOutputComponent("smMod",        "OPTSIGMAMEAN", "modifier for all sigma mean");
}

Metainference::Metainference(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao), 
doscale_(false),
scale_mu_(0),
scale_sigma_(-1),
scale_min_(1),
scale_max_(-1),
Dscale_(-1),
sm_mod_(1.),
ndata_(getNumberOfArguments()),
random(3),
MCsteps_(1), 
MCstride_(1), 
MCaccept_(0), 
MCtrial_(0),
MCchunksize_(0),
write_stride_(0),
do_reweight(false),
do_optsigmamean_(0),
atoms(plumed.getAtoms())
{
  // set up replica stuff 
  master = (comm.Get_rank()==0);
  if(master) {
    nrep_    = multi_sim_comm.Get_size();
    replica_ = multi_sim_comm.Get_rank();
  } else {
    nrep_    = 0;
    replica_ = 0;
  }
  comm.Sum(&nrep_,1);
  comm.Sum(&replica_,1);

  // reweight implies a different number of arguments (the latest one must always be the bias) 
  parseFlag("REWEIGHT", do_reweight); 
  narg = getNumberOfArguments();
  if(do_reweight) narg--;

  parseVector("PARAMETERS",parameters);
  if(parameters.size()!=static_cast<unsigned>(narg)&&!parameters.empty())
    error("Size of PARAMETERS array should be either 0 or the same as of the number of arguments in ARG1");

  vector<Value*> arg2;
  parseArgumentList("PARARG",arg2);
  if(!arg2.empty()) {
    if(parameters.size()>0) error("It is not possible to use PARARG and PARAMETERS together");
    if(arg2.size()!=narg) error("Size of PARARG array should be the same as number for arguments in ARG");
    for(unsigned i=0;i<arg2.size();i++){
      parameters.push_back(arg2[i]->get()); 
      if(arg2[i]->hasDerivatives()==true) error("PARARG can only accept arguments without derivatives");
    }
  }

  if(parameters.size()!=narg) 
    error("PARARG or PARAMETERS arrays should include the same number of elements as the arguments in ARG");

  string stringa_noise;
  parse("NOISETYPE",stringa_noise);
  if(stringa_noise=="GAUSS")           noise_type_ = GAUSS; 
  else if(stringa_noise=="MGAUSS")     noise_type_ = MGAUSS;
  else if(stringa_noise=="OUTLIERS")   noise_type_ = OUTLIERS;
  else if(stringa_noise=="MOUTLIERS")  noise_type_ = MOUTLIERS;
  else error("Unknown noise type!"); 

  parse("WRITE_STRIDE",write_stride_);
  string status_file_name_;
  parse("STATUS_FILE",status_file_name_);
  if(status_file_name_=="") status_file_name_ = "MISTATUS"+getLabel();

  string stringa_optsigma;
  parse("OPTSIGMAMEAN", stringa_optsigma);
  if(stringa_optsigma=="NONE")      do_optsigmamean_=0;
  else if(stringa_optsigma=="SEM")  do_optsigmamean_=1;
  else if(stringa_optsigma=="FULL") do_optsigmamean_=2;

  parseFlag("SCALEDATA", doscale_);
  if(doscale_) {
    string stringa_noise;
    parse("SCALE_PRIOR",stringa_noise);
    if(stringa_noise=="GAUSSIAN")  scale_prior_ = SC_GAUSS; 
    else if(stringa_noise=="FLAT") scale_prior_ = SC_FLAT;
    else error("Unknown SCALE_PRIOR type!");
    parse("SCALE0",scale_);
    if(scale_prior_==SC_GAUSS) {
      parse("SCALE_SIGMA",scale_sigma_);
      scale_mu_=scale_;
    } else {
      parse("SCALE_MIN",scale_min_);
      parse("SCALE_MAX",scale_max_);
    }
    if(scale_prior_==SC_GAUSS&&scale_sigma_<0.) 
      error("SCALE_SIGMA must be set when using SCALE_PRIOR=GAUSS");
    if(scale_prior_==SC_FLAT&&scale_max_<scale_min_) 
      error("SCALE_MAX and SCALE_MIN must be set when using SCALE_PRIOR=FLAT");
    parse("DSCALE",Dscale_);
    if(Dscale_<0) {
      if(scale_prior_==SC_FLAT) Dscale_ = 0.05*(scale_max_ - scale_min_);
    }
  } else {
    scale_=1.0;
  }

  vector<double> readsigma;
  parseVector("SIGMA0",readsigma);
  if((noise_type_!=MGAUSS&&noise_type_!=MOUTLIERS)&&readsigma.size()>1) 
    error("If you want to use more than one SIGMA you should use NOISETYPE=MGAUSS|MOUTLIERS");
  if(noise_type_==MGAUSS||noise_type_==MOUTLIERS) {
    if(readsigma.size()==narg) {
      sigma_.resize(narg);
      sigma_=readsigma;
    } else if(readsigma.size()==1) {
      sigma_.resize(narg,readsigma[0]);
    } else {
      error("SIGMA0 can accept either one single value or as many values as the number of arguments (with NOISETYPE=MGAUSS|MOUTLIERS)");
    } 
  } else sigma_.resize(1, readsigma[0]);

  double read_smin_;
  parse("SIGMA_MIN",read_smin_);
  sigma_min_.resize(sigma_.size(),read_smin_);
  double read_smax_;
  parse("SIGMA_MAX",read_smax_);
  sigma_max_.resize(sigma_.size(),read_smax_);
  double read_dsigma_=-1.;
  parse("DSIGMA",read_dsigma_);
  if(read_dsigma_<0) read_dsigma_ = 0.05*(read_smax_ - read_smin_);
  Dsigma_.resize(sigma_.size(),read_dsigma_);

  // monte carlo stuff
  parse("MC_STEPS",MCsteps_);
  parse("MC_STRIDE",MCstride_);
  parse("MC_CHUNKSIZE", MCchunksize_);
  // adjust for multiple-time steps
  MCstride_ *= getStride();
  // get temperature
  double temp=0.0;
  parse("TEMP",temp);
  if(temp>0.0) kbt_=plumed.getAtoms().getKBoltzmann()*temp;
  else kbt_=plumed.getAtoms().getKbT();

  // variance is always the size of narg
  variance_.resize(narg,0);
  // while sigma_mean_ has the same size of sigma
  vector<double> read_sigma_mean_;
  parseVector("SIGMA_MEAN0",read_sigma_mean_);
  if(!do_optsigmamean_ && read_sigma_mean_.size()==0 && !getRestart()) 
    error("If you don't use OPTSIGMAMEAN and you are not RESTARTING then you MUST SET SIGMA_MEAN0");

  if(noise_type_==MGAUSS||noise_type_==MOUTLIERS) {
    if(read_sigma_mean_.size()==narg) {
      sigma_mean_.resize(narg);
      sigma_mean_=read_sigma_mean_;
    } else if(read_sigma_mean_.size()==1) {
      sigma_mean_.resize(narg,read_sigma_mean_[0]);
    } else if(read_sigma_mean_.size()==0) {
      sigma_mean_.resize(narg,0.000001);
    } else {
      error("SIGMA_MEAN0 can accept either one single value or as many values as the arguments (with NOISETYPE=MGAUSS|MOUTLIERS)");
    }
    for(unsigned i=0;i<narg;i++) if(sigma_mean_[i]>0) variance_[i] = sigma_mean_[i]*sigma_mean_[i]*static_cast<double>(nrep_);
  } else {
    if(read_sigma_mean_.size()==1) {
      sigma_mean_.resize(1, read_sigma_mean_[0]);
    } else if(read_sigma_mean_.size()==0) {
      sigma_mean_.resize(narg,0.000001);
    } else {
      error("If you want to use more than one SIGMA_MEAN0 you should use NOISETYPE=MGAUSS|MOUTLIERS");
    }
    for(unsigned i=0;i<narg;i++) variance_[i] = sigma_mean_[0]*sigma_mean_[0]*static_cast<double>(nrep_);
  } 

  // sigma mean optimisation
  if(do_optsigmamean_==2) {
    max_force_=3000.;
    parse("MAX_FORCE", max_force_);
    max_force_ *= max_force_;
    sm_mod_=1.0;
    parse("SIGMA_MEAN_MOD0", sm_mod_);
    sm_mod_min_=1.0;
    parse("SIGMA_MEAN_MOD_MIN", sm_mod_min_);
    sm_mod_max_=sqrt(10.);
    parse("SIGMA_MEAN_MOD_MAX", sm_mod_max_);
    Dsm_mod_=0.01;
    parse("DSIGMA_MEAN_MOD", Dsm_mod_);
  }

  checkRead();

  IFile restart_sfile;
  restart_sfile.link(*this);
  if(getRestart()&&restart_sfile.FileExist(status_file_name_)) {
    restart_sfile.open(status_file_name_);
    log.printf("  Restarting from %s\n", status_file_name_.c_str());
    double dummy;
    if(restart_sfile.scanField("time",dummy)){
      for(unsigned i=0;i<variance_.size();++i) {
        std::string msg;
        Tools::convert(i,msg);
        restart_sfile.scanField("variance_"+msg,variance_[i]);
      }
      for(unsigned i=0;i<sigma_.size();++i) {
        std::string msg;
        Tools::convert(i,msg);
        restart_sfile.scanField("sigma_"+msg,sigma_[i]);
      }
      if(doscale_) restart_sfile.scanField("scale0_",scale_);
      if(do_optsigmamean_==2) {
        restart_sfile.scanField("sigma_mean_mod0",sm_mod_);
      }
    }
    restart_sfile.scanField();
    restart_sfile.close();
    /* set sigma_mean from variance */
    if(noise_type_==MGAUSS||noise_type_==MOUTLIERS) {
      for(unsigned i=0;i<variance_.size();++i) sigma_mean_[i] = sqrt(variance_[i]/static_cast<double>(nrep_));
    } else {
      double s_v = *max_element(variance_.begin(), variance_.end());
      sigma_mean_[0] = sqrt(s_v/static_cast<double>(nrep_));
    }
  }

  switch(noise_type_) {
    case GAUSS:
      log.printf("  with gaussian noise and a single noise parameter for all the data\n");
      break;
    case MGAUSS:
      log.printf("  with gaussian noise and a noise parameter for each data point\n");
      break;
    case OUTLIERS:
      log.printf("  with long tailed gaussian noise and a single noise parameter for all the data\n");
      break;
    case MOUTLIERS:
      log.printf("  with long tailed gaussian noise and a noise parameter for each data point\n");
      break;
  }

  if(doscale_) {
    log.printf("  sampling a common scaling factor with:\n");
    log.printf("    initial scale parameter %f\n",scale_);
    log.printf("    minimum scale parameter %f\n",scale_min_);
    log.printf("    maximum scale parameter %f\n",scale_max_);
    log.printf("    maximum MC move of scale parameter %f\n",Dscale_);
  }

  log.printf("  number of experimental data points %u\n",narg);
  log.printf("  number of replicas %u\n",nrep_);
  log.printf("  initial data uncertainties");
  for(unsigned i=0;i<sigma_.size();++i) log.printf(" %f", sigma_[i]);
  log.printf("\n");
  log.printf("  minimum data uncertainty %f\n",read_smin_);
  log.printf("  maximum data uncertainty %f\n",read_smax_);
  log.printf("  maximum MC move of data uncertainty %f\n",read_dsigma_);
  log.printf("  temperature of the system %f\n",kbt_);
  log.printf("  MC steps %u\n",MCsteps_);
  log.printf("  MC stride %u\n",MCstride_);
  log.printf("  initial standard errors of the mean");
  for(unsigned i=0;i<sigma_mean_.size();++i) log.printf(" %f", sigma_mean_[i]);
  log.printf("\n");

  if(do_reweight) {
    addComponent("MetaDf");
    componentIsNotPeriodic("MetaDf");
    valueRSigmaMean=getPntrToComponent("MetaDf");
    addComponent("weight");
    componentIsNotPeriodic("weight");
    valueRSigmaMean=getPntrToComponent("weight");
  }

  if(doscale_) { 
    addComponent("scale");  
    componentIsNotPeriodic("scale");
    valueScale=getPntrToComponent("scale");
  }

  addComponent("accept");
  componentIsNotPeriodic("accept");
  valueAccept=getPntrToComponent("accept");
  
  addComponent("rewSigmaMean");
  componentIsNotPeriodic("rewSigmaMean");
  valueRSigmaMean=getPntrToComponent("rewSigmaMean");

  if(do_optsigmamean_==2) {
    addComponent("maxForceMD");
    componentIsNotPeriodic("maxForceMD");
    valueMaxForceMD=getPntrToComponent("maxForceMD");
    addComponent("smMod");
    componentIsNotPeriodic("smMod");
    valueSMmod=getPntrToComponent("smMod");
  }

  if(noise_type_==MGAUSS||noise_type_==MOUTLIERS) {
    for(unsigned i=0;i<sigma_mean_.size();++i){
      std::string num; Tools::convert(i,num);
      addComponent("sigmaMean_"+num); componentIsNotPeriodic("sigmaMean_"+num);
      valueSigmaMean.push_back(getPntrToComponent("sigmaMean_"+num));
      getPntrToComponent("sigmaMean_"+num)->set(sigma_mean_[i]);
      addComponent("sigma_"+num); componentIsNotPeriodic("sigma_"+num);
      valueSigma.push_back(getPntrToComponent("sigma_"+num));
      getPntrToComponent("sigma_"+num)->set(sigma_[i]);
    }
  } else {
    addComponent("sigmaMean"); componentIsNotPeriodic("sigmaMean");
    valueSigmaMean.push_back(getPntrToComponent("sigmaMean"));
    getPntrToComponent("sigmaMean")->set(sigma_mean_[0]);
    addComponent("sigma"); componentIsNotPeriodic("sigma");
    valueSigma.push_back(getPntrToComponent("sigma"));
    getPntrToComponent("sigma")->set(sigma_[0]);
  }

  // initialize random seed
  unsigned iseed;
  if(master) iseed = time(NULL)+replica_;
  else iseed = 0;     
  comm.Sum(&iseed, 1);
  random[0].setSeed(-iseed);
  // Random chunk
  if(master) iseed = time(NULL)+replica_;
  else iseed = 0;     
  comm.Sum(&iseed, 1);
  random[2].setSeed(-iseed);
  if(doscale_) {
    // in this case we want the same seed everywhere
    iseed = time(NULL);
    if(master) multi_sim_comm.Bcast(iseed,0);
    comm.Bcast(iseed,0);
    random[1].setSeed(-iseed);
  }

  // outfile stuff
  if(write_stride_>0) {
    sfile_.link(*this);
    sfile_.open(status_file_name_);
  }

  log<<"  Bibliography "<<plumed.cite("Bonomi, Camilloni, Cavalli, Vendruscolo, Sci. Adv. 2, e150117 (2016)");
  if(do_reweight) log<<plumed.cite("Bonomi, Camilloni, Vendruscolo, Sci. Rep. 6, 31232 (2016)");
  log<<"\n";
}

Metainference::~Metainference()
{
  if(sfile_.isOpen()) sfile_.close();
}

double Metainference::getEnergySP(const vector<double> &mean, const vector<double> &sigma, const double scale, const double modifier)
{
  const double sm2 = scale*scale*sigma_mean_[0]*sigma_mean_[0]*modifier*modifier;
  const double ss = sigma[0]*sigma[0] + sm2;

  double ene = 0.0;
  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(ene)
  {
    #pragma omp for reduction( + : ene)
    for(unsigned i=0;i<narg;++i){
      const double dev = scale*mean[i]-parameters[i]; 
      const double a2 = 0.5*dev*dev + ss;
      ene += std::log( 2.0 * a2 / ( 1.0 - exp(- a2 / sm2) ) );
    }
  }
  // add one single Jeffrey's prior and one normalisation per data point
  ene += std::log(sigma[0]+sigma_mean_[0]*modifier) - static_cast<double>(ndata_)*0.5*std::log(2./(M_PI*M_PI)*ss);
  return kbt_ * ene;
}

double Metainference::getEnergySPE(const vector<double> &mean, const vector<double> &sigma, const double scale, const double modifier)
{
  const double mul_sm2 = scale*scale*modifier*modifier;
  const double m2 = modifier*modifier;
  double ene = 0.0;
  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(ene)
  {
    #pragma omp for reduction( + : ene)
    for(unsigned i=0;i<narg;++i){
      const double sm2 = mul_sm2*sigma_mean_[i]*sigma_mean_[i];
      const double ss  = sigma[i]*sigma[i] + sm2;
      const double sss = sigma[i]*sigma[i] + sigma_mean_[i]*sigma_mean_[i]*m2;
      const double dev = scale*mean[i]-parameters[i]; 
      const double a2  = 0.5*dev*dev + ss;
      // deviation + jeffreys + normalisation
      // ene += 0.5*std::log(sss) - 0.5*std::log(sqrt2_div_pi*ss);
      // this is equivalent to (but with less log)
      const double add = 0.5*std::log(0.5*M_PI*M_PI*sss/ss); 
      ene += std::log( 2.0 * a2 / ( 1.0 - exp(- a2 / sm2) ) ) + add;
    }
  }
  // add normalization and Jeffrey's prior
  return kbt_ * ene;
}

double Metainference::getEnergyGJ(const vector<double> &mean, const vector<double> &sigma, const double scale, const double modifier)
{
  const double inv_s2 = 1./(sigma[0]*sigma[0] + modifier*modifier*scale*scale*sigma_mean_[0]*sigma_mean_[0]);
  const double inv_sss = 1./(sigma[0]*sigma[0] + modifier*modifier*sigma_mean_[0]*sigma_mean_[0]);

  double ene = 0.0;
  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(ene)
  {
    #pragma omp for reduction( + : ene)
    for(unsigned i=0;i<narg;++i){
      double dev = scale*mean[i]-parameters[i];
      ene += 0.5*dev*dev*inv_s2;
    }
  }
  // add Jeffrey's prior in case one sigma for all data points + one normalisation per datapoint
  // scale doens't enter in the jeffrey, because we don't want the scale to be biased towards zero
  ene += -0.5*std::log(inv_sss) - 0.5*narg*std::log(inv_s2*2.*M_PI);;
  return kbt_ * ene;
}

double Metainference::getEnergyGJE(const vector<double> &mean, const vector<double> &sigma, const double scale, const double modifier)
{
  const double scale2 = scale * scale;
  const double mu2 = modifier*modifier;
  const double mul_sm2 = mu2*scale2;

  double ene = 0.0;
  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(ene)
  { 
    #pragma omp for reduction( + : ene)
    for(unsigned i=0;i<narg;++i){
      const double sigma2 = sigma[i] * sigma[i];
      const double sigma_mean2 = sigma_mean_[i] * sigma_mean_[i];
      const double sss = sigma2 + mul_sm2*sigma_mean2;
      const double ss = sigma2 + mu2*sigma_mean2;
      double dev = scale*mean[i]-parameters[i];
      // deviation + normalisation + jeffrey
      // ene += 0.5*dev*dev/ss + 0.5*std::log(sss*2.*M_PI) + 0.5*std::log(ss);
      // scale doens't enter in the jeffrey, because we don't want the scale to be biased towards zero
      // this is equivalent to (but with less log to calculate)
      ene += 0.5*dev*dev/ss + 0.5*std::log(2.*M_PI*sss*ss);
    }
  }
  return kbt_ * ene;
}

void Metainference::doMonteCarlo(const vector<double> &mean_, const double modifier)
{
  // calculate old energy with the updated coordinates
  double old_energy;
  switch(noise_type_) {
    case GAUSS:
      old_energy = getEnergyGJ(mean_,sigma_,scale_,modifier);
      break;
    case MGAUSS:
      old_energy = getEnergyGJE(mean_,sigma_,scale_,modifier);
      break;
    case OUTLIERS:
      old_energy = getEnergySP(mean_,sigma_,scale_,modifier);
      break;
    case MOUTLIERS:
      old_energy = getEnergySPE(mean_,sigma_,scale_,modifier);
      break;
  }

  // Create vector of random sigma indices
  vector<unsigned> indices;
  if (MCchunksize_ > 0) {
      for (unsigned j=0; j<sigma_.size(); j++) {
          indices.push_back(j);
      }
      random[2].Shuffle(indices);
  }
  bool breaknow = false;

  // cycle on MC steps 
  for(unsigned i=0;i<MCsteps_;++i){
 
    // propose move for scale
    double new_scale = scale_;
    if(doscale_) {
      MCtrial_++;
      if(scale_prior_==SC_FLAT) {
        const double r1 = random[1].Gaussian();
        const double ds1 = Dscale_*r1;
        new_scale += ds1;
        // check boundaries
        if(new_scale > scale_max_){new_scale = 2.0 * scale_max_ - new_scale;}
        if(new_scale < scale_min_){new_scale = 2.0 * scale_min_ - new_scale;}
      } else {
        const double r1 = random[1].Gaussian();
        const double ds1 = 0.5*(scale_mu_-new_scale)+scale_sigma_*exp(1)/M_PI*r1;
        new_scale += ds1;
      }

      // calculate new energy
      double new_energy;
      switch(noise_type_) {
        case GAUSS:
          new_energy = getEnergyGJ(mean_,sigma_,new_scale,modifier);
          break;
        case MGAUSS:
          new_energy = getEnergyGJE(mean_,sigma_,new_scale,modifier);
          break;
        case OUTLIERS:
          new_energy = getEnergySP(mean_,sigma_,new_scale,modifier);
          break;
        case MOUTLIERS:
          new_energy = getEnergySPE(mean_,sigma_,new_scale,modifier);
          break;
      }
      // for the scale we need to consider the total energy
      vector<double> totenergies(2);
      if(master) {
        totenergies[0] = old_energy;
        totenergies[1] = new_energy;
        if(nrep_>1) multi_sim_comm.Sum(totenergies);
      } else {
        totenergies[0] = 0;
        totenergies[1] = 0;
      }
      comm.Sum(totenergies);

      // accept or reject
      const double delta = ( totenergies[1] - totenergies[0] ) / kbt_;
      // if delta is negative always accept move
      if( delta <= 0.0 ){
        old_energy = new_energy;
        scale_ = new_scale;
        MCaccept_++;
      // otherwise extract random number
      } else {
        double s = random[1].RandU01();
        if( s < exp(-delta) ){
          old_energy = new_energy;
          scale_ = new_scale;
          MCaccept_++;
        }
      }
    }
  
    // propose move for sigma
    MCtrial_++;

    vector<double> new_sigma(sigma_.size());
    new_sigma = sigma_;

    // change MCchunksize_ sigmas
    if (MCchunksize_ > 0) {
        if ((MCchunksize_ * i) >= sigma_.size()) {
            // This means we are not moving any sigma, so we should break immediately
            breaknow = true;
        }

        // change random sigmas
        for(unsigned j=0;j<MCchunksize_;j++) {
            const unsigned shuffle_index = j + MCchunksize_ * i;
            if (shuffle_index >= sigma_.size()) {
                // Going any further will segfault but we should still evaluate the sigmas we changed
                break;
            }
            const unsigned index = indices[shuffle_index];
            const double r2 = random[0].Gaussian();
            const double ds2 = Dsigma_[index]*r2;
            new_sigma[index] = sigma_[index] + ds2;
            // check boundaries
            if(new_sigma[index] > sigma_max_[index]){new_sigma[index] = 2.0 * sigma_max_[index] - new_sigma[index];}
            if(new_sigma[index] < sigma_min_[index]){new_sigma[index] = 2.0 * sigma_min_[index] - new_sigma[index];}
        }
    } else {
        // change all sigmas
        for(unsigned j=0;j<sigma_.size();j++) {
            const double r2 = random[0].Gaussian();
            const double ds2 = Dsigma_[j]*r2;
            new_sigma[j] = sigma_[j] + ds2;
            // check boundaries
            if(new_sigma[j] > sigma_max_[j]){new_sigma[j] = 2.0 * sigma_max_[j] - new_sigma[j];}
            if(new_sigma[j] < sigma_min_[j]){new_sigma[j] = 2.0 * sigma_min_[j] - new_sigma[j];}
        }
    }

    if (breaknow) {
        // We didnt move any sigmas, so no sense in evaluating anything
        break;
    }

 
    // calculate new energy
    double new_energy;
    switch(noise_type_) {
      case GAUSS:
        new_energy = getEnergyGJ(mean_,new_sigma,scale_,modifier);
        break;
      case MGAUSS:
        new_energy = getEnergyGJE(mean_,new_sigma,scale_,modifier);
        break;
      case OUTLIERS:
        new_energy = getEnergySP(mean_,new_sigma,scale_,modifier);
        break;
      case MOUTLIERS:
        new_energy = getEnergySPE(mean_,new_sigma,scale_,modifier);
        break;
    }

    // accept or reject
    const double delta = ( new_energy - old_energy ) / kbt_;
    // if delta is negative always accept move
    if( delta <= 0.0 ){
      old_energy = new_energy;
      sigma_ = new_sigma;
      MCaccept_++;
    // otherwise extract random number
    } else {
      const double s = random[0].RandU01();
      if( s < exp(-delta) ){
        old_energy = new_energy;
        sigma_ = new_sigma;
        MCaccept_++;
      }
    }
 
  }
  /* save the result of the sampling */
  double accept = static_cast<double>(MCaccept_) / static_cast<double>(MCtrial_);
  valueAccept->set(accept);
  if(doscale_) valueScale->set(scale_);
  for(unsigned i=0; i<sigma_.size(); i++) valueSigma[i]->set(sigma_[i]);
}

/* 
   In the following energy-force functions we don't add the normalisation and the jeffreys priors
   because they are not needed for the forces, the correct MetaInference energy is the one calculated
   in the Monte-Carlo 
*/

double Metainference::getEnergyForceSP(const vector<double> &mean, const double fact, const double modifier)
{
  const double sm2 = modifier*modifier*sigma_mean_[0]*sigma_mean_[0]; 
  const double ss = sigma_[0]*sigma_[0] + sm2;
  vector<double> f(narg+1,0);
  
  if(master){
    double omp_ene=0.;
    #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(omp_ene)
    { 
      #pragma omp for reduction( + : omp_ene)
      for(unsigned i=0;i<narg;++i){
        const double dev = scale_*mean[i]-parameters[i]; 
        const double a2 = 0.5*dev*dev + ss;
        const double t = exp(-a2/sm2);
        const double dt = 1./t;
        const double it = 1./(1.-t);
        const double dit = 1./(1.-dt);
        omp_ene += std::log(2.*a2*it);
        f[i] = -scale_*dev*(dit/sm2 + 1./a2);
      }
    }
    f[narg] = omp_ene;
    // collect contribution to forces and energy from other replicas
    if(nrep_>1) multi_sim_comm.Sum(&f[0],narg+1);
  }
  // intra-replica summation
  comm.Sum(&f[0],narg+1);

  const double ene = f[narg];
  double w_tmp = 0.;
  for(unsigned i=0; i<narg; ++i) {
    setOutputForce(i, kbt_ * fact * f[i]);
    w_tmp += fact*(getArgument(i) - mean[i])*f[i];
  }

  if(do_reweight) {
    setOutputForce(narg, -w_tmp);
    getPntrToComponent("MetaDf")->set(w_tmp);
    getPntrToComponent("weight")->set(fact);
  }

  return kbt_*ene;
}

double Metainference::getEnergyForceSPE(const vector<double> &mean, const double fact, const double modifier)
{
  const double mul_sm2 = modifier*modifier;
  vector<double> f(narg+1,0);

  if(master){
    double omp_ene = 0;
    #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(omp_ene)
    { 
      #pragma omp for reduction( + : omp_ene)
      for(unsigned i=0;i<narg;++i){
        const double sm2 = mul_sm2*sigma_mean_[i]*sigma_mean_[i]; 
        const double ss  = sigma_[i]*sigma_[i] + sm2;
        const double dev = scale_*mean[i]-parameters[i]; 
        const double a2  = 0.5*dev*dev + ss;
        const double t   = exp(-a2/sm2);
        const double dt  = 1./t;
        const double it  = 1./(1.-t);
        const double dit = 1./(1.-dt);
        omp_ene += std::log(2.*a2*it);
        f[i] = -scale_*dev*(dit/sm2 + 1./a2);
      }
    }
    f[narg] = omp_ene;
    // collect contribution to forces and energy from other replicas
    if(nrep_>1) multi_sim_comm.Sum(&f[0],narg+1);
  }
  comm.Sum(&f[0],narg+1);

  const double ene = f[narg];
  double w_tmp = 0.;
  for(unsigned i=0; i<narg; ++i) {
    setOutputForce(i, kbt_ * fact * f[i]);
    w_tmp += fact*(getArgument(i) - mean[i])*f[i];
  }

  if(do_reweight) {
    setOutputForce(narg, -w_tmp);
    getPntrToComponent("MetaDf")->set(w_tmp);
    getPntrToComponent("weight")->set(fact);
  }

  return kbt_*ene;
}

double Metainference::getEnergyForceGJ(const vector<double> &mean, const double fact, const double modifier)
{
  double inv_s2=0.;
  if(master) {
    const double ss = sigma_[0]*sigma_[0] + modifier*modifier*sigma_mean_[0]*sigma_mean_[0];
    inv_s2 = 1.0/ss;
    if(nrep_>1) multi_sim_comm.Sum(inv_s2);
  } 
  comm.Sum(inv_s2);  

  double ene = 0.0;
  double w_tmp = 0.;
  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(ene,w_tmp)
  { 
    #pragma omp for reduction( + : ene,w_tmp)
    for(unsigned i=0;i<narg;++i){
      const double dev = scale_*mean[i]-parameters[i];
      const double mult = fact*dev*scale_*inv_s2;
      ene += 0.5*dev*dev*inv_s2;
      setOutputForce(i, -kbt_*mult);
      w_tmp += (getArgument(i)-mean[i])*mult;
    }
  }
  if(do_reweight) {
    setOutputForce(narg, -w_tmp);
    getPntrToComponent("MetaDf")->set(w_tmp);
    getPntrToComponent("weight")->set(fact);
  }

  return kbt_*ene;
}

double Metainference::getEnergyForceGJE(const vector<double> &mean, const double fact, const double modifier)
{
  const double sm_m2 = modifier*modifier;
  vector<double> inv_s2(sigma_.size());

  if(master) {
    for(unsigned i=0;i<sigma_.size(); ++i) {
      const double ss = sigma_[i]*sigma_[i] + sm_m2*sigma_mean_[i]*sigma_mean_[i];
      inv_s2[i] = 1.0/ss;
    }
    if(nrep_>1) multi_sim_comm.Sum(&inv_s2[0],sigma_.size());
  } else { 
    for(unsigned i=0;i<sigma_.size(); ++i) inv_s2[i] = 0.;
  }
  comm.Sum(&inv_s2[0],sigma_.size());  
  
  double ene = 0.0;
  double w_tmp = 0.;

  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(ene,w_tmp)
  { 
    #pragma omp for reduction( + : ene,w_tmp)
    for(unsigned i=0;i<narg;++i){
      const double dev  = scale_*mean[i]-parameters[i];
      const double mult = fact*dev*scale_*inv_s2[i];
      ene += 0.5*dev*dev*inv_s2[i];
      setOutputForce(i, -kbt_*mult);
      w_tmp += (getArgument(i)-mean[i])*mult;
    }
  }

  if(do_reweight) {
    setOutputForce(narg, -w_tmp);
    getPntrToComponent("MetaDf")->set(w_tmp);
    getPntrToComponent("weight")->set(fact);
  }

  return kbt_*ene;
}

void Metainference::calculate()
{
  double norm = 0.0;
  double fact = 0.0;
  double idof = 1.0;
  const double dnrep = static_cast<double>(nrep_);

  if(do_reweight){
    // calculate the weights either from BIAS 
    vector<double> bias(nrep_,0);
    if(master){
      bias[replica_] = getArgument(narg); 
      if(nrep_>1) multi_sim_comm.Sum(&bias[0], nrep_);  
    }
    comm.Sum(&bias[0], nrep_);
    const double maxbias = *(std::max_element(bias.begin(), bias.end()));
    double n2=0.;
    for(unsigned i=0; i<nrep_; ++i){
      bias[i] = exp((bias[i]-maxbias)/kbt_); 
      norm += bias[i];
      n2 += bias[i]*bias[i];
    }
    fact = bias[replica_]/norm;
    idof = sqrt(1./(1. - n2/(norm*norm)));
  } else {
    // or arithmetic ones
    norm = dnrep; 
    fact = 1.0/norm; 
  }

  vector<double> mean(narg,0);
  // calculate the mean 
  if(master) {
    for(unsigned i=0;i<narg;++i) mean[i] = fact*getArgument(i); 
    if(nrep_>1) multi_sim_comm.Sum(&mean[0], narg);
  }
  comm.Sum(&mean[0], narg);

  if(do_optsigmamean_>0) {
    /* this is SIGMA_MEAN before the corrections due to the #DOF and the SCALING */
    vector<double> v_moment(narg,0);
    if(master) {
      for(unsigned i=0;i<narg;++i) { 
        double tmp  = getArgument(i)-mean[i];
        v_moment[i] = fact*tmp*tmp;
      }
      if(nrep_>1) multi_sim_comm.Sum(&v_moment[0], narg);
    }
    comm.Sum(&v_moment[0], narg);
 
    const double sq_dnrep = sqrt(dnrep);
    bool sm_update = false;
    for(unsigned i=0;i<narg;++i) { 
      if(v_moment[i]>variance_[i]) {
        sm_update = true;
        variance_[i] = v_moment[i];
        if(noise_type_==MGAUSS||noise_type_==MOUTLIERS) {
          double s_v = sqrt(variance_[i]);
          if(sigma_max_[i] < s_v) {
            Dsigma_[i] *= s_v/sigma_max_[i]; 
            sigma_max_[i] = s_v;
          }
          sigma_mean_[i] = s_v/sq_dnrep;
          valueSigmaMean[i]->set(sigma_mean_[i]);
        }
      }
    }
    if(sm_update&&(noise_type_==GAUSS||noise_type_==OUTLIERS)) {
      double s_v = sqrt(*max_element(variance_.begin(), variance_.end()));
      if(sigma_max_[0] < s_v) {
        Dsigma_[0] *= s_v/sigma_max_[0]; 
        sigma_max_[0] = s_v;
      }
      sigma_mean_[0] = s_v/sq_dnrep;
      valueSigmaMean[0]->set(sigma_mean_[0]);
    }
  }

  // correct sigma_mean for the weighted average effect
  double sigma_mean_modifier = idof;

  /* MONTE CARLO */
  const long int step = getStep();
  if(step%MCstride_==0&&!getExchangeStep()) doMonteCarlo(mean, sigma_mean_modifier);

  /* fix sigma_mean_ for the scaling factor */
  sigma_mean_modifier *= scale_;
  /* fix sigma_mean_ for the effect of large forces */
  if(do_optsigmamean_==2) sigma_mean_modifier *= sm_mod_;

  valueRSigmaMean->set(sigma_mean_modifier);
  
  // calculate bias and forces
  double ene = 0; 
  switch(noise_type_) {
    case GAUSS:
      ene = getEnergyForceGJ(mean, fact, sigma_mean_modifier);
      break;
    case MGAUSS:
      ene = getEnergyForceGJE(mean, fact, sigma_mean_modifier);
      break;
    case OUTLIERS:
      ene = getEnergyForceSP(mean, fact, sigma_mean_modifier);
      break;
    case MOUTLIERS:
      ene = getEnergyForceSPE(mean, fact, sigma_mean_modifier);
      break;
  }

  // set value of the bias
  setBias(ene);
}

void Metainference::writeStatus()
{
  sfile_.rewind();
  sfile_.printField("time",getTimeStep()*getStep());
  for(unsigned i=0;i<variance_.size();++i) {
    std::string msg;
    Tools::convert(i,msg);
    sfile_.printField("variance_"+msg,variance_[i]);
  }
  for(unsigned i=0;i<sigma_.size();++i) {
    std::string msg;
    Tools::convert(i,msg);
    sfile_.printField("sigma_"+msg,sigma_[i]);
  }
  if(doscale_) {
    sfile_.printField("scale0_",scale_);
  }
  if(do_optsigmamean_==2) {
    sfile_.printField("sigma_mean_mod0",sm_mod_);
  }
  sfile_.printField();
  sfile_.flush();
}

void Metainference::update() {
  if(do_optsigmamean_==2) {
    const double EPS = 0.1;
    // Get max force of whole system
    vector<Vector> md_forces;
    vector<Vector> plumed_forces;
    vector<double> masses;

    atoms.getLocalMDForces(md_forces);
    atoms.getLocalForces(plumed_forces);
    atoms.getLocalMasses(masses);

    vector<double> allforces_md;
    allforces_md.reserve(md_forces.size());

    for(unsigned i = 0; i < plumed_forces.size(); ++i) {
      const double pf2 = plumed_forces[i].modulo2();
      // we are only interested in forces plumed has an effect on
      if(pf2 > EPS && masses[i] > EPS ) {
        const double invm2 = 1./(masses[i]*masses[i]);
        allforces_md.push_back(md_forces[i].modulo2()*invm2);
      }
    }

    vector<double> fmax_tmp_md(comm.Get_size(),0);
    vector<double> nfmax_md(nrep_, 0);

    /* each local thread should look for the maximum force and number of violations */
    if(allforces_md.size()>0) {
      fmax_tmp_md[comm.Get_rank()] = *max_element(allforces_md.begin(), allforces_md.end());
      for(unsigned i = 0; i < allforces_md.size(); ++i) {
        if(allforces_md[i] > max_force_) {
          nfmax_md[replica_] += allforces_md[i]/max_force_;
        }
      }
    }
    // the largest forces are shared among the local threads but not over the replicas 
    comm.Sum(fmax_tmp_md);
    // these are the largest forces for a specific replica
    const double fmax_md = *max_element(fmax_tmp_md.begin(), fmax_tmp_md.end());

    // the number of violations is summed up over the local thread and over the replicas 
    comm.Sum(nfmax_md);
    if(master && nrep_ > 1) multi_sim_comm.Sum(nfmax_md);
    comm.Bcast(&nfmax_md[0], nrep_, 0);
    
    const double nnfmax_md = (*max_element(nfmax_md.begin(), nfmax_md.end()))*nrep_;

    if( nnfmax_md == 0) {
      sm_mod_ -= Dsm_mod_ * 0.01 * std::log(sm_mod_/sm_mod_min_);
      if(sm_mod_<sm_mod_min_) sm_mod_=sm_mod_min_;
    } else {
      const double sm_mod_new = sm_mod_ + Dsm_mod_ * std::log(nnfmax_md+1.);
      if(sm_mod_new > sm_mod_max_) {
        sm_mod_ = sm_mod_max_;
      } else {
        sm_mod_ = sm_mod_new;
      }
    }

    valueSMmod->set(sm_mod_);
    valueMaxForceMD->set(sqrt(fmax_md));
  }

  // write status file
  const long int step = getStep();
  if(write_stride_>0&& (step%write_stride_==0 || getCPT()) ) writeStatus();
}

}
}

