/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2023 The plumed team
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
#include "MetainferenceBase.h"
#include "tools/File.h"
#include <chrono>
#include <numeric>

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

namespace PLMD {
namespace isdb {

void MetainferenceBase::registerKeywords( Keywords& keys ) {
  Action::registerKeywords(keys);
  ActionAtomistic::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  componentsAreNotOptional(keys);
  keys.use("ARG");
  keys.addFlag("DOSCORE",false,"activate metainference");
  keys.addFlag("NOENSEMBLE",false,"don't perform any replica-averaging");
  keys.addFlag("REWEIGHT",false,"simple REWEIGHT using the ARG as energy");
  keys.add("optional","AVERAGING", "Stride for calculation of averaged weights and sigma_mean");
  keys.add("compulsory","NOISETYPE","MGAUSS","functional form of the noise (GAUSS,MGAUSS,OUTLIERS,MOUTLIERS,GENERIC)");
  keys.add("compulsory","LIKELIHOOD","GAUSS","the likelihood for the GENERIC metainference model, GAUSS or LOGN");
  keys.add("compulsory","DFTILDE","0.1","fraction of sigma_mean used to evolve ftilde");
  keys.addFlag("SCALEDATA",false,"Set to TRUE if you want to sample a scaling factor common to all values and replicas");
  keys.add("compulsory","SCALE0","1.0","initial value of the scaling factor");
  keys.add("compulsory","SCALE_PRIOR","FLAT","either FLAT or GAUSSIAN");
  keys.add("optional","SCALE_MIN","minimum value of the scaling factor");
  keys.add("optional","SCALE_MAX","maximum value of the scaling factor");
  keys.add("optional","DSCALE","maximum MC move of the scaling factor");
  keys.addFlag("ADDOFFSET",false,"Set to TRUE if you want to sample an offset common to all values and replicas");
  keys.add("compulsory","OFFSET0","0.0","initial value of the offset");
  keys.add("compulsory","OFFSET_PRIOR","FLAT","either FLAT or GAUSSIAN");
  keys.add("optional","OFFSET_MIN","minimum value of the offset");
  keys.add("optional","OFFSET_MAX","maximum value of the offset");
  keys.add("optional","DOFFSET","maximum MC move of the offset");
  keys.add("optional","REGRES_ZERO","stride for regression with zero offset");
  keys.add("compulsory","SIGMA0","1.0","initial value of the uncertainty parameter");
  keys.add("compulsory","SIGMA_MIN","0.0","minimum value of the uncertainty parameter");
  keys.add("compulsory","SIGMA_MAX","10.","maximum value of the uncertainty parameter");
  keys.add("optional","DSIGMA","maximum MC move of the uncertainty parameter");
  keys.add("compulsory","OPTSIGMAMEAN","NONE","Set to NONE/SEM to manually set sigma mean, or to estimate it on the fly");
  keys.add("optional","SIGMA_MEAN0","starting value for the uncertainty in the mean estimate");
  keys.add("optional","SIGMA_MAX_STEPS", "Number of steps used to optimise SIGMA_MAX, before that the SIGMA_MAX value is used");
  keys.add("optional","TEMP","the system temperature - this is only needed if code doesn't pass the temperature to plumed");
  keys.add("optional","MC_STEPS","number of MC steps");
  keys.add("optional","MC_CHUNKSIZE","MC chunksize");
  keys.add("optional","STATUS_FILE","write a file with all the data useful for restart/continuation of Metainference");
  keys.add("compulsory","WRITE_STRIDE","10000","write the status to a file every N steps, this can be used for restart/continuation");
  keys.add("optional","SELECTOR","name of selector");
  keys.add("optional","NSELECT","range of values for selector [0, N-1]");
  keys.use("RESTART");
  keys.addOutputComponent("score",        "default",      "the Metainference score");
  keys.addOutputComponent("sigma",        "default",      "uncertainty parameter");
  keys.addOutputComponent("sigmaMean",    "default",      "uncertainty in the mean estimate");
  keys.addOutputComponent("neff",         "default",      "effective number of replicas");
  keys.addOutputComponent("acceptSigma",  "default",      "MC acceptance for sigma values");
  keys.addOutputComponent("acceptScale",  "SCALEDATA",    "MC acceptance for scale value");
  keys.addOutputComponent("acceptFT",     "GENERIC",      "MC acceptance for general metainference f tilde value");
  keys.addOutputComponent("weight",       "REWEIGHT",     "weights of the weighted average");
  keys.addOutputComponent("biasDer",      "REWEIGHT",     "derivatives with respect to the bias");
  keys.addOutputComponent("scale",        "SCALEDATA",    "scale parameter");
  keys.addOutputComponent("offset",       "ADDOFFSET",    "offset parameter");
  keys.addOutputComponent("ftilde",       "GENERIC",      "ensemble average estimator");
}

MetainferenceBase::MetainferenceBase(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionWithArguments(ao),
  ActionWithValue(ao),
  doscore_(false),
  write_stride_(0),
  narg(0),
  doscale_(false),
  scale_(1.),
  scale_mu_(0),
  scale_min_(1),
  scale_max_(-1),
  Dscale_(-1),
  dooffset_(false),
  offset_(0.),
  offset_mu_(0),
  offset_min_(1),
  offset_max_(-1),
  Doffset_(-1),
  doregres_zero_(false),
  nregres_zero_(0),
  Dftilde_(0.1),
  random(3),
  MCsteps_(1),
  MCaccept_(0),
  MCacceptScale_(0),
  MCacceptFT_(0),
  MCtrial_(0),
  MCchunksize_(0),
  firstTime(true),
  do_reweight_(false),
  do_optsigmamean_(0),
  nsel_(1),
  iselect(0),
  optsigmamean_stride_(0),
  N_optimized_step_(0),
  optimized_step_(0),
  sigmamax_opt_done_(false),
  decay_w_(1.)
{
  parseFlag("DOSCORE", doscore_);

  bool noensemble = false;
  parseFlag("NOENSEMBLE", noensemble);

  // set up replica stuff
  master = (comm.Get_rank()==0);
  if(master) {
    nrep_    = multi_sim_comm.Get_size();
    replica_ = multi_sim_comm.Get_rank();
    if(noensemble) nrep_ = 1;
  } else {
    nrep_    = 0;
    replica_ = 0;
  }
  comm.Sum(&nrep_,1);
  comm.Sum(&replica_,1);

  parse("SELECTOR", selector_);
  parse("NSELECT", nsel_);
  // do checks
  if(selector_.length()>0 && nsel_<=1) error("With SELECTOR active, NSELECT must be greater than 1");
  if(selector_.length()==0 && nsel_>1) error("With NSELECT greater than 1, you must specify SELECTOR");

  // initialise firstTimeW
  firstTimeW.resize(nsel_, true);

  // reweight implies a different number of arguments (the latest one must always be the bias)
  parseFlag("REWEIGHT", do_reweight_);
  if(do_reweight_&&getNumberOfArguments()!=1) error("To REWEIGHT one must provide one single bias as an argument");
  if(do_reweight_&&nrep_<2) error("REWEIGHT can only be used in parallel with 2 or more replicas");
  if(!getRestart()) average_weights_.resize(nsel_, std::vector<double> (nrep_, 1./static_cast<double>(nrep_)));
  else average_weights_.resize(nsel_, std::vector<double> (nrep_, 0.));

  unsigned averaging=0;
  parse("AVERAGING", averaging);
  if(averaging>0) {
    decay_w_ = 1./static_cast<double> (averaging);
    optsigmamean_stride_ = averaging;
  }

  std::string stringa_noise;
  parse("NOISETYPE",stringa_noise);
  if(stringa_noise=="GAUSS")           noise_type_ = GAUSS;
  else if(stringa_noise=="MGAUSS")     noise_type_ = MGAUSS;
  else if(stringa_noise=="OUTLIERS")   noise_type_ = OUTLIERS;
  else if(stringa_noise=="MOUTLIERS")  noise_type_ = MOUTLIERS;
  else if(stringa_noise=="GENERIC")    noise_type_ = GENERIC;
  else error("Unknown noise type!");

  if(noise_type_== GENERIC) {
    std::string stringa_like;
    parse("LIKELIHOOD",stringa_like);
    if(stringa_like=="GAUSS") gen_likelihood_ = LIKE_GAUSS;
    else if(stringa_like=="LOGN") gen_likelihood_ = LIKE_LOGN;
    else error("Unknown likelihood type!");

    parse("DFTILDE",Dftilde_);
  }

  parse("WRITE_STRIDE",write_stride_);
  parse("STATUS_FILE",status_file_name_);
  if(status_file_name_=="") status_file_name_ = "MISTATUS"+getLabel();
  else                      status_file_name_ = status_file_name_+getLabel();

  std::string stringa_optsigma;
  parse("OPTSIGMAMEAN", stringa_optsigma);
  if(stringa_optsigma=="NONE")      do_optsigmamean_=0;
  else if(stringa_optsigma=="SEM")  do_optsigmamean_=1;
  else if(stringa_optsigma=="SEM_MAX")  do_optsigmamean_=2;

  unsigned aver_max_steps=0;
  parse("SIGMA_MAX_STEPS", aver_max_steps);
  if(aver_max_steps==0&&do_optsigmamean_==2) aver_max_steps=averaging*2000;
  if(aver_max_steps>0&&do_optsigmamean_<2) error("SIGMA_MAX_STEPS can only be used together with OPTSIGMAMEAN=SEM_MAX");
  if(aver_max_steps>0&&do_optsigmamean_==2) N_optimized_step_=aver_max_steps;
  if(aver_max_steps>0&&aver_max_steps<averaging) error("SIGMA_MAX_STEPS must be greater than AVERAGING");

  std::vector<double> read_sigma_mean_;
  parseVector("SIGMA_MEAN0",read_sigma_mean_);
  if(do_optsigmamean_==0 && read_sigma_mean_.size()==0 && !getRestart() && doscore_)
    error("If you don't use OPTSIGMAMEAN and you are not RESTARTING then you MUST SET SIGMA_MEAN0");

  if(noise_type_==MGAUSS||noise_type_==MOUTLIERS||noise_type_==GENERIC) {
    if(read_sigma_mean_.size()>0) {
      sigma_mean2_.resize(read_sigma_mean_.size());
      for(unsigned i=0; i<read_sigma_mean_.size(); i++) sigma_mean2_[i]=read_sigma_mean_[i]*read_sigma_mean_[i];
    } else {
      sigma_mean2_.resize(1,0.000001);
    }
  } else {
    if(read_sigma_mean_.size()==1) {
      sigma_mean2_.resize(1, read_sigma_mean_[0]*read_sigma_mean_[0]);
    } else if(read_sigma_mean_.size()==0) {
      sigma_mean2_.resize(1, 0.000001);
    } else {
      error("If you want to use more than one SIGMA_MEAN0 you should use NOISETYPE=MGAUSS|MOUTLIERS");
    }
  }

  parseFlag("SCALEDATA", doscale_);
  if(doscale_) {
    std::string stringa_noise;
    parse("SCALE_PRIOR",stringa_noise);
    if(stringa_noise=="GAUSSIAN")  scale_prior_ = SC_GAUSS;
    else if(stringa_noise=="FLAT") scale_prior_ = SC_FLAT;
    else error("Unknown SCALE_PRIOR type!");
    parse("SCALE0",scale_);
    parse("DSCALE",Dscale_);
    if(Dscale_<0.) error("DSCALE must be set when using SCALEDATA");
    if(scale_prior_==SC_GAUSS) {
      scale_mu_=scale_;
    } else {
      parse("SCALE_MIN",scale_min_);
      parse("SCALE_MAX",scale_max_);
      if(scale_max_<scale_min_) error("SCALE_MAX and SCALE_MIN must be set when using SCALE_PRIOR=FLAT");
    }
  }

  parseFlag("ADDOFFSET", dooffset_);
  if(dooffset_) {
    std::string stringa_noise;
    parse("OFFSET_PRIOR",stringa_noise);
    if(stringa_noise=="GAUSSIAN")  offset_prior_ = SC_GAUSS;
    else if(stringa_noise=="FLAT") offset_prior_ = SC_FLAT;
    else error("Unknown OFFSET_PRIOR type!");
    parse("OFFSET0",offset_);
    parse("DOFFSET",Doffset_);
    if(offset_prior_==SC_GAUSS) {
      offset_mu_=offset_;
      if(Doffset_<0.) error("DOFFSET must be set when using OFFSET_PRIOR=GAUSS");
    } else {
      parse("OFFSET_MIN",offset_min_);
      parse("OFFSET_MAX",offset_max_);
      if(Doffset_<0) Doffset_ = 0.05*(offset_max_ - offset_min_);
      if(offset_max_<offset_min_) error("OFFSET_MAX and OFFSET_MIN must be set when using OFFSET_PRIOR=FLAT");
    }
  }

  // regression with zero intercept
  parse("REGRES_ZERO", nregres_zero_);
  if(nregres_zero_>0) {
    // set flag
    doregres_zero_=true;
    // check if already sampling scale and offset
    if(doscale_)  error("REGRES_ZERO and SCALEDATA are mutually exclusive");
    if(dooffset_) error("REGRES_ZERO and ADDOFFSET are mutually exclusive");
  }

  std::vector<double> readsigma;
  parseVector("SIGMA0",readsigma);
  if((noise_type_!=MGAUSS&&noise_type_!=MOUTLIERS&&noise_type_!=GENERIC)&&readsigma.size()>1)
    error("If you want to use more than one SIGMA you should use NOISETYPE=MGAUSS|MOUTLIERS|GENERIC");
  if(noise_type_==MGAUSS||noise_type_==MOUTLIERS||noise_type_==GENERIC) {
    sigma_.resize(readsigma.size());
    sigma_=readsigma;
  } else sigma_.resize(1, readsigma[0]);

  std::vector<double> readsigma_min;
  parseVector("SIGMA_MIN",readsigma_min);
  if((noise_type_!=MGAUSS&&noise_type_!=MOUTLIERS&&noise_type_!=GENERIC)&&readsigma_min.size()>1)
    error("If you want to use more than one SIGMA you should use NOISETYPE=MGAUSS|MOUTLIERS|GENERIC");
  if(noise_type_==MGAUSS||noise_type_==MOUTLIERS||noise_type_==GENERIC) {
    sigma_min_.resize(readsigma_min.size());
    sigma_min_=readsigma_min;
  } else sigma_min_.resize(1, readsigma_min[0]);

  std::vector<double> readsigma_max;
  parseVector("SIGMA_MAX",readsigma_max);
  if((noise_type_!=MGAUSS&&noise_type_!=MOUTLIERS&&noise_type_!=GENERIC)&&readsigma_max.size()>1)
    error("If you want to use more than one SIGMA you should use NOISETYPE=MGAUSS|MOUTLIERS|GENERIC");
  if(noise_type_==MGAUSS||noise_type_==MOUTLIERS||noise_type_==GENERIC) {
    sigma_max_.resize(readsigma_max.size());
    sigma_max_=readsigma_max;
  } else sigma_max_.resize(1, readsigma_max[0]);

  if(sigma_max_.size()!=sigma_min_.size()) error("The number of values for SIGMA_MIN and SIGMA_MAX must be the same");

  std::vector<double> read_dsigma;
  parseVector("DSIGMA",read_dsigma);
  if((noise_type_!=MGAUSS&&noise_type_!=MOUTLIERS&&noise_type_!=GENERIC)&&readsigma_max.size()>1)
    error("If you want to use more than one SIGMA you should use NOISETYPE=MGAUSS|MOUTLIERS|GENERIC");
  if(read_dsigma.size()>0) {
    Dsigma_.resize(read_dsigma.size());
    Dsigma_=read_dsigma;
  } else {
    Dsigma_.resize(sigma_max_.size(), -1.);
    /* in this case Dsigma is initialised after reading the restart file if present */
  }

  // monte carlo stuff
  parse("MC_STEPS",MCsteps_);
  parse("MC_CHUNKSIZE", MCchunksize_);
  // get temperature
  double temp=0.0;
  parse("TEMP",temp);
  if(temp>0.0) kbt_=plumed.getAtoms().getKBoltzmann()*temp;
  else kbt_=plumed.getAtoms().getKbT();
  if(kbt_==0.0&&doscore_) error("Unless the MD engine passes the temperature to plumed, you must specify it using TEMP");

  // initialize random seed
  unsigned iseed;
  if(master) {
    auto ts = std::chrono::time_point_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now()).time_since_epoch().count();
    iseed = static_cast<unsigned>(ts)+replica_;
  } else {
    iseed = 0;
  }
  comm.Sum(&iseed, 1);
  // this is used for ftilde and sigma both the move and the acceptance
  // this is different for each replica
  random[0].setSeed(-iseed);
  if(doscale_||dooffset_) {
    // in this case we want the same seed everywhere
    auto ts = std::chrono::time_point_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now()).time_since_epoch().count();
    iseed = static_cast<unsigned>(ts);
    if(master&&nrep_>1) multi_sim_comm.Bcast(iseed,0);
    comm.Bcast(iseed,0);
    // this is used for scale and offset sampling and acceptance
    random[1].setSeed(-iseed);
  }
  // this is used for random chunk of sigmas, and it is different for each replica
  if(master) {
    auto ts = std::chrono::time_point_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now()).time_since_epoch().count();
    iseed = static_cast<unsigned>(ts)+replica_;
  } else {
    iseed = 0;
  }
  comm.Sum(&iseed, 1);
  random[2].setSeed(-iseed);

  // outfile stuff
  if(write_stride_>0&&doscore_) {
    sfile_.link(*this);
    sfile_.open(status_file_name_);
  }

}

MetainferenceBase::~MetainferenceBase()
{
  if(sfile_.isOpen()) sfile_.close();
}

void MetainferenceBase::Initialise(const unsigned input)
{
  setNarg(input);
  if(narg!=parameters.size()) {
    std::string num1; Tools::convert(parameters.size(),num1);
    std::string num2; Tools::convert(narg,num2);
    std::string msg = "The number of experimental values " + num1 +" must be the same of the calculated values " + num2;
    error(msg);
  }

  // resize std::vector for sigma_mean history
  sigma_mean2_last_.resize(nsel_);
  for(unsigned i=0; i<nsel_; i++) sigma_mean2_last_[i].resize(narg);
  if(noise_type_==MGAUSS||noise_type_==MOUTLIERS||noise_type_==GENERIC) {
    if(sigma_mean2_.size()==1) {
      double tmp = sigma_mean2_[0];
      sigma_mean2_.resize(narg, tmp);
    } else if(sigma_mean2_.size()>1&&sigma_mean2_.size()!=narg) {
      error("SIGMA_MEAN0 can accept either one single value or as many values as the number of arguments (with NOISETYPE=MGAUSS|MOUTLIERS|GENERIC)");
    }
    // set the initial value for the history
    for(unsigned i=0; i<nsel_; i++) for(unsigned j=0; j<narg; j++) sigma_mean2_last_[i][j].push_back(sigma_mean2_[j]);
  } else {
    // set the initial value for the history
    for(unsigned i=0; i<nsel_; i++) for(unsigned j=0; j<narg; j++) sigma_mean2_last_[i][j].push_back(sigma_mean2_[0]);
  }

  // set sigma_bias
  if(noise_type_==MGAUSS||noise_type_==MOUTLIERS||noise_type_==GENERIC) {
    if(sigma_.size()==1) {
      double tmp = sigma_[0];
      sigma_.resize(narg, tmp);
    } else if(sigma_.size()>1&&sigma_.size()!=narg) {
      error("SIGMA0 can accept either one single value or as many values as the number of arguments (with NOISETYPE=MGAUSS|MOUTLIERS|GENERIC)");
    }
    if(sigma_min_.size()==1) {
      double tmp = sigma_min_[0];
      sigma_min_.resize(narg, tmp);
    } else if(sigma_min_.size()>1&&sigma_min_.size()!=narg) {
      error("SIGMA_MIN can accept either one single value or as many values as the number of arguments (with NOISETYPE=MGAUSS|MOUTLIERS|GENERIC)");
    }
    if(sigma_max_.size()==1) {
      double tmp = sigma_max_[0];
      sigma_max_.resize(narg, tmp);
    } else if(sigma_max_.size()>1&&sigma_max_.size()!=narg) {
      error("SIGMA_MAX can accept either one single value or as many values as the number of arguments (with NOISETYPE=MGAUSS|MOUTLIERS|GENERIC)");
    }
    if(Dsigma_.size()==1) {
      double tmp = Dsigma_[0];
      Dsigma_.resize(narg, tmp);
    } else if(Dsigma_.size()>1&&Dsigma_.size()!=narg) {
      error("DSIGMA can accept either one single value or as many values as the number of arguments (with NOISETYPE=MGAUSS|MOUTLIERS|GENERIC)");
    }
  }

  sigma_max_est_.resize(sigma_max_.size(), 0.);

  IFile restart_sfile;
  restart_sfile.link(*this);
  if(getRestart()&&restart_sfile.FileExist(status_file_name_)) {
    firstTime = false;
    for(unsigned i=0; i<nsel_; i++) firstTimeW[i] = false;
    restart_sfile.open(status_file_name_);
    log.printf("  Restarting from %s\n", status_file_name_.c_str());
    double dummy;
    if(restart_sfile.scanField("time",dummy)) {
      // check for syncronisation
      std::vector<double> dummy_time(nrep_,0);
      if(master&&nrep_>1) {
        dummy_time[replica_] = dummy;
        multi_sim_comm.Sum(dummy_time);
      }
      comm.Sum(dummy_time);
      for(unsigned i=1; i<nrep_; i++) {
        std::string msg = "METAINFERENCE restart files " + status_file_name_ + "  are not in sync";
        if(dummy_time[i]!=dummy_time[0]) plumed_merror(msg);
      }
      // nsel
      for(unsigned i=0; i<sigma_mean2_last_.size(); i++) {
        std::string msg_i;
        Tools::convert(i,msg_i);
        // narg
        if(noise_type_==MGAUSS||noise_type_==MOUTLIERS||noise_type_==GENERIC) {
          for(unsigned j=0; j<narg; ++j) {
            std::string msg_j;
            Tools::convert(j,msg_j);
            std::string msg = msg_i+"_"+msg_j;
            double read_sm;
            restart_sfile.scanField("sigmaMean_"+msg,read_sm);
            sigma_mean2_last_[i][j][0] = read_sm*read_sm;
          }
        }
        if(noise_type_==GAUSS||noise_type_==OUTLIERS) {
          double read_sm;
          std::string msg_j;
          Tools::convert(0,msg_j);
          std::string msg = msg_i+"_"+msg_j;
          restart_sfile.scanField("sigmaMean_"+msg,read_sm);
          for(unsigned j=0; j<narg; j++) sigma_mean2_last_[i][j][0] = read_sm*read_sm;
        }
      }

      for(unsigned i=0; i<sigma_.size(); ++i) {
        std::string msg;
        Tools::convert(i,msg);
        restart_sfile.scanField("sigma_"+msg,sigma_[i]);
      }
      for(unsigned i=0; i<sigma_max_.size(); ++i) {
        std::string msg;
        Tools::convert(i,msg);
        restart_sfile.scanField("sigma_max_"+msg,sigma_max_[i]);
        sigmamax_opt_done_=true;
      }
      if(noise_type_==GENERIC) {
        for(unsigned i=0; i<ftilde_.size(); ++i) {
          std::string msg;
          Tools::convert(i,msg);
          restart_sfile.scanField("ftilde_"+msg,ftilde_[i]);
        }
      }
      restart_sfile.scanField("scale0_",scale_);
      restart_sfile.scanField("offset0_",offset_);

      for(unsigned i=0; i<nsel_; i++) {
        std::string msg;
        Tools::convert(i,msg);
        double tmp_w;
        restart_sfile.scanField("weight_"+msg,tmp_w);
        if(master) {
          average_weights_[i][replica_] = tmp_w;
          if(nrep_>1) multi_sim_comm.Sum(&average_weights_[i][0], nrep_);
        }
        comm.Sum(&average_weights_[i][0], nrep_);
      }

    }
    restart_sfile.scanField();
    restart_sfile.close();
  }

  /* If DSIGMA is not yet initialised do it now */
  for(unsigned i=0; i<sigma_max_.size(); i++) if(Dsigma_[i]==-1) Dsigma_[i] = 0.05*(sigma_max_[i] - sigma_min_[i]);

  addComponentWithDerivatives("score");
  componentIsNotPeriodic("score");
  valueScore=getPntrToComponent("score");

  if(do_reweight_) {
    addComponent("biasDer");
    componentIsNotPeriodic("biasDer");
    addComponent("weight");
    componentIsNotPeriodic("weight");
  }

  addComponent("neff");
  componentIsNotPeriodic("neff");

  if(doscale_ || doregres_zero_) {
    addComponent("scale");
    componentIsNotPeriodic("scale");
    valueScale=getPntrToComponent("scale");
  }

  if(dooffset_) {
    addComponent("offset");
    componentIsNotPeriodic("offset");
    valueOffset=getPntrToComponent("offset");
  }

  if(dooffset_||doscale_) {
    addComponent("acceptScale");
    componentIsNotPeriodic("acceptScale");
    valueAcceptScale=getPntrToComponent("acceptScale");
  }

  if(noise_type_==GENERIC) {
    addComponent("acceptFT");
    componentIsNotPeriodic("acceptFT");
    valueAcceptFT=getPntrToComponent("acceptFT");
  }

  addComponent("acceptSigma");
  componentIsNotPeriodic("acceptSigma");
  valueAccept=getPntrToComponent("acceptSigma");

  if(noise_type_==MGAUSS||noise_type_==MOUTLIERS||noise_type_==GENERIC) {
    for(unsigned i=0; i<sigma_mean2_.size(); ++i) {
      std::string num; Tools::convert(i,num);
      addComponent("sigmaMean-"+num); componentIsNotPeriodic("sigmaMean-"+num);
      valueSigmaMean.push_back(getPntrToComponent("sigmaMean-"+num));
      getPntrToComponent("sigmaMean-"+num)->set(std::sqrt(sigma_mean2_[i]));
      addComponent("sigma-"+num); componentIsNotPeriodic("sigma-"+num);
      valueSigma.push_back(getPntrToComponent("sigma-"+num));
      getPntrToComponent("sigma-"+num)->set(sigma_[i]);
      if(noise_type_==GENERIC) {
        addComponent("ftilde-"+num); componentIsNotPeriodic("ftilde-"+num);
        valueFtilde.push_back(getPntrToComponent("ftilde-"+num));
      }
    }
  } else {
    addComponent("sigmaMean"); componentIsNotPeriodic("sigmaMean");
    valueSigmaMean.push_back(getPntrToComponent("sigmaMean"));
    getPntrToComponent("sigmaMean")->set(std::sqrt(sigma_mean2_[0]));
    addComponent("sigma"); componentIsNotPeriodic("sigma");
    valueSigma.push_back(getPntrToComponent("sigma"));
    getPntrToComponent("sigma")->set(sigma_[0]);
  }

  switch(noise_type_) {
  case GENERIC:
    log.printf("  with general metainference ");
    if(gen_likelihood_==LIKE_GAUSS) log.printf(" and a gaussian likelihood\n");
    else if(gen_likelihood_==LIKE_LOGN) log.printf(" and a log-normal likelihood\n");
    log.printf("  ensemble average parameter sampled with a step %lf of sigma_mean\n", Dftilde_);
    break;
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
    // check that the scale value is the same for all replicas
    std::vector<double> dummy_scale(nrep_,0);
    if(master&&nrep_>1) {
      dummy_scale[replica_] = scale_;
      multi_sim_comm.Sum(dummy_scale);
    }
    comm.Sum(dummy_scale);
    for(unsigned i=1; i<nrep_; i++) {
      std::string msg = "The SCALE value must be the same for all replicas: check your input or restart file";
      if(dummy_scale[i]!=dummy_scale[0]) plumed_merror(msg);
    }
    log.printf("  sampling a common scaling factor with:\n");
    log.printf("    initial scale parameter %f\n",scale_);
    if(scale_prior_==SC_GAUSS) {
      log.printf("    gaussian prior with mean %f and width %f\n",scale_mu_,Dscale_);
    }
    if(scale_prior_==SC_FLAT) {
      log.printf("    flat prior between %f - %f\n",scale_min_,scale_max_);
      log.printf("    maximum MC move of scale parameter %f\n",Dscale_);
    }
  }

  if(dooffset_) {
    // check that the offset value is the same for all replicas
    std::vector<double> dummy_offset(nrep_,0);
    if(master&&nrep_>1) {
      dummy_offset[replica_] = offset_;
      multi_sim_comm.Sum(dummy_offset);
    }
    comm.Sum(dummy_offset);
    for(unsigned i=1; i<nrep_; i++) {
      std::string msg = "The OFFSET value must be the same for all replicas: check your input or restart file";
      if(dummy_offset[i]!=dummy_offset[0]) plumed_merror(msg);
    }
    log.printf("  sampling a common offset with:\n");
    log.printf("    initial offset parameter %f\n",offset_);
    if(offset_prior_==SC_GAUSS) {
      log.printf("    gaussian prior with mean %f and width %f\n",offset_mu_,Doffset_);
    }
    if(offset_prior_==SC_FLAT) {
      log.printf("    flat prior between %f - %f\n",offset_min_,offset_max_);
      log.printf("    maximum MC move of offset parameter %f\n",Doffset_);
    }
  }

  if(doregres_zero_)
    log.printf("  doing regression with zero intercept with stride: %d\n", nregres_zero_);

  log.printf("  number of experimental data points %u\n",narg);
  log.printf("  number of replicas %u\n",nrep_);
  log.printf("  initial data uncertainties");
  for(unsigned i=0; i<sigma_.size(); ++i) log.printf(" %f", sigma_[i]);
  log.printf("\n");
  log.printf("  minimum data uncertainties");
  for(unsigned i=0; i<sigma_.size(); ++i) log.printf(" %f",sigma_min_[i]);
  log.printf("\n");
  log.printf("  maximum data uncertainties");
  for(unsigned i=0; i<sigma_.size(); ++i) log.printf(" %f",sigma_max_[i]);
  log.printf("\n");
  log.printf("  maximum MC move of data uncertainties");
  for(unsigned i=0; i<sigma_.size(); ++i) log.printf(" %f",Dsigma_[i]);
  log.printf("\n");
  log.printf("  temperature of the system %f\n",kbt_);
  log.printf("  MC steps %u\n",MCsteps_);
  log.printf("  initial standard errors of the mean");
  for(unsigned i=0; i<sigma_mean2_.size(); ++i) log.printf(" %f", std::sqrt(sigma_mean2_[i]));
  log.printf("\n");

  //resize the number of metainference derivatives and the number of back-calculated data
  metader_.resize(narg, 0.);
  calc_data_.resize(narg, 0.);

  log<<"  Bibliography "<<plumed.cite("Bonomi, Camilloni, Cavalli, Vendruscolo, Sci. Adv. 2, e150117 (2016)");
  if(do_reweight_) log<<plumed.cite("Bonomi, Camilloni, Vendruscolo, Sci. Rep. 6, 31232 (2016)");
  if(do_optsigmamean_>0) log<<plumed.cite("Loehr, Jussupow, Camilloni, J. Chem. Phys. 146, 165102 (2017)");
  log<<plumed.cite("Bonomi, Camilloni, Bioinformatics, 33, 3999 (2017)");
  log<<"\n";
}

void MetainferenceBase::Selector()
{
  iselect = 0;
  // set the value of selector for  REM-like stuff
  if(selector_.length()>0) iselect = static_cast<unsigned>(plumed.passMap[selector_]);
}

double MetainferenceBase::getEnergySP(const std::vector<double> &mean, const std::vector<double> &sigma,
                                      const double scale, const double offset)
{
  const double scale2 = scale*scale;
  const double sm2    = sigma_mean2_[0];
  const double ss2    = sigma[0]*sigma[0] + scale2*sm2;
  const double sss    = sigma[0]*sigma[0] + sm2;

  double ene = 0.0;
  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(ene)
  {
    #pragma omp for reduction( + : ene)
    for(unsigned i=0; i<narg; ++i) {
      const double dev = scale*mean[i]-parameters[i]+offset;
      const double a2 = 0.5*dev*dev + ss2;
      if(sm2 > 0.0) {
        ene += std::log(2.0*a2/(1.0-std::exp(-a2/sm2)));
      }
      else {
        ene += std::log(2.0*a2);
      }
    }
  }
  // add one single Jeffrey's prior and one normalisation per data point
  ene += 0.5*std::log(sss) + static_cast<double>(narg)*0.5*std::log(0.5*M_PI*M_PI/ss2);
  if(doscale_ || doregres_zero_) ene += 0.5*std::log(sss);
  if(dooffset_) ene += 0.5*std::log(sss);
  return kbt_ * ene;
}

double MetainferenceBase::getEnergySPE(const std::vector<double> &mean, const std::vector<double> &sigma,
                                       const double scale, const double offset)
{
  const double scale2 = scale*scale;
  double ene = 0.0;
  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(ene)
  {
    #pragma omp for reduction( + : ene)
    for(unsigned i=0; i<narg; ++i) {
      const double sm2 = sigma_mean2_[i];
      const double ss2 = sigma[i]*sigma[i] + scale2*sm2;
      const double sss = sigma[i]*sigma[i] + sm2;
      const double dev = scale*mean[i]-parameters[i]+offset;
      const double a2  = 0.5*dev*dev + ss2;
      if(sm2 > 0.0) {
        ene += 0.5*std::log(sss) + 0.5*std::log(0.5*M_PI*M_PI/ss2) + std::log(2.0*a2/(1.0-std::exp(-a2/sm2)));
      }
      else {
        ene += 0.5*std::log(sss) + 0.5*std::log(0.5*M_PI*M_PI/ss2) + std::log(2.0*a2);
      }
      if(doscale_ || doregres_zero_)  ene += 0.5*std::log(sss);
      if(dooffset_) ene += 0.5*std::log(sss);
    }
  }
  return kbt_ * ene;
}

double MetainferenceBase::getEnergyMIGEN(const std::vector<double> &mean, const std::vector<double> &ftilde, const std::vector<double> &sigma,
    const double scale, const double offset)
{
  double ene = 0.0;
  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(ene)
  {
    #pragma omp for reduction( + : ene)
    for(unsigned i=0; i<narg; ++i) {
      const double inv_sb2  = 1./(sigma[i]*sigma[i]);
      const double inv_sm2  = 1./sigma_mean2_[i];
      double devb = 0;
      if(gen_likelihood_==LIKE_GAUSS)     devb = scale*ftilde[i]-parameters[i]+offset;
      else if(gen_likelihood_==LIKE_LOGN) devb = std::log(scale*ftilde[i]/parameters[i]);
      double devm = mean[i] - ftilde[i];
      // deviation + normalisation + jeffrey
      double normb = 0.;
      if(gen_likelihood_==LIKE_GAUSS)     normb = -0.5*std::log(0.5/M_PI*inv_sb2);
      else if(gen_likelihood_==LIKE_LOGN) normb = -0.5*std::log(0.5/M_PI*inv_sb2/(parameters[i]*parameters[i]));
      const double normm         = -0.5*std::log(0.5/M_PI*inv_sm2);
      const double jeffreys      = -0.5*std::log(2.*inv_sb2);
      ene += 0.5*devb*devb*inv_sb2 + 0.5*devm*devm*inv_sm2 + normb + normm + jeffreys;
      if(doscale_ || doregres_zero_)  ene += jeffreys;
      if(dooffset_) ene += jeffreys;
    }
  }
  return kbt_ * ene;
}

double MetainferenceBase::getEnergyGJ(const std::vector<double> &mean, const std::vector<double> &sigma,
                                      const double scale, const double offset)
{
  const double scale2  = scale*scale;
  const double inv_s2  = 1./(sigma[0]*sigma[0] + scale2*sigma_mean2_[0]);
  const double inv_sss = 1./(sigma[0]*sigma[0] + sigma_mean2_[0]);

  double ene = 0.0;
  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(ene)
  {
    #pragma omp for reduction( + : ene)
    for(unsigned i=0; i<narg; ++i) {
      double dev = scale*mean[i]-parameters[i]+offset;
      ene += 0.5*dev*dev*inv_s2;
    }
  }
  const double normalisation = -0.5*std::log(0.5/M_PI*inv_s2);
  const double jeffreys = -0.5*std::log(2.*inv_sss);
  // add Jeffrey's prior in case one sigma for all data points + one normalisation per datapoint
  ene += jeffreys + static_cast<double>(narg)*normalisation;
  if(doscale_ || doregres_zero_)  ene += jeffreys;
  if(dooffset_) ene += jeffreys;

  return kbt_ * ene;
}

double MetainferenceBase::getEnergyGJE(const std::vector<double> &mean, const std::vector<double> &sigma,
                                       const double scale, const double offset)
{
  const double scale2 = scale*scale;

  double ene = 0.0;
  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(ene)
  {
    #pragma omp for reduction( + : ene)
    for(unsigned i=0; i<narg; ++i) {
      const double inv_s2  = 1./(sigma[i]*sigma[i] + scale2*sigma_mean2_[i]);
      const double inv_sss = 1./(sigma[i]*sigma[i] + sigma_mean2_[i]);
      double dev = scale*mean[i]-parameters[i]+offset;
      // deviation + normalisation + jeffrey
      const double normalisation = -0.5*std::log(0.5/M_PI*inv_s2);
      const double jeffreys      = -0.5*std::log(2.*inv_sss);
      ene += 0.5*dev*dev*inv_s2 + normalisation + jeffreys;
      if(doscale_ || doregres_zero_)  ene += jeffreys;
      if(dooffset_) ene += jeffreys;
    }
  }
  return kbt_ * ene;
}

void MetainferenceBase::moveTilde(const std::vector<double> &mean_, double &old_energy)
{
  std::vector<double> new_ftilde(sigma_.size());
  new_ftilde = ftilde_;

  // change all tildes
  for(unsigned j=0; j<sigma_.size(); j++) {
    const double r3 = random[0].Gaussian();
    const double ds3 = Dftilde_*std::sqrt(sigma_mean2_[j])*r3;
    new_ftilde[j] = ftilde_[j] + ds3;
  }
  // calculate new energy
  double new_energy = getEnergyMIGEN(mean_,new_ftilde,sigma_,scale_,offset_);

  // accept or reject
  const double delta = ( new_energy - old_energy ) / kbt_;
  // if delta is negative always accept move
  if( delta <= 0.0 ) {
    old_energy = new_energy;
    ftilde_ = new_ftilde;
    MCacceptFT_++;
    // otherwise extract random number
  } else {
    const double s = random[0].RandU01();
    if( s < std::exp(-delta) ) {
      old_energy = new_energy;
      ftilde_ = new_ftilde;
      MCacceptFT_++;
    }
  }
}

void MetainferenceBase::moveScaleOffset(const std::vector<double> &mean_, double &old_energy)
{
  double new_scale = scale_;

  if(doscale_) {
    if(scale_prior_==SC_FLAT) {
      const double r1 = random[1].Gaussian();
      const double ds1 = Dscale_*r1;
      new_scale += ds1;
      // check boundaries
      if(new_scale > scale_max_) {new_scale = 2.0 * scale_max_ - new_scale;}
      if(new_scale < scale_min_) {new_scale = 2.0 * scale_min_ - new_scale;}
    } else {
      const double r1 = random[1].Gaussian();
      const double ds1 = 0.5*(scale_mu_-new_scale)+Dscale_*std::exp(1)/M_PI*r1;
      new_scale += ds1;
    }
  }

  double new_offset = offset_;

  if(dooffset_) {
    if(offset_prior_==SC_FLAT) {
      const double r1 = random[1].Gaussian();
      const double ds1 = Doffset_*r1;
      new_offset += ds1;
      // check boundaries
      if(new_offset > offset_max_) {new_offset = 2.0 * offset_max_ - new_offset;}
      if(new_offset < offset_min_) {new_offset = 2.0 * offset_min_ - new_offset;}
    } else {
      const double r1 = random[1].Gaussian();
      const double ds1 = 0.5*(offset_mu_-new_offset)+Doffset_*std::exp(1)/M_PI*r1;
      new_offset += ds1;
    }
  }

  // calculate new energy
  double new_energy = 0.;

  switch(noise_type_) {
  case GAUSS:
    new_energy = getEnergyGJ(mean_,sigma_,new_scale,new_offset);
    break;
  case MGAUSS:
    new_energy = getEnergyGJE(mean_,sigma_,new_scale,new_offset);
    break;
  case OUTLIERS:
    new_energy = getEnergySP(mean_,sigma_,new_scale,new_offset);
    break;
  case MOUTLIERS:
    new_energy = getEnergySPE(mean_,sigma_,new_scale,new_offset);
    break;
  case GENERIC:
    new_energy = getEnergyMIGEN(mean_,ftilde_,sigma_,new_scale,new_offset);
    break;
  }

  // for the scale/offset we need to consider the total energy
  std::vector<double> totenergies(2);
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
  if( delta <= 0.0 ) {
    old_energy = new_energy;
    scale_ = new_scale;
    offset_ = new_offset;
    MCacceptScale_++;
    // otherwise extract random number
  } else {
    double s = random[1].RandU01();
    if( s < std::exp(-delta) ) {
      old_energy = new_energy;
      scale_ = new_scale;
      offset_ = new_offset;
      MCacceptScale_++;
    }
  }
}

void MetainferenceBase::moveSigmas(const std::vector<double> &mean_, double &old_energy, const unsigned i, const std::vector<unsigned> &indices, bool &breaknow)
{
  std::vector<double> new_sigma(sigma_.size());
  new_sigma = sigma_;

  // change MCchunksize_ sigmas
  if (MCchunksize_ > 0) {
    if ((MCchunksize_ * i) >= sigma_.size()) {
      // This means we are not moving any sigma, so we should break immediately
      breaknow = true;
    }

    // change random sigmas
    for(unsigned j=0; j<MCchunksize_; j++) {
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
      if(new_sigma[index] > sigma_max_[index]) {new_sigma[index] = 2.0 * sigma_max_[index] - new_sigma[index];}
      if(new_sigma[index] < sigma_min_[index]) {new_sigma[index] = 2.0 * sigma_min_[index] - new_sigma[index];}
    }
  } else {
    // change all sigmas
    for(unsigned j=0; j<sigma_.size(); j++) {
      const double r2 = random[0].Gaussian();
      const double ds2 = Dsigma_[j]*r2;
      new_sigma[j] = sigma_[j] + ds2;
      // check boundaries
      if(new_sigma[j] > sigma_max_[j]) {new_sigma[j] = 2.0 * sigma_max_[j] - new_sigma[j];}
      if(new_sigma[j] < sigma_min_[j]) {new_sigma[j] = 2.0 * sigma_min_[j] - new_sigma[j];}
    }
  }

  if (breaknow) {
    // We didnt move any sigmas, so no sense in evaluating anything
    return;
  }

  // calculate new energy
  double new_energy = 0.;
  switch(noise_type_) {
  case GAUSS:
    new_energy = getEnergyGJ(mean_,new_sigma,scale_,offset_);
    break;
  case MGAUSS:
    new_energy = getEnergyGJE(mean_,new_sigma,scale_,offset_);
    break;
  case OUTLIERS:
    new_energy = getEnergySP(mean_,new_sigma,scale_,offset_);
    break;
  case MOUTLIERS:
    new_energy = getEnergySPE(mean_,new_sigma,scale_,offset_);
    break;
  case GENERIC:
    new_energy = getEnergyMIGEN(mean_,ftilde_,new_sigma,scale_,offset_);
    break;
  }

  // accept or reject
  const double delta = ( new_energy - old_energy ) / kbt_;
  // if delta is negative always accept move
  if( delta <= 0.0 ) {
    old_energy = new_energy;
    sigma_ = new_sigma;
    MCaccept_++;
    // otherwise extract random number
  } else {
    const double s = random[0].RandU01();
    if( s < std::exp(-delta) ) {
      old_energy = new_energy;
      sigma_ = new_sigma;
      MCaccept_++;
    }
  }
}

double MetainferenceBase::doMonteCarlo(const std::vector<double> &mean_)
{
  // calculate old energy with the updated coordinates
  double old_energy=0.;

  switch(noise_type_) {
  case GAUSS:
    old_energy = getEnergyGJ(mean_,sigma_,scale_,offset_);
    break;
  case MGAUSS:
    old_energy = getEnergyGJE(mean_,sigma_,scale_,offset_);
    break;
  case OUTLIERS:
    old_energy = getEnergySP(mean_,sigma_,scale_,offset_);
    break;
  case MOUTLIERS:
    old_energy = getEnergySPE(mean_,sigma_,scale_,offset_);
    break;
  case GENERIC:
    old_energy = getEnergyMIGEN(mean_,ftilde_,sigma_,scale_,offset_);
    break;
  }

  // do not run MC if this is a replica-exchange trial
  if(!getExchangeStep()) {

    // Create std::vector of random sigma indices
    std::vector<unsigned> indices;
    if (MCchunksize_ > 0) {
      for (unsigned j=0; j<sigma_.size(); j++) {
        indices.push_back(j);
      }
      random[2].Shuffle(indices);
    }
    bool breaknow = false;

    // cycle on MC steps
    for(unsigned i=0; i<MCsteps_; ++i) {
      MCtrial_++;
      // propose move for ftilde
      if(noise_type_==GENERIC) moveTilde(mean_, old_energy);
      // propose move for scale and/or offset
      if(doscale_||dooffset_) moveScaleOffset(mean_, old_energy);
      // propose move for sigma
      moveSigmas(mean_, old_energy, i, indices, breaknow);
      // exit from the loop if this is the case
      if(breaknow) break;
    }

    /* save the result of the sampling */
    /* ftilde */
    if(noise_type_==GENERIC) {
      double accept = static_cast<double>(MCacceptFT_) / static_cast<double>(MCtrial_);
      valueAcceptFT->set(accept);
      for(unsigned i=0; i<sigma_.size(); i++) valueFtilde[i]->set(ftilde_[i]);
    }
    /* scale and offset */
    if(doscale_ || doregres_zero_) valueScale->set(scale_);
    if(dooffset_) valueOffset->set(offset_);
    if(doscale_||dooffset_) {
      double accept = static_cast<double>(MCacceptScale_) / static_cast<double>(MCtrial_);
      valueAcceptScale->set(accept);
    }
    /* sigmas */
    for(unsigned i=0; i<sigma_.size(); i++) valueSigma[i]->set(sigma_[i]);
    double accept = static_cast<double>(MCaccept_) / static_cast<double>(MCtrial_);
    valueAccept->set(accept);
  }

  // here we sum the score over the replicas to get the full metainference score that we save as a bias
  if(master) {
    if(nrep_>1) multi_sim_comm.Sum(old_energy);
  } else {
    old_energy=0;
  }
  comm.Sum(old_energy);

  // this is the energy with current coordinates and parameters
  return old_energy;
}

/*
   In the following energy-force functions we don't add the normalisation and the jeffreys priors
   because they are not needed for the forces, the correct MetaInference energy is the one calculated
   in the Monte-Carlo
*/

void MetainferenceBase::getEnergyForceSP(const std::vector<double> &mean, const std::vector<double> &dmean_x,
    const std::vector<double> &dmean_b)
{
  const double scale2 = scale_*scale_;
  const double sm2    = sigma_mean2_[0];
  const double ss2    = sigma_[0]*sigma_[0] + scale2*sm2;
  std::vector<double> f(narg,0);

  if(master) {
    #pragma omp parallel num_threads(OpenMP::getNumThreads())
    {
      #pragma omp for
      for(unsigned i=0; i<narg; ++i) {
        const double dev = scale_*mean[i]-parameters[i]+offset_;
        const double a2 = 0.5*dev*dev + ss2;
        if(sm2 > 0.0) {
          const double t = std::exp(-a2/sm2);
          const double dt = 1./t;
          const double dit = 1./(1.-dt);
          f[i] = -scale_*dev*(dit/sm2 + 1./a2);
        }
        else {
          f[i] = -scale_*dev*(1./a2);
        }
      }
    }
    // collect contribution to forces and energy from other replicas
    if(nrep_>1) multi_sim_comm.Sum(&f[0],narg);
  }
  // intra-replica summation
  comm.Sum(&f[0],narg);

  double w_tmp = 0.;
  for(unsigned i=0; i<narg; ++i) {
    setMetaDer(i, -kbt_*f[i]*dmean_x[i]);
    w_tmp += kbt_*f[i]*dmean_b[i];
  }

  if(do_reweight_) {
    setArgDerivatives(valueScore, -w_tmp);
    getPntrToComponent("biasDer")->set(-w_tmp);
  }
}

void MetainferenceBase::getEnergyForceSPE(const std::vector<double> &mean, const std::vector<double> &dmean_x,
    const std::vector<double> &dmean_b)
{
  const double scale2 = scale_*scale_;
  std::vector<double> f(narg,0);

  if(master) {
    #pragma omp parallel num_threads(OpenMP::getNumThreads())
    {
      #pragma omp for
      for(unsigned i=0; i<narg; ++i) {
        const double sm2 = sigma_mean2_[i];
        const double ss2 = sigma_[i]*sigma_[i] + scale2*sm2;
        const double dev = scale_*mean[i]-parameters[i]+offset_;
        const double a2  = 0.5*dev*dev + ss2;
        if(sm2 > 0.0) {
          const double t   = std::exp(-a2/sm2);
          const double dt  = 1./t;
          const double dit = 1./(1.-dt);
          f[i] = -scale_*dev*(dit/sm2 + 1./a2);
        }
        else {
          f[i] = -scale_*dev*(1./a2);
        }
      }
    }
    // collect contribution to forces and energy from other replicas
    if(nrep_>1) multi_sim_comm.Sum(&f[0],narg);
  }
  comm.Sum(&f[0],narg);

  double w_tmp = 0.;
  for(unsigned i=0; i<narg; ++i) {
    setMetaDer(i, -kbt_ * dmean_x[i] * f[i]);
    w_tmp += kbt_ * dmean_b[i] *f[i];
  }

  if(do_reweight_) {
    setArgDerivatives(valueScore, -w_tmp);
    getPntrToComponent("biasDer")->set(-w_tmp);
  }
}

void MetainferenceBase::getEnergyForceGJ(const std::vector<double> &mean, const std::vector<double> &dmean_x,
    const std::vector<double> &dmean_b)
{
  const double scale2 = scale_*scale_;
  double inv_s2=0.;

  if(master) {
    inv_s2 = 1./(sigma_[0]*sigma_[0] + scale2*sigma_mean2_[0]);
    if(nrep_>1) multi_sim_comm.Sum(inv_s2);
  }
  comm.Sum(inv_s2);

  double w_tmp = 0.;
  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(w_tmp)
  {
    #pragma omp for reduction( + : w_tmp)
    for(unsigned i=0; i<narg; ++i) {
      const double dev = scale_*mean[i]-parameters[i]+offset_;
      const double mult = dev*scale_*inv_s2;
      setMetaDer(i, kbt_*dmean_x[i]*mult);
      w_tmp += kbt_*dmean_b[i]*mult;
    }
  }

  if(do_reweight_) {
    setArgDerivatives(valueScore, w_tmp);
    getPntrToComponent("biasDer")->set(w_tmp);
  }
}

void MetainferenceBase::getEnergyForceGJE(const std::vector<double> &mean, const std::vector<double> &dmean_x,
    const std::vector<double> &dmean_b)
{
  const double scale2 = scale_*scale_;
  std::vector<double> inv_s2(sigma_.size(),0.);

  if(master) {
    for(unsigned i=0; i<sigma_.size(); ++i) inv_s2[i] = 1./(sigma_[i]*sigma_[i] + scale2*sigma_mean2_[i]);
    if(nrep_>1) multi_sim_comm.Sum(&inv_s2[0],sigma_.size());
  }
  comm.Sum(&inv_s2[0],sigma_.size());

  double w_tmp = 0.;
  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(w_tmp)
  {
    #pragma omp for reduction( + : w_tmp)
    for(unsigned i=0; i<narg; ++i) {
      const double dev  = scale_*mean[i]-parameters[i]+offset_;
      const double mult = dev*scale_*inv_s2[i];
      setMetaDer(i, kbt_*dmean_x[i]*mult);
      w_tmp += kbt_*dmean_b[i]*mult;
    }
  }

  if(do_reweight_) {
    setArgDerivatives(valueScore, w_tmp);
    getPntrToComponent("biasDer")->set(w_tmp);
  }
}

void MetainferenceBase::getEnergyForceMIGEN(const std::vector<double> &mean, const std::vector<double> &dmean_x, const std::vector<double> &dmean_b)
{
  std::vector<double> inv_s2(sigma_.size(),0.);
  std::vector<double> dev(sigma_.size(),0.);
  std::vector<double> dev2(sigma_.size(),0.);

  for(unsigned i=0; i<sigma_.size(); ++i) {
    inv_s2[i]   = 1./sigma_mean2_[i];
    if(master) {
      dev[i]  = (mean[i]-ftilde_[i]);
      dev2[i] = dev[i]*dev[i];
    }
  }
  if(master&&nrep_>1) {
    multi_sim_comm.Sum(&dev[0],dev.size());
    multi_sim_comm.Sum(&dev2[0],dev2.size());
  }
  comm.Sum(&dev[0],dev.size());
  comm.Sum(&dev2[0],dev2.size());

  double dene_b = 0.;
  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(dene_b)
  {
    #pragma omp for reduction( + : dene_b)
    for(unsigned i=0; i<narg; ++i) {
      const double dene_x  = kbt_*inv_s2[i]*dmean_x[i]*dev[i];
      dene_b += kbt_*inv_s2[i]*dmean_b[i]*dev[i];
      setMetaDer(i, dene_x);
    }
  }

  if(do_reweight_) {
    setArgDerivatives(valueScore, dene_b);
    getPntrToComponent("biasDer")->set(dene_b);
  }
}

void MetainferenceBase::get_weights(double &weight, double &norm, double &neff)
{
  const double dnrep = static_cast<double>(nrep_);
  // calculate the weights either from BIAS
  if(do_reweight_) {
    std::vector<double> bias(nrep_,0);
    if(master) {
      bias[replica_] = getArgument(0);
      if(nrep_>1) multi_sim_comm.Sum(&bias[0], nrep_);
    }
    comm.Sum(&bias[0], nrep_);

    // accumulate weights
    if(!firstTimeW[iselect]) {
      for(unsigned i=0; i<nrep_; ++i) {
        const double delta=bias[i]-average_weights_[iselect][i];
        average_weights_[iselect][i]+=decay_w_*delta;
      }
    } else {
      firstTimeW[iselect] = false;
      for(unsigned i=0; i<nrep_; ++i) average_weights_[iselect][i] = bias[i];
    }

    // set average back into bias and set norm to one
    const double maxbias = *(std::max_element(average_weights_[iselect].begin(), average_weights_[iselect].end()));
    for(unsigned i=0; i<nrep_; ++i) bias[i] = std::exp((average_weights_[iselect][i]-maxbias)/kbt_);
    // set local weight, norm and weight variance
    weight = bias[replica_];
    double w2=0.;
    for(unsigned i=0; i<nrep_; ++i) {
      w2 += bias[i]*bias[i];
      norm += bias[i];
    }
    neff = norm*norm/w2;
    getPntrToComponent("weight")->set(weight/norm);
  } else {
    // or arithmetic ones
    neff = dnrep;
    weight = 1.0;
    norm = dnrep;
  }
  getPntrToComponent("neff")->set(neff);
}

void MetainferenceBase::get_sigma_mean(const double weight, const double norm, const double neff, const std::vector<double> &mean)
{
  const double dnrep    = static_cast<double>(nrep_);
  std::vector<double> sigma_mean2_tmp(sigma_mean2_.size());

  if(do_optsigmamean_>0) {
    // remove first entry of the history std::vector
    if(sigma_mean2_last_[iselect][0].size()==optsigmamean_stride_&&optsigmamean_stride_>0)
      for(unsigned i=0; i<narg; ++i) sigma_mean2_last_[iselect][i].erase(sigma_mean2_last_[iselect][i].begin());
    /* this is the current estimate of sigma mean for each argument
       there is one of this per argument in any case  because it is
       the maximum among these to be used in case of GAUSS/OUTLIER */
    std::vector<double> sigma_mean2_now(narg,0);
    if(master) {
      for(unsigned i=0; i<narg; ++i) sigma_mean2_now[i] = weight*(calc_data_[i]-mean[i])*(calc_data_[i]-mean[i]);
      if(nrep_>1) multi_sim_comm.Sum(&sigma_mean2_now[0], narg);
    }
    comm.Sum(&sigma_mean2_now[0], narg);
    for(unsigned i=0; i<narg; ++i) sigma_mean2_now[i] *= 1.0/(neff-1.)/norm;

    // add sigma_mean2 to history
    if(optsigmamean_stride_>0) {
      for(unsigned i=0; i<narg; ++i) sigma_mean2_last_[iselect][i].push_back(sigma_mean2_now[i]);
    } else {
      for(unsigned i=0; i<narg; ++i) if(sigma_mean2_now[i] > sigma_mean2_last_[iselect][i][0]) sigma_mean2_last_[iselect][i][0] = sigma_mean2_now[i];
    }

    if(noise_type_==MGAUSS||noise_type_==MOUTLIERS||noise_type_==GENERIC) {
      for(unsigned i=0; i<narg; ++i) {
        /* set to the maximum in history std::vector */
        sigma_mean2_tmp[i] = *max_element(sigma_mean2_last_[iselect][i].begin(), sigma_mean2_last_[iselect][i].end());
        /* the standard error of the mean */
        valueSigmaMean[i]->set(std::sqrt(sigma_mean2_tmp[i]));
        if(noise_type_==GENERIC) {
          sigma_min_[i] = std::sqrt(sigma_mean2_tmp[i]);
          if(sigma_[i] < sigma_min_[i]) sigma_[i] = sigma_min_[i];
        }
      }
    } else if(noise_type_==GAUSS||noise_type_==OUTLIERS) {
      // find maximum for each data point
      std::vector <double> max_values;
      for(unsigned i=0; i<narg; ++i) max_values.push_back(*max_element(sigma_mean2_last_[iselect][i].begin(), sigma_mean2_last_[iselect][i].end()));
      // find maximum across data points
      const double max_now = *max_element(max_values.begin(), max_values.end());
      // set new value
      sigma_mean2_tmp[0] = max_now;
      valueSigmaMean[0]->set(std::sqrt(sigma_mean2_tmp[0]));
    }
    // endif sigma mean optimization
    // start sigma max optimization
    if(do_optsigmamean_>1&&!sigmamax_opt_done_) {
      for(unsigned i=0; i<sigma_max_.size(); i++) {
        if(sigma_max_est_[i]<sigma_mean2_tmp[i]&&optimized_step_>optsigmamean_stride_) sigma_max_est_[i]=sigma_mean2_tmp[i];
        // ready to set once and for all the value of sigma_max
        if(optimized_step_==N_optimized_step_) {
          sigmamax_opt_done_=true;
          for(unsigned i=0; i<sigma_max_.size(); i++) {
            sigma_max_[i]=std::sqrt(sigma_max_est_[i]*dnrep);
            Dsigma_[i] = 0.05*(sigma_max_[i] - sigma_min_[i]);
            if(sigma_[i]>sigma_max_[i]) sigma_[i]=sigma_max_[i];
          }
        }
      }
      optimized_step_++;
    }
    // end sigma max optimization
  } else {
    if(noise_type_==MGAUSS||noise_type_==MOUTLIERS||noise_type_==GENERIC) {
      for(unsigned i=0; i<narg; ++i) {
        sigma_mean2_tmp[i] = sigma_mean2_last_[iselect][i][0];
        valueSigmaMean[i]->set(std::sqrt(sigma_mean2_tmp[i]));
      }
    } else if(noise_type_==GAUSS||noise_type_==OUTLIERS) {
      sigma_mean2_tmp[0] = sigma_mean2_last_[iselect][0][0];
      valueSigmaMean[0]->set(std::sqrt(sigma_mean2_tmp[0]));
    }
  }

  sigma_mean2_ = sigma_mean2_tmp;
}

void MetainferenceBase::replica_averaging(const double weight, const double norm, std::vector<double> &mean, std::vector<double> &dmean_b)
{
  if(master) {
    for(unsigned i=0; i<narg; ++i) mean[i] = weight/norm*calc_data_[i];
    if(nrep_>1) multi_sim_comm.Sum(&mean[0], narg);
  }
  comm.Sum(&mean[0], narg);
  // set the derivative of the mean with respect to the bias
  for(unsigned i=0; i<narg; ++i) dmean_b[i] = weight/norm/kbt_*(calc_data_[i]-mean[i])*decay_w_;

  // this is only for generic metainference
  if(firstTime) {ftilde_ = mean; firstTime = false;}
}

void MetainferenceBase::do_regression_zero(const std::vector<double> &mean)
{
// parameters[i] = scale_ * mean[i]: find scale_ with linear regression
  double num = 0.0;
  double den = 0.0;
  for(unsigned i=0; i<parameters.size(); ++i) {
    num += mean[i] * parameters[i];
    den += mean[i] * mean[i];
  }
  if(den>0) {
    scale_ = num / den;
  } else {
    scale_ = 1.0;
  }
}

double MetainferenceBase::getScore()
{
  /* Metainference */
  /* 1) collect weights */
  double weight = 0.;
  double neff = 0.;
  double norm = 0.;
  get_weights(weight, norm, neff);

  /* 2) calculate average */
  std::vector<double> mean(getNarg(),0);
  // this is the derivative of the mean with respect to the argument
  std::vector<double> dmean_x(getNarg(),weight/norm);
  // this is the derivative of the mean with respect to the bias
  std::vector<double> dmean_b(getNarg(),0);
  // calculate it
  replica_averaging(weight, norm, mean, dmean_b);

  /* 3) calculates parameters */
  get_sigma_mean(weight, norm, neff, mean);

  // in case of regression with zero intercept, calculate scale
  if(doregres_zero_ && getStep()%nregres_zero_==0) do_regression_zero(mean);

  /* 4) run monte carlo */
  double ene = doMonteCarlo(mean);

  // calculate bias and forces
  switch(noise_type_) {
  case GAUSS:
    getEnergyForceGJ(mean, dmean_x, dmean_b);
    break;
  case MGAUSS:
    getEnergyForceGJE(mean, dmean_x, dmean_b);
    break;
  case OUTLIERS:
    getEnergyForceSP(mean, dmean_x, dmean_b);
    break;
  case MOUTLIERS:
    getEnergyForceSPE(mean, dmean_x, dmean_b);
    break;
  case GENERIC:
    getEnergyForceMIGEN(mean, dmean_x, dmean_b);
    break;
  }

  return ene;
}

void MetainferenceBase::writeStatus()
{
  if(!doscore_) return;
  sfile_.rewind();
  sfile_.printField("time",getTimeStep()*getStep());
  //nsel
  for(unsigned i=0; i<sigma_mean2_last_.size(); i++) {
    std::string msg_i,msg_j;
    Tools::convert(i,msg_i);
    std::vector <double> max_values;
    //narg
    for(unsigned j=0; j<narg; ++j) {
      Tools::convert(j,msg_j);
      std::string msg = msg_i+"_"+msg_j;
      if(noise_type_==MGAUSS||noise_type_==MOUTLIERS||noise_type_==GENERIC) {
        sfile_.printField("sigmaMean_"+msg,std::sqrt(*max_element(sigma_mean2_last_[i][j].begin(), sigma_mean2_last_[i][j].end())));
      } else {
        // find maximum for each data point
        max_values.push_back(*max_element(sigma_mean2_last_[i][j].begin(), sigma_mean2_last_[i][j].end()));
      }
    }
    if(noise_type_==GAUSS||noise_type_==OUTLIERS) {
      // find maximum across data points
      const double max_now = std::sqrt(*max_element(max_values.begin(), max_values.end()));
      Tools::convert(0,msg_j);
      std::string msg = msg_i+"_"+msg_j;
      sfile_.printField("sigmaMean_"+msg, max_now);
    }
  }
  for(unsigned i=0; i<sigma_.size(); ++i) {
    std::string msg;
    Tools::convert(i,msg);
    sfile_.printField("sigma_"+msg,sigma_[i]);
  }
  for(unsigned i=0; i<sigma_max_.size(); ++i) {
    std::string msg;
    Tools::convert(i,msg);
    sfile_.printField("sigma_max_"+msg,sigma_max_[i]);
  }
  if(noise_type_==GENERIC) {
    for(unsigned i=0; i<ftilde_.size(); ++i) {
      std::string msg;
      Tools::convert(i,msg);
      sfile_.printField("ftilde_"+msg,ftilde_[i]);
    }
  }
  sfile_.printField("scale0_",scale_);
  sfile_.printField("offset0_",offset_);
  for(unsigned i=0; i<average_weights_.size(); i++) {
    std::string msg_i;
    Tools::convert(i,msg_i);
    sfile_.printField("weight_"+msg_i,average_weights_[i][replica_]);
  }
  sfile_.printField();
  sfile_.flush();
}

}
}

