/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
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
#include "bias/ActionRegister.h"
#include "core/Atoms.h"
#include "core/PlumedMain.h"
#include <fstream>

using namespace std;

namespace PLMD {
namespace isdb {

//+PLUMEDOC BIAS CALIBEREXTERNAL
/*
Add a time-dependent, harmonic restraint on one or more variables.


*/
//+ENDPLUMEDOC


class CaliberExternal : public bias::Bias {
public:
  explicit CaliberExternal(const ActionOptions&);
  void calculate();
  static void registerKeywords( Keywords& keys );
private:
  vector<double> time;
  vector<double> var;
  vector<double> min;
  double kappa;
  double deltat;
  double mult;
  bool     firstTime_;
  bool     adaptive_;
  bool     master;  
  unsigned replica_;
  unsigned nrep_;
  vector<double> sigma_mean2_;
  // optimize sigma mean
  vector < vector <double> > sigma_mean2_last_;
  unsigned optsigmamean_stride_;

  void get_sigma_mean(const double fact, const vector<double> &mean);
  void replica_averaging(const double fact, vector<double> &mean);
};

PLUMED_REGISTER_ACTION(CaliberExternal,"CALIBER")

void CaliberExternal::registerKeywords( Keywords& keys ) {
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","FILE","the name of the file containing the CV values in function of time");
  keys.add("compulsory","KAPPA","the array of force constants");
  keys.addFlag("ADAPTIVE",false,"abmd");
  keys.addOutputComponent("x0","default","the instantaneous value of the center of the potential");
}

CaliberExternal::CaliberExternal(const ActionOptions&ao):
  PLUMED_BIAS_INIT(ao),
  firstTime_(true),
  adaptive_(false)
{
  string filename;
  double tempT, tempVar;
  parse("KAPPA",mult);
  parse("FILE",filename);
  if( filename.length()==0 ) error("No external variable file was specified");
  parseFlag("ADAPTIVE",adaptive_);

  checkRead();

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
  const unsigned narg = getNumberOfArguments();
  min.resize(narg);
  sigma_mean2_.resize(narg,0.00001);
  sigma_mean2_last_.resize(narg);
  for(unsigned j=0; j<narg; j++) sigma_mean2_last_[j].push_back(sigma_mean2_[j]);
  optsigmamean_stride_ = 10;
  log.printf("  External variable from file %s\n",filename.c_str());

// read varfile
  std::ifstream varfile(filename.c_str());

  while (!varfile.eof()) {
    varfile >> tempT >> tempVar;
    time.push_back(tempT);
    var.push_back(tempVar);
  }

  varfile.close();

  deltat = time[1] - time[0];

  addComponent("x0"); componentIsNotPeriodic("x0");
  addComponent("kappa"); componentIsNotPeriodic("kappa");
  addComponent("mean"); componentIsNotPeriodic("mean");
  addComponent("min"); componentIsNotPeriodic("min");
}



void CaliberExternal::get_sigma_mean(const double fact, const vector<double> &mean)
{
  const unsigned narg = getNumberOfArguments();
  const double dnrep    = static_cast<double>(nrep_);
  vector<double> sigma_mean2_tmp(sigma_mean2_.size(), 0.);

  // remove first entry of the history vector
  if(sigma_mean2_last_[0].size()==optsigmamean_stride_&&optsigmamean_stride_>0)
    for(unsigned i=0; i<narg; ++i) sigma_mean2_last_[i].erase(sigma_mean2_last_[i].begin());
  /* this is the current estimate of sigma mean for each argument
     there is one of this per argument in any case  because it is
     the maximum among these to be used in case of GAUSS/OUTLIER */
  vector<double> sigma_mean2_now(narg,0);
  if(master) {
    for(unsigned i=0; i<narg; ++i) {
      double tmp  = getArgument(i)-mean[i];
      sigma_mean2_now[i] = fact*tmp*tmp;
    }
    if(nrep_>1) multi_sim_comm.Sum(&sigma_mean2_now[0], narg);
  }
  comm.Sum(&sigma_mean2_now[0], narg);
  for(unsigned i=0; i<narg; ++i) sigma_mean2_now[i] /= dnrep;

  // add sigma_mean2 to history
  for(unsigned i=0; i<narg; ++i) sigma_mean2_last_[i].push_back(sigma_mean2_now[i]);

  for(unsigned i=0; i<narg; ++i) {
    /* set to the maximum in history vector */
    sigma_mean2_tmp[i] = *max_element(sigma_mean2_last_[i].begin(), sigma_mean2_last_[i].end());
  }
  // endif sigma optimization

  sigma_mean2_ = sigma_mean2_tmp;
}

void CaliberExternal::replica_averaging(const double fact, vector<double> &mean)
{
  const unsigned narg = getNumberOfArguments();
  if(master) {
    for(unsigned i=0; i<narg; ++i) mean[i] = fact*getArgument(i);
    if(nrep_>1) multi_sim_comm.Sum(&mean[0], narg);
  }
  comm.Sum(&mean[0], narg);
}


void CaliberExternal::calculate() {
  long int now=getStep();
  double x0, dnow = static_cast<double>(now);
  const unsigned narg = getNumberOfArguments();
  int tindex;


  tindex = now/static_cast<int>(deltat);

  x0 = var[tindex]   * (1. - (dnow - time[tindex])/(time[tindex+1] - time[tindex]) ) + \
       var[tindex+1] * (1. - (time[tindex+1] - dnow)/(time[tindex+1] - time[tindex]) );


  const double dnrep    = static_cast<double>(nrep_);
  double fact = 1.0/dnrep;

  // calculate the mean
  vector<double> mean(narg,0);
  // this is the derivative of the mean with respect to the argument
  vector<double> dmean_x(narg,fact);
  // calculate it
  replica_averaging(fact, mean);
  // get kappa
  get_sigma_mean(fact, mean);

  if(firstTime_) {
    if(adaptive_) for(unsigned i=0; i<narg; ++i) min[i] = mean[i];
    firstTime_=false;
  }

  double ene=0;
  if(!adaptive_){
   for(unsigned i=0; i<narg; ++i) {
    kappa = mult*dnrep/sigma_mean2_[i];
    const double cv=difference(i,x0,mean[i]); // this gives: getArgument(i) - x0
    const double f=-kappa*cv*dmean_x[i];
    setOutputForce(i,f);
    ene+=0.5*kappa*cv*cv;
    getPntrToComponent("kappa")->set(kappa);
    getPntrToComponent("mean")->set(mean[i]);
   };
  } else {
   for(unsigned i=0; i<narg; ++i) {
    if(mean[i]>min[i]) min[i]=mean[i];
    else {
     kappa = mult*dnrep/sigma_mean2_[i];
     const double cv=difference(i,min[i],mean[i]); // this gives: getArgument(i) - x0
     const double f=-kappa*cv*dmean_x[i];
     setOutputForce(i,f);
     ene+=0.5*kappa*cv*cv;
     getPntrToComponent("kappa")->set(kappa);
     getPntrToComponent("mean")->set(mean[i]);
     getPntrToComponent("min")->set(min[i]);
    }
   }
  }

  // we will need to add back the calculation of the work
  setBias(ene);
  getPntrToComponent("x0")->set(x0);
}

}
}


