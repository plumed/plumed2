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

//+PLUMEDOC BIAS CALIBER
/*
Add a time-dependent, harmonic restraint on one or more variables.

\par Examples

In the following example

*/
//+ENDPLUMEDOC


class Caliber : public bias::Bias {
public:
  explicit Caliber(const ActionOptions&);
  void calculate();
  static void registerKeywords( Keywords& keys );
private:
  vector<double> time;
  vector< vector<double> > var;
  double   deltat;
  double   mult;
  bool     master;
  unsigned replica_;
  unsigned nrep_;
  unsigned optsigmamean_stride_;
  vector<double> sigma_mean2_;
  vector< vector<double> > sigma_mean2_last_;

  void get_sigma_mean(const double fact, const vector<double> &mean);
  void replica_averaging(const double fact, vector<double> &mean);
};

PLUMED_REGISTER_ACTION(Caliber,"CALIBER")

void Caliber::registerKeywords( Keywords& keys ) {
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","FILE","the name of the file containing the time-resolved values");
  keys.add("compulsory","KAPPA","a force constant, this can be use to scale a constant estimanted on-the-fly using AVERAGING");
  keys.add("optional","AVERAGING", "Stride for calculation of the optimum kappa, if 0 only KAPPA is used.");
  keys.addOutputComponent("x0","default","the instantaneous value of the center of the potential");
  keys.addOutputComponent("mean","default","the current average value of the calculated observable");
  keys.addOutputComponent("kappa","default","the current force constant");
}

Caliber::Caliber(const ActionOptions&ao):
  PLUMED_BIAS_INIT(ao),
  mult(0),
  optsigmamean_stride_(0)
{
  parse("KAPPA",mult);
  string filename;
  parse("FILE",filename);
  if( filename.length()==0 ) error("No external variable file was specified");
  unsigned averaging=0;
  parse("AVERAGING", averaging);
  if(averaging>0) optsigmamean_stride_ = averaging;

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
  sigma_mean2_.resize(narg,1);
  sigma_mean2_last_.resize(narg);
  for(unsigned j=0; j<narg; j++) sigma_mean2_last_[j].push_back(0.000001);

  log.printf("  Time resolved data from file %s\n",filename.c_str());
  std::ifstream varfile(filename.c_str());
  var.resize(narg);
  while (!varfile.eof()) {
    double tempT, tempVar;
    varfile >> tempT;
    time.push_back(tempT);
    for(unsigned i=0;i<narg;i++) {
      varfile >> tempVar;
      var[i].push_back(tempVar);
    }
  }
  varfile.close();

  deltat = time[1] - time[0];

  for(unsigned i=0;i<narg;i++) {
    std::string num; Tools::convert(i,num);
    addComponent("sigmaMean_"+num); componentIsNotPeriodic("sigmaMean_"+num);
    addComponent("x0_"+num); componentIsNotPeriodic("x0_"+num);
    addComponent("kappa_"+num); componentIsNotPeriodic("kappa_"+num);
    addComponent("mean_"+num); componentIsNotPeriodic("mean_"+num);
  }
}

void Caliber::get_sigma_mean(const double fact, const vector<double> &mean)
{
  const unsigned narg = getNumberOfArguments();
  const double dnrep = static_cast<double>(nrep_);

  if(sigma_mean2_last_[0].size()==optsigmamean_stride_) for(unsigned i=0; i<narg; ++i) sigma_mean2_last_[i].erase(sigma_mean2_last_[i].begin());
  vector<double> sigma_mean2_now(narg,0);
  if(master) {
    for(unsigned i=0; i<narg; ++i) {
      double tmp = getArgument(i)-mean[i];
      sigma_mean2_now[i] = fact*tmp*tmp;
    }
    if(nrep_>1) multi_sim_comm.Sum(&sigma_mean2_now[0], narg);
  }
  comm.Sum(&sigma_mean2_now[0], narg);

  for(unsigned i=0; i<narg; ++i) {
    sigma_mean2_last_[i].push_back(sigma_mean2_now[i]/dnrep);
    sigma_mean2_[i] = *max_element(sigma_mean2_last_[i].begin(), sigma_mean2_last_[i].end());
  }
}

void Caliber::replica_averaging(const double fact, vector<double> &mean)
{
  const unsigned narg = getNumberOfArguments();
  if(master) {
    for(unsigned i=0; i<narg; ++i) mean[i] = fact*getArgument(i);
    if(nrep_>1) multi_sim_comm.Sum(&mean[0], narg);
  }
  comm.Sum(&mean[0], narg);
}


void Caliber::calculate()
{
  const long int now = getStep(); 
  const double dnow = static_cast<double>(now);
  const unsigned narg = getNumberOfArguments();
  const int tindex = now/static_cast<int>(deltat);
  const double dnrep = static_cast<double>(nrep_);
  double fact = 1.0/dnrep;

  vector<double> mean(narg,0);
  vector<double> dmean_x(narg,fact);
  replica_averaging(fact, mean);
  if(optsigmamean_stride_>0) get_sigma_mean(fact, mean);

  double ene=0;
  for(unsigned i=0; i<narg; ++i) {
    double x0 = var[i][tindex]   * (1. - (dnow - time[tindex])/(time[tindex+1] - time[tindex]) ) + \
                var[i][tindex+1] * (1. - (time[tindex+1] - dnow)/(time[tindex+1] - time[tindex]) );
    double kappa = mult*dnrep/sigma_mean2_[i];
    const double cv=difference(i,x0,mean[i]);
    const double f=-kappa*cv*dmean_x[i];
    setOutputForce(i,f);
    ene+=0.5*kappa*cv*cv;
    std::string num; Tools::convert(i,num);
    getPntrToComponent("kappa_"+num)->set(kappa);
    getPntrToComponent("x0_"+num)->set(x0);
    getPntrToComponent("mean_"+num)->set(mean[i]);
  }

  setBias(ene);
}

}
}


