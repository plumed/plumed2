/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2018-2020 The plumed team
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

//+PLUMEDOC ISDB_BIAS CALIBER
/*
Add a time-dependent, harmonic restraint on one or more variables.

This allows implementing a maximum caliber restraint on one or more experimental time series by replica-averaged restrained simulations.
See \cite Capelli:2018jt .

The time resolved experiments are read from a text file and intermediate values are obtained by splines.

\par Examples

In the following example a restraint is applied on the time evolution of a saxs spectrum

\plumedfile
MOLINFO STRUCTURE=first.pdb

# Define saxs variable
SAXS ...
LABEL=saxs
ATOMISTIC
ATOMS=1-436
QVALUE1=0.02 # Q-value at which calculate the scattering
QVALUE2=0.0808
QVALUE3=0.1264
QVALUE4=0.1568
QVALUE5=0.172
QVALUE6=0.1872
QVALUE7=0.2176
QVALUE8=0.2328
QVALUE9=0.248
QVALUE10=0.2632
QVALUE11=0.2936
QVALUE12=0.3088
QVALUE13=0.324
QVALUE14=0.3544
QVALUE15=0.4
... SAXS


#define the caliber restraint
CALIBER ...
  ARG=(saxs\.q_.*)
  FILE=expsaxs.dat
  KAPPA=10
  LABEL=cal0
  STRIDE=10
  REGRES_ZERO=200
  AVERAGING=200
... CALIBER
\endplumedfile

In particular the file expsaxs.dat contains the time traces for the 15 intensities at the selected scattering lengths, organized as time, q_1, etc.
The strength of the bias is automatically evaluated from the standard error of the mean over AVERAGING steps and multiplied by KAPPA. This is useful when working with multiple experimental data
Because \ref SAXS is usually defined in a manner that is irrespective of a scaling factor the scaling is evaluated from a linear fit every REGRES_ZERO step. Alternatively it can be given as a fixed constant as SCALE.
The bias is here applied every tenth step.

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
  vector< vector<double> > dvar;
  double   mult;
  double   scale_;
  bool     master;
  unsigned replica_;
  unsigned nrep_;
  // scale and offset regression
  bool doregres_zero_;
  int  nregres_zero_;
  // force constant
  unsigned optsigmamean_stride_;
  vector<double> sigma_mean2_;
  vector< vector<double> > sigma_mean2_last_;
  vector<Value*> x0comp;
  vector<Value*> kcomp;
  vector<Value*> mcomp;
  Value* valueScale;

  void get_sigma_mean(const double fact, const vector<double> &mean);
  void replica_averaging(const double fact, vector<double> &mean);
  double getSpline(const unsigned iarg);
  void do_regression_zero(const vector<double> &mean);
};

PLUMED_REGISTER_ACTION(Caliber,"CALIBER")

void Caliber::registerKeywords( Keywords& keys ) {
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.addFlag("NOENSEMBLE",false,"don't perform any replica-averaging");
  keys.add("compulsory","FILE","the name of the file containing the time-resolved values");
  keys.add("compulsory","KAPPA","a force constant, this can be use to scale a constant estimated on-the-fly using AVERAGING");
  keys.add("optional","AVERAGING", "Stride for calculation of the optimum kappa, if 0 only KAPPA is used.");
  keys.add("compulsory","TSCALE","1.0","Apply a time scaling on the experimental time scale");
  keys.add("compulsory","SCALE","1.0","Apply a constant scaling on the data provided as arguments");
  keys.add("optional","REGRES_ZERO","stride for regression with zero offset");
  keys.addOutputComponent("x0","default","the instantaneous value of the center of the potential");
  keys.addOutputComponent("mean","default","the current average value of the calculated observable");
  keys.addOutputComponent("kappa","default","the current force constant");
  keys.addOutputComponent("scale","REGRES_ZERO","the current scaling constant");
}

Caliber::Caliber(const ActionOptions&ao):
  PLUMED_BIAS_INIT(ao),
  mult(0),
  scale_(1),
  doregres_zero_(false),
  nregres_zero_(0),
  optsigmamean_stride_(0)
{
  parse("KAPPA",mult);
  string filename;
  parse("FILE",filename);
  if( filename.length()==0 ) error("No external variable file was specified");
  unsigned averaging=0;
  parse("AVERAGING", averaging);
  if(averaging>0) optsigmamean_stride_ = averaging;
  double tscale=1.0;
  parse("TSCALE", tscale);
  if(tscale<=0.) error("The time scale factor must be greater than 0.");
  parse("SCALE", scale_);
  if(scale_==0.) error("The time scale factor cannot be 0.");
  // regression with zero intercept
  parse("REGRES_ZERO", nregres_zero_);
  if(nregres_zero_>0) {
    // set flag
    doregres_zero_=true;
    log.printf("  doing regression with zero intercept with stride: %d\n", nregres_zero_);
  }


  bool noensemble = false;
  parseFlag("NOENSEMBLE", noensemble);

  checkRead();

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

  const unsigned narg = getNumberOfArguments();
  sigma_mean2_.resize(narg,1);
  sigma_mean2_last_.resize(narg);
  for(unsigned j=0; j<narg; j++) sigma_mean2_last_[j].push_back(0.000001);

  log.printf("  Time resolved data from file %s\n",filename.c_str());
  std::ifstream varfile(filename.c_str());
  if(varfile.fail()) error("Cannot open "+filename);
  var.resize(narg);
  dvar.resize(narg);
  while (!varfile.eof()) {
    double tempT, tempVar;
    varfile >> tempT;
    time.push_back(tempT/tscale);
    for(unsigned i=0; i<narg; i++) {
      varfile >> tempVar;
      var[i].push_back(tempVar);
    }
  }
  varfile.close();

  const double deltat = time[1] - time[0];
  for(unsigned i=0; i<narg; i++) {
    for(unsigned j=0; j<var[i].size(); j++) {
      if(j==0) dvar[i].push_back((var[i][j+1] - var[i][j])/(deltat));
      else if(j==var[i].size()-1) dvar[i].push_back((var[i][j] - var[i][j-1])/(deltat));
      else dvar[i].push_back((var[i][j+1] - var[i][j-1])/(2.*deltat));
    }
  }

  for(unsigned i=0; i<narg; i++) {
    std::string num; Tools::convert(i,num);
    addComponent("x0-"+num); componentIsNotPeriodic("x0-"+num); x0comp.push_back(getPntrToComponent("x0-"+num));
    addComponent("kappa-"+num); componentIsNotPeriodic("kappa-"+num); kcomp.push_back(getPntrToComponent("kappa-"+num));
    addComponent("mean-"+num); componentIsNotPeriodic("mean-"+num); mcomp.push_back(getPntrToComponent("mean-"+num));
  }

  if(doregres_zero_) {
    addComponent("scale");
    componentIsNotPeriodic("scale");
    valueScale=getPntrToComponent("scale");
  }

  log<<"  Bibliography "<<plumed.cite("Capelli, Tiana, Camilloni, J Chem Phys, 148, 184114");
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

double Caliber::getSpline(const unsigned iarg)
{
  const double deltat = time[1] - time[0];
  const int tindex = static_cast<int>(getTime()/deltat);

  unsigned start, end;
  start=tindex;
  if(tindex+1<var[iarg].size()) end=tindex+2;
  else end=var[iarg].size();

  double value=0;
  for(unsigned ipoint=start; ipoint<end; ++ipoint) {
    double grid=var[iarg][ipoint];
    double dder=dvar[iarg][ipoint];
    double yy=0.;
    if(fabs(grid)>0.0000001) yy=-dder/grid;

    int x0=1;
    if(ipoint==tindex) x0=0;

    double X=fabs((getTime()-time[tindex])/deltat-(double)x0);
    double X2=X*X;
    double X3=X2*X;
    double C=(1.0-3.0*X2+2.0*X3) - (x0?-1.0:1.0)*yy*(X-2.0*X2+X3)*deltat;

    value+=grid*C;
  }
  return value;
}

void Caliber::do_regression_zero(const vector<double> &mean)
{
// parameters[i] = scale_ * mean[i]: find scale_ with linear regression
  double num = 0.0;
  double den = 0.0;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    num += mean[i] * getSpline(i);
    den += mean[i] * mean[i];
  }
  if(den>0) {
    scale_ = num / den;
  } else {
    scale_ = 1.0;
  }
}

void Caliber::calculate()
{
  const unsigned narg = getNumberOfArguments();
  const double dnrep = static_cast<double>(nrep_);
  const double fact = 1.0/dnrep;

  vector<double> mean(narg,0);
  vector<double> dmean_x(narg,fact);
  replica_averaging(fact, mean);
  if(optsigmamean_stride_>0) get_sigma_mean(fact, mean);

  // in case of regression with zero intercept, calculate scale
  if(doregres_zero_ && getStep()%nregres_zero_==0) do_regression_zero(mean);

  double ene=0;
  for(unsigned i=0; i<narg; ++i) {
    const double x0 = getSpline(i);
    const double kappa = mult*dnrep/sigma_mean2_[i];
    const double cv=difference(i,x0,scale_*mean[i]);
    const double f=-kappa*cv*dmean_x[i]/scale_;
    setOutputForce(i,f);
    ene+=0.5*kappa*cv*cv;
    x0comp[i]->set(x0);
    kcomp[i]->set(kappa);
    mcomp[i]->set(mean[i]);
  }

  if(doregres_zero_) valueScale->set(scale_);

  setBias(ene);
}

}
}


