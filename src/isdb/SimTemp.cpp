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

*/
#include "bias/Bias.h"
#include "bias/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "core/Value.h"
#include "tools/File.h"
#include "tools/Random.h"
#include <cmath>
#include <ctime>

using namespace std;

namespace PLMD {
namespace isdb {

//+PLUMEDOC ISDB_BIAS SIMTEMP
/*

*/
//+ENDPLUMEDOC

class SimTemp : public bias::Bias
{
  // gamma parameter
  vector<double> gamma_;
  vector<double> w_;
  vector<double> ave_ene_;
  vector<double> count_;
  vector<unsigned> ensemble_;
  unsigned nores_;
  // status
  unsigned int init_;
  unsigned int statusstride_;
  string       statusfilename_;
  OFile        statusfile_;
  bool         first_status_;
  // temperature in kbt
  double kbt_;
  // Monte Carlo stuff
  unsigned MCsteps_;
  unsigned MCstride_;
  long int MCfirst_;
  long unsigned MCaccgamma_;
  // replica stuff
  unsigned nrep_;
  unsigned replica_;
  unsigned st_rep_;
  unsigned st_up_;
  // selector
  string selector_;

  // Monte Carlo
  void doMonteCarlo(unsigned ig, double oldE, vector<double> args, vector<double> bargs);
  unsigned proposeMove(unsigned x, unsigned xmin, unsigned xmax);
  bool doAccept(double oldE, double newE);
  // read and print status
  void read_status();
  void print_status(long int step);

public:
  explicit SimTemp(const ActionOptions&);
  ~SimTemp();
  void calculate();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(SimTemp,"SIMTEMP")

void SimTemp::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  componentsAreNotOptional(keys);
  useCustomisableComponents(keys);
  keys.use("ARG");
  keys.add("compulsory","SELECTOR", "name of the SELECTOR used for rescaling");
  keys.add("compulsory","MAX_RESCALE","maximum value of rescaling");
  keys.add("compulsory","NBIN","number of bins for gamma grid");
  keys.add("compulsory","STATUS_STRIDE", "stride for writing status");
  keys.add("compulsory","STATUS_FILE", "file name for status");
  keys.add("optional","TEMP", "temperature");
  keys.add("optional","ENSEMBLE", "list of ensemble bias (from 1 to N ) to always rescale");
  keys.add("optional","NOT_RESCALED", "these last N arguments will not be rescaled");
  keys.add("optional","MC_STEPS","number of MC steps");
  keys.add("optional","MC_STRIDE","MC stride");
  keys.add("optional","INIT", "initialize averages");
  keys.addOutputComponent("ig",  "default","gamma parameter");
  keys.addOutputComponent("accgamma","default","MC acceptance for gamma");
  keys.addOutputComponent("strep","default","replica doing ST");
  keys.addOutputComponent("w","default","ST weights");
}

SimTemp::SimTemp(const ActionOptions&ao):
  PLUMED_BIAS_INIT(ao),
  nores_(0), init_(0), first_status_(true),
  MCsteps_(1), MCstride_(1), MCfirst_(-1), MCaccgamma_(0)
{
  // set up replica stuff
  if(comm.Get_rank()==0) {
    nrep_    = multi_sim_comm.Get_size();
    replica_ = multi_sim_comm.Get_rank();
  } else {
    nrep_    = 0;
    replica_ = 0;
  }
  // local communication
  comm.Sum(&nrep_,1);
  comm.Sum(&replica_,1);

  // selector name
  parse("SELECTOR", selector_);

  // number of bins for gamma ladder
  unsigned nbin;
  parse("NBIN", nbin);

  // initialize averages
  parse("INIT", init_);
  init_ *= nbin;

  // number of terms not to be rescaled
  parse("NOT_RESCALED", nores_);
  if(nores_>0 && nores_!=nbin) error("The number of non rescaled arguments must be equal to either 0 or the number of bins");

  // maximum value of rescale
  double max_rescale;
  parse("MAX_RESCALE", max_rescale);

  // allocate gamma grid and set weights/ave energy to zero
  for(unsigned i=0; i<nbin; ++i) {
    // weights
    w_.push_back(0.0);
    // average energy
    ave_ene_.push_back(0.0);
    // counter
    count_.push_back(0.0);
    // gamma ladder
    double gamma = exp( static_cast<double>(i) / static_cast<double>(nbin-1) * std::log(max_rescale) );
    gamma_.push_back(gamma);
    // add ST weight component
    std::string num; Tools::convert(i,num);
    addComponent("w"+num); componentIsNotPeriodic("w"+num);
  }

  // print status to file
  parse("STATUS_STRIDE", statusstride_);
  parse("STATUS_FILE",   statusfilename_);

  // list of ensemble bias to always rescale
  vector<unsigned> ens;
  parseVector("ENSEMBLE", ens);
  // initialize ensemble to zero
  for(unsigned i=0; i<getNumberOfArguments(); ++i) ensemble_.push_back(0);
  // set ensemble bias to one
  for(unsigned i=0; i<ens.size(); ++i) ensemble_[ens[i]-1] = 1;

  // monte carlo stuff
  parse("MC_STEPS",MCsteps_);
  parse("MC_STRIDE",MCstride_);
  // adjust for multiple-time steps
  MCstride_ *= getStride();

  // get temperature
  double temp=0.0;
  parse("TEMP",temp);
  if(temp>0.0) kbt_=plumed.getAtoms().getKBoltzmann()*temp;
  else kbt_=plumed.getAtoms().getKbT();

  checkRead();

  log.printf("  temperature of the system in energy unit %f\n",kbt_);
  log.printf("  name of the SELECTOR use for this action %s\n",selector_.c_str());
  log.printf("  number of bins in grid %u\n",nbin);
  log.printf("  number of arguments that will not be rescaled %u\n",nores_);
  if(nrep_>1) log.printf("  number of ensemble bias %u\n",ens.size());
  log.printf("  stride for printing to status file %u\n",statusstride_);
  log.printf("  name of status file : %s\n",statusfilename_.c_str());
  log.printf("  number of replicas %u\n",nrep_);
  log.printf("  id of the replica %u\n",replica_);
  log.printf("  number of MC steps %d\n",MCsteps_);
  log.printf("  do MC every %d steps\n", MCstride_);
  if(init_>0) log.printf("  initialize averages for %d steps\n", init_);
  log.printf("\n");

  // add components
  addComponent("ig");       componentIsNotPeriodic("ig");
  addComponent("accgamma"); componentIsNotPeriodic("accgamma");
  addComponent("strep");    componentIsNotPeriodic("strep");

  // initialize random seed
  srand (time(NULL));

  // read status file if restarting
  if(getRestart()) read_status();
}

SimTemp::~SimTemp()
{
  statusfile_.close();
}

void SimTemp::read_status()
{
  double MDtime, st_rep, st_up;
// open file
  IFile *ifile = new IFile();
  ifile->link(*this);
  if(ifile->FileExist(statusfilename_)) {
    ifile->open(statusfilename_);
    // read all the lines, store last value of weights
    while(ifile->scanField("MD_time",MDtime)) {
      for(unsigned i=0; i<w_.size(); ++i) {
        // convert i to string
        std::string num; Tools::convert(i,num);
        // read entries
        ifile->scanField("ave"+num, ave_ene_[i]);
        ifile->scanField("n"+num, count_[i]);
        ifile->scanField("w"+num, w_[i]);
      }
      ifile->scanField("st_rep", st_rep);
      st_rep_ = static_cast<unsigned>(st_rep);
      ifile->scanField("st_up", st_up);
      st_up_ = static_cast<unsigned>(st_up);
      ifile->scanField("ig", plumed.passMap[selector_]);
      // new line
      ifile->scanField();
    }
    ifile->close();
  } else {
    error("Cannot find status file "+statusfilename_+"\n");
  }
  delete ifile;
}

unsigned SimTemp::proposeMove(unsigned x, unsigned xmin, unsigned xmax)
{
  int xmin_i = static_cast<int>(xmin);
  int xmax_i = static_cast<int>(xmax);
  int dx;
  int r = rand() % 2;
  if( r == 0 ) dx = +1;
  else         dx = -1;
// new index, integer
  int x_new = static_cast<int>(x) + dx;
// check boundaries
  if(x_new >= xmax_i) x_new = xmax_i-1;
  if(x_new <  xmin_i) x_new = xmin_i;
  return static_cast<unsigned>(x_new);
}

bool SimTemp::doAccept(double oldE, double newE)
{
  bool accept = false;
  // calculate delta energy
  double delta = ( newE - oldE ) / kbt_;
  // if delta is negative always accept move
  if( delta < 0.0 ) {
    accept = true;
  } else {
    // otherwise extract random number
    double s = static_cast<double>(rand()) / RAND_MAX;
    if( s < exp(-delta) ) { accept = true; }
  }
  return accept;
}

void SimTemp::doMonteCarlo(unsigned ig, double oldE,
                           vector<double> args, vector<double> bargs)
{
// only the active replica changes igamma
  if(replica_ == st_rep_) {
    // cycle on MC steps
    for(unsigned i=0; i<MCsteps_; ++i) {
      // propose move in ig
      unsigned new_ig = proposeMove(ig, 0, gamma_.size());
      // calculate new energy
      double newE = 0.0;
      for(unsigned j=0; j<args.size(); ++j) {
        // calculate energy term
        double fact = 1.0/gamma_[new_ig] - 1.0;
        newE += args[j] * fact;
      }
      // add contributions from ST weights
      double oldB = -w_[ig];
      double newB = -w_[new_ig];
      // and, in case, from non-rescaled terms (i.e., metadynamics bias...)
      if(bargs.size()>0) {
        oldB += bargs[ig];
        newB += bargs[new_ig];
      }
      // accept or reject
      bool accept = doAccept(oldE+oldB, newE+newB);
      if(accept) {
        ig = new_ig;
        oldE = newE;
        MCaccgamma_++;
      }
    }
  }
// if not active replica, or not master node of each replica
  if(replica_ != st_rep_ || comm.Get_rank()!=0) {
    ig = 0;
    MCaccgamma_ = 0;
  }
// send values of gamma/acceptance to all replicas
  if(comm.Get_rank()==0) {
    multi_sim_comm.Sum(&ig, 1);
    multi_sim_comm.Sum(&MCaccgamma_, 1);
  }
// local communication
  comm.Sum(&ig, 1);
  comm.Sum(&MCaccgamma_, 1);

// set the value of gamma into passMap
  plumed.passMap[selector_]=static_cast<double>(ig);
}

void SimTemp::print_status(long int step)
{
// if first time open the file
  if(first_status_) {
    first_status_ = false;
    statusfile_.link(*this);
    statusfile_.open(statusfilename_);
    statusfile_.setHeavyFlush();
    statusfile_.fmtField("%30.5f");
  }

// write fields
  double MDtime = static_cast<double>(step)*getTimeStep();
  statusfile_.printField("MD_time", MDtime);
  for(unsigned i=0; i<w_.size(); ++i) {
    // convert i to string
    std::string num; Tools::convert(i,num);
    // print entry
    statusfile_.printField("ave"+num, ave_ene_[i]);
    statusfile_.printField("n"+num, count_[i]);
    statusfile_.printField("w"+num, w_[i]);
  }
  double st_rep = static_cast<double>(st_rep_);
  statusfile_.printField("st_rep", st_rep);
  double st_up = static_cast<double>(st_up_);
  statusfile_.printField("st_up", st_up);
  statusfile_.printField("ig", plumed.passMap[selector_]);
  statusfile_.printField();
}

void SimTemp::calculate()
{
  // get step
  long int step = getStep();

  // get the current value of the selector
  unsigned ig = static_cast<unsigned>(plumed.passMap[selector_]);

  // if in the initialization process
  if(step < init_) {
    // calculate segment
    double seg = static_cast<double>(step) / static_cast<double>(init_) * static_cast<double>(w_.size());
    // get index
    ig = static_cast<unsigned>(floor(seg));
    // reset selector
    plumed.passMap[selector_] = static_cast<double>(ig);
    // replica 0 is doing the initialization
    st_rep_ = 0;
    st_up_ = 1;
  }

  // after initialization reset ig to zero
  if(step == init_) {
    ig = 0;
    plumed.passMap[selector_] = static_cast<double>(ig);
  }

  // now separate terms that should be rescaled
  vector<double> args;
  if(getNumberOfArguments()-nores_>0) args.resize(getNumberOfArguments()-nores_);
  for(unsigned i=0; i<args.size(); ++i)  args[i]  = getArgument(i);
  // and terms that should not
  vector<double> bargs;
  if(nores_>0) bargs.resize(nores_);
  for(unsigned i=0; i<bargs.size(); ++i) bargs[i] = getArgument(i+args.size());

  // calculate energy and forces, only on rescaled terms and active replica
  double ene = 0.0;
  for(unsigned i=0; i<args.size(); ++i) {
    // active replica or ensemble bias
    if(replica_ == st_rep_ || ensemble_[i]==1) {
      // calculate energy term
      double fact = 1.0/gamma_[ig] - 1.0;
      ene += args[i] * fact;
      // add force
      setOutputForce(i, -fact);
    } else {
      setOutputForce(i, 0.0);
    }
  }

  // set force to zero on the other arguments
  for(unsigned i=0; i<bargs.size(); ++i) setOutputForce(i+args.size(), 0.0);

  // if active replica...
  if(replica_ == st_rep_) {
    // calculate sum of rescaled terms
    double sum_ene = 0.0;
    for(unsigned i=0; i<args.size(); ++i) sum_ene += args[i];
    // update average
    ave_ene_[ig] = count_[ig] / ( count_[ig] + 1.0 ) *  ave_ene_[ig] + sum_ene / ( count_[ig] + 1.0 );
    // and increase counter
    count_[ig] += 1.0;
    // update ST weights
    if(ig>0) {
      w_[ig] = w_[ig-1] + (1.0/gamma_[ig]-1.0/gamma_[ig-1]) * ( ave_ene_[ig] + ave_ene_[ig-1] ) / 2.0;
    }
    if(ig<gamma_.size()-1) {
      w_[ig+1] = w_[ig] + (1.0/gamma_[ig+1]-1.0/gamma_[ig]) * ( ave_ene_[ig+1] + ave_ene_[ig] ) / 2.0;
    }
  }
  // put to zero if not active replica, or not master node of each replica
  if(replica_ != st_rep_ || comm.Get_rank()!=0) {
    for(unsigned i=0; i<w_.size(); ++i) {
      w_[i] = 0.0;
      ave_ene_[i] = 0.0;
      count_[i] = 0.0;
    }
  }
  // communicate stuff to non-active replicas
  if(comm.Get_rank()==0) {
    multi_sim_comm.Sum(&w_[0], w_.size());
    multi_sim_comm.Sum(&ave_ene_[0], ave_ene_.size());
    multi_sim_comm.Sum(&count_[0], count_.size());
  }
  // and local sharing
  comm.Sum(&w_[0], w_.size());
  comm.Sum(&ave_ene_[0], ave_ene_.size());
  comm.Sum(&count_[0], count_.size());

  // set value of the bias
  setBias(ene);
  // set values of gamma
  getPntrToComponent("ig")->set(ig);
  // set values of weights
  for(unsigned i=0; i<w_.size(); ++i) {
    // convert i to string
    std::string num; Tools::convert(i,num);
    // set
    getPntrToComponent("w"+num)->set(w_[i]);
  }
  // set value of replica doing ST
  getPntrToComponent("strep")->set(st_rep_);

  // after initialization process
  if(step >= init_) {
    // do MC at the right time step
    if(step%MCstride_==0&&!getExchangeStep()) doMonteCarlo(ig, ene, args, bargs);
    // calculate gamma acceptance
    if(MCfirst_==-1) MCfirst_=step;
    double MCtrials = std::floor(static_cast<double>(step-MCfirst_) / static_cast<double>(MCstride_))+1.0;
    double accgamma = static_cast<double>(MCaccgamma_) / static_cast<double>(MCsteps_) / MCtrials;
    getPntrToComponent("accgamma")->set(accgamma);
  }

  // print status
  if(step%statusstride_==0) print_status(step);

  // time to change replica to perform ST?
  if(step >= init_) {
    // recover updated gamma selector
    ig = static_cast<unsigned>(plumed.passMap[selector_]);
    // reaching top of the ladder, when going up
    if(ig==w_.size()-1 && st_up_==1) st_up_ = 0;
    // reaching bottom of the ladder, when going down
    if(ig==0 && st_up_==0) {
      // increment replica index
      st_rep_ += 1;
      // if getting to the end of the chain of replicas, restart
      if(st_rep_==nrep_) st_rep_ = 0;
      // reset up/down flag
      st_up_ = 1;
    }
  }
}


}
}

