/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2020 The plumed team
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
#include <ctime>

namespace PLMD {
namespace isdb {

//+PLUMEDOC ISDB_BIAS RESCALE
/*
Scales the value of an another action, being a Collective Variable or a Bias.

The rescaling factor is determined by a parameter defined on a logarithmic grid of dimension NBIN in the range
from 1 to MAX_RESCALE. The current value of the rescaling parameter is stored and shared across
other actions using a \ref SELECTOR. A Monte Carlo procedure is used to update the value
of the rescaling factor every MC_STRIDE steps of molecular dynamics. Well-tempered metadynamics, defined by the
parameters W0 and BIASFACTOR, is used to enhance the sampling in the space of the rescaling factor.
The well-tempered metadynamics bias potential is written to the file BFILE every BSTRIDE steps and read
when restarting the simulation using the directive \ref RESTART.

\note
Additional arguments not to be scaled, one for each bin in the rescaling parameter ladder, can be
provided at the end of the ARG list. The number of such arguments is specified by the option NOT_RESCALED.
These arguments will be not be scaled, but they will be
considered as bias potentials and used in the computation of the Metropolis
acceptance probability when proposing a move in the rescaling parameter. See example below.

\note
If PLUMED is running in a multiple-replica framework (for example using the -multi option in GROMACS),
the arguments will be summed across replicas, unless the NOT_SHARED option is used. Also, the value of the
\ref SELECTOR will be shared and thus will be the same in all replicas.

\par Examples

In this example we use \ref RESCALE to implement a simulated-tempering like approach.
The total potential energy of the system is scaled by a parameter defined on a logarithmic grid
of 5 bins in the range from 1 to 1.5.
A well-tempered metadynamics bias potential is used to ensure diffusion in the space of the rescaling
parameter.

\plumedfile
ene: ENERGY

SELECTOR NAME=GAMMA VALUE=0

RESCALE ...
LABEL=res ARG=ene TEMP=300
SELECTOR=GAMMA MAX_RESCALE=1.5 NBIN=5
W0=1000 BIASFACTOR=100.0 BSTRIDE=2000 BFILE=bias.dat
...

PRINT FILE=COLVAR ARG=* STRIDE=100
\endplumedfile

In this second example, we add to the simulated-tempering approach introduced above
one Parallel Bias metadynamics simulation (see \ref PBMETAD) for each value of the rescaling parameter.
At each moment of the simulation, only one of the \ref PBMETAD
actions is activated, based on the current value of the associated \ref SELECTOR.
The \ref PBMETAD bias potentials are not scaled, but just used in the calculation of
the Metropolis acceptance probability when proposing a move in the rescaling parameter.

\plumedfile
ene: ENERGY
d: DISTANCE ATOMS=1,2

SELECTOR NAME=GAMMA VALUE=0

pbmetad0: PBMETAD ARG=d SELECTOR=GAMMA SELECTOR_ID=0 SIGMA=0.1 PACE=500 HEIGHT=1 BIASFACTOR=8 FILE=HILLS.0
pbmetad1: PBMETAD ARG=d SELECTOR=GAMMA SELECTOR_ID=1 SIGMA=0.1 PACE=500 HEIGHT=1 BIASFACTOR=8 FILE=HILLS.1
pbmetad2: PBMETAD ARG=d SELECTOR=GAMMA SELECTOR_ID=2 SIGMA=0.1 PACE=500 HEIGHT=1 BIASFACTOR=8 FILE=HILLS.2
pbmetad3: PBMETAD ARG=d SELECTOR=GAMMA SELECTOR_ID=3 SIGMA=0.1 PACE=500 HEIGHT=1 BIASFACTOR=8 FILE=HILLS.3
pbmetad4: PBMETAD ARG=d SELECTOR=GAMMA SELECTOR_ID=4 SIGMA=0.1 PACE=500 HEIGHT=1 BIASFACTOR=8 FILE=HILLS.4

RESCALE ...
LABEL=res TEMP=300
ARG=ene,pbmetad0.bias,pbmetad1.bias,pbmetad2.bias,pbmetad3.bias,pbmetad4.bias
SELECTOR=GAMMA MAX_RESCALE=1.5 NOT_RESCALED=5 NBIN=5
W0=1000 BIASFACTOR=100.0 BSTRIDE=2000 BFILE=bias.dat
...

PRINT FILE=COLVAR ARG=* STRIDE=100
\endplumedfile



*/
//+ENDPLUMEDOC

class Rescale : public bias::Bias
{
  // gamma parameter
  std::vector<double> gamma_;
  double         w0_;
  double         biasf_;
  std::vector<double> bias_;
  std::vector<double> expo_;
  std::vector<unsigned> shared_;
  unsigned nores_;
  // bias
  unsigned int   Biasstride_;
  unsigned int   Biaspace_;
  std::string         Biasfilename_;
  bool           first_bias_;
  OFile          Biasfile_;
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
  // selector
  std::string selector_;

  // Monte Carlo
  void doMonteCarlo(unsigned igamma, double oldE, std::vector<double> args, std::vector<double> bargs);
  unsigned proposeMove(unsigned x, unsigned xmin, unsigned xmax);
  bool doAccept(double oldE, double newE);
  // read and print bias
  void read_bias();
  void print_bias(long int step);

public:
  explicit Rescale(const ActionOptions&);
  ~Rescale();
  void calculate();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(Rescale,"RESCALE")

void Rescale::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","TEMP","temperature");
  keys.add("compulsory","SELECTOR", "name of the SELECTOR used for rescaling");
  keys.add("compulsory","MAX_RESCALE","maximum values for rescaling");
  keys.add("compulsory","NBIN","number of bins for gamma grid");
  keys.add("compulsory","W0", "initial bias height");
  keys.add("compulsory","BIASFACTOR", "bias factor");
  keys.add("compulsory","BSTRIDE", "stride for writing bias");
  keys.add("compulsory","BFILE", "file name for bias");
  keys.add("optional","NOT_SHARED",   "list of arguments (from 1 to N) not summed across replicas");
  keys.add("optional","NOT_RESCALED", "these last N arguments will not be scaled");
  keys.add("optional","MC_STEPS","number of MC steps");
  keys.add("optional","MC_STRIDE","MC stride");
  keys.add("optional","PACE", "Pace for adding bias, in MC stride unit");
  componentsAreNotOptional(keys);
  keys.addOutputComponent("igamma",  "default","gamma parameter");
  keys.addOutputComponent("accgamma","default","MC acceptance for gamma");
  keys.addOutputComponent("wtbias",  "default","well-tempered bias");
}

Rescale::Rescale(const ActionOptions&ao):
  PLUMED_BIAS_INIT(ao),
  nores_(0), Biaspace_(1), first_bias_(true),
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
  comm.Sum(&nrep_,1);
  comm.Sum(&replica_,1);

  // wt-parameters
  parse("W0", w0_);
  parse("BIASFACTOR", biasf_);

  // selector name
  parse("SELECTOR", selector_);

  // number of bins for gamma ladder
  unsigned nbin;
  parse("NBIN", nbin);

  // number of bias
  parse("NOT_RESCALED", nores_);
  if(nores_>0 && nores_!=nbin) error("The number of non scaled arguments must be equal to either 0 or the number of bins");

  // maximum value of rescale
  std::vector<double> max_rescale;
  parseVector("MAX_RESCALE", max_rescale);
  // check dimension of max_rescale
  if(max_rescale.size()!=(getNumberOfArguments()-nores_))
    error("Size of MAX_RESCALE array must be equal to the number of arguments that will to be scaled");

  // calculate exponents
  double igamma_max = static_cast<double>(nbin);
  for(unsigned i=0; i<max_rescale.size(); ++i)
    expo_.push_back(std::log(max_rescale[i])/std::log(igamma_max));

  // allocate gamma grid and set bias to zero
  for(unsigned i=0; i<nbin; ++i) {
    // bias grid
    bias_.push_back(0.0);
    // gamma ladder
    double gamma = std::exp( static_cast<double>(i) / static_cast<double>(nbin-1) * std::log(igamma_max) );
    gamma_.push_back(gamma);
  }
  // print bias to file
  parse("BSTRIDE", Biasstride_);
  parse("BFILE",   Biasfilename_);

  // create vectors of shared arguments
  // by default they are all shared
  for(unsigned i=0; i<getNumberOfArguments(); ++i) shared_.push_back(1);
  // share across replicas or not
  std::vector<unsigned> not_shared;
  parseVector("NOT_SHARED", not_shared);
  // and change the non-shared
  for(unsigned i=0; i<not_shared.size(); ++i) {
    if((not_shared[i]-1)>=(getNumberOfArguments()-nores_) && nrep_>1)
      error("NOT_RESCALED args must always be shared when using multiple replicas");
    if((not_shared[i]-1)>=getNumberOfArguments())
      error("NOT_SHARED args should be lower than total number of arguments");
    shared_[not_shared[i]-1] = 0;
  }

  // monte carlo stuff
  parse("MC_STEPS",MCsteps_);
  parse("MC_STRIDE",MCstride_);
  // adjust for multiple-time steps
  MCstride_ *= getStride();
  // read bias deposition pace
  parse("PACE", Biaspace_);
  // multiply by MCstride
  Biaspace_ *= MCstride_;

  // get temperature
  double temp=0.0;
  parse("TEMP",temp);
  if(temp>0.0) kbt_=plumed.getAtoms().getKBoltzmann()*temp;
  else kbt_=plumed.getAtoms().getKbT();

  checkRead();

  log.printf("  temperature of the system in energy unit %f\n",kbt_);
  log.printf("  name of the SELECTOR use for this action %s\n",selector_.c_str());
  log.printf("  number of bins in grid %u\n",nbin);
  log.printf("  number of arguments that will not be scaled %u\n",nores_);
  if(nrep_>1) log<<"  number of arguments that will not be summed across replicas "<<not_shared.size()<<"\n";
  log.printf("  biasfactor %f\n",biasf_);
  log.printf("  initial hills height %f\n",w0_);
  log.printf("  stride to write bias to file %u\n",Biasstride_);
  log.printf("  write bias to file : %s\n",Biasfilename_.c_str());
  log.printf("  number of replicas %u\n",nrep_);
  log.printf("  number of MC steps %d\n",MCsteps_);
  log.printf("  do MC every %d steps\n", MCstride_);
  log.printf("\n");

  log << " Bibliography" << plumed.cite("Bonomi, Camilloni, Bioinformatics, 33, 3999 (2017)") << "\n";


  // add components
  addComponent("igamma");   componentIsNotPeriodic("igamma");
  addComponent("accgamma"); componentIsNotPeriodic("accgamma");
  addComponent("wtbias");   componentIsNotPeriodic("wtbias");

  // initialize random seed
  srand (time(NULL));

  // read bias if restarting
  if(getRestart()) read_bias();
}

Rescale::~Rescale()
{
  Biasfile_.close();
}

void Rescale::read_bias()
{
// open file
  auto ifile=Tools::make_unique<IFile>();
  ifile->link(*this);
  if(ifile->FileExist(Biasfilename_)) {
    ifile->open(Biasfilename_);
    // read all the lines, store last value of bias
    double MDtime;
    while(ifile->scanField("MD_time",MDtime)) {
      for(unsigned i=0; i<bias_.size(); ++i) {
        // convert i to string
        std::stringstream ss;
        ss << i;
        // label
        std::string label = "b" + ss.str();
        // read entry
        ifile->scanField(label, bias_[i]);
      }
      // new line
      ifile->scanField();
    }
    ifile->close();
  } else {
    error("Cannot find bias file "+Biasfilename_+"\n");
  }
}

unsigned Rescale::proposeMove(unsigned x, unsigned xmin, unsigned xmax)
{
  int xmin_i = static_cast<int>(xmin);
  int xmax_i = static_cast<int>(xmax);
  int dx;
  int r = rand() % 2;
  if( r % 2 == 0 ) dx = +1;
  else             dx = -1;
// new index, integer
  int x_new = static_cast<int>(x) + dx;
// check boundaries
  if(x_new >= xmax_i) x_new = xmax_i-1;
  if(x_new <  xmin_i) x_new = xmin_i;
  return static_cast<unsigned>(x_new);
}

bool Rescale::doAccept(double oldE, double newE)
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
    if( s < std::exp(-delta) ) { accept = true; }
  }
  return accept;
}

void Rescale::doMonteCarlo(unsigned igamma, double oldE,
                           std::vector<double> args, std::vector<double> bargs)
{
  double oldB, newB;

// cycle on MC steps
  for(unsigned i=0; i<MCsteps_; ++i) {
    // propose move in igamma
    unsigned new_igamma = proposeMove(igamma, 0, gamma_.size());
    // calculate new energy
    double newE = 0.0;
    for(unsigned j=0; j<args.size(); ++j) {
      // calculate energy term
      double fact = 1.0/pow(gamma_[new_igamma], expo_[j]) - 1.0;
      newE += args[j] * fact;
    }
    // calculate contributions from non-rescaled terms
    if(bargs.size()>0) {
      oldB = bias_[igamma]+bargs[igamma];
      newB = bias_[new_igamma]+bargs[new_igamma];
    } else {
      oldB = bias_[igamma];
      newB = bias_[new_igamma];
    }
    // accept or reject
    bool accept = doAccept(oldE+oldB, newE+newB);
    if(accept) {
      igamma = new_igamma;
      oldE = newE;
      MCaccgamma_++;
    }
  }
// send values of gamma to all replicas
  if(comm.Get_rank()==0) {
    if(multi_sim_comm.Get_rank()!=0) igamma = 0;
    multi_sim_comm.Sum(&igamma, 1);
  } else {
    igamma = 0;
  }
// local communication
  comm.Sum(&igamma, 1);

// set the value of gamma into passMap
  plumed.passMap[selector_]=static_cast<double>(igamma);
}

void Rescale::print_bias(long int step)
{
// if first time open the file
  if(first_bias_) {
    first_bias_ = false;
    Biasfile_.link(*this);
    Biasfile_.open(Biasfilename_);
    Biasfile_.setHeavyFlush();
    Biasfile_.fmtField("%30.5f");
  }

// write fields
  double MDtime = static_cast<double>(step)*getTimeStep();
  Biasfile_.printField("MD_time", MDtime);
  for(unsigned i=0; i<bias_.size(); ++i) {
    // convert i to string
    std::stringstream ss;
    ss << i;
    // label
    std::string label = "b" + ss.str();
    // print entry
    Biasfile_.printField(label, bias_[i]);
  }
  Biasfile_.printField();
}

void Rescale::calculate()
{
  // get the current value of the selector
  unsigned igamma = static_cast<unsigned>(plumed.passMap[selector_]);

  // collect data from other replicas
  std::vector<double> all_args(getNumberOfArguments(), 0.0);
  // first calculate arguments
  for(unsigned i=0; i<all_args.size(); ++i) {
    double arg = getArgument(i);
    // sum shared arguments across replicas
    if(shared_[i]==1) {
      if(comm.Get_rank()==0) multi_sim_comm.Sum(arg);
      else                   arg = 0.0;
      if(comm.Get_size()>1)  comm.Sum(arg);
    }
    // put into all_args
    all_args[i] = arg;
  }

  // now separate terms that should be rescaled
  std::vector<double> args;
  if(getNumberOfArguments()-nores_>0) args.resize(getNumberOfArguments()-nores_);
  for(unsigned i=0; i<args.size(); ++i)  args[i]  = all_args[i];
  // and terms that should not
  std::vector<double> bargs;
  if(nores_>0) bargs.resize(nores_);
  for(unsigned i=0; i<bargs.size(); ++i) bargs[i] = all_args[i+args.size()];

  // calculate energy and forces, only on rescaled terms
  double ene = 0.0;
  for(unsigned i=0; i<args.size(); ++i) {
    // calculate energy term
    double fact = 1.0/pow(gamma_[igamma], expo_[i]) - 1.0;
    ene += args[i] * fact;
    // add force
    setOutputForce(i, -fact);
  }

  // set to zero on the others
  for(unsigned i=0; i<bargs.size(); ++i) setOutputForce(i+args.size(), 0.0);

  // set value of the bias
  setBias(ene);
  // set value of the wt-bias
  getPntrToComponent("wtbias")->set(bias_[igamma]);
  // set values of gamma
  getPntrToComponent("igamma")->set(igamma);
  // get time step
  long int step = getStep();
  if(MCfirst_==-1) MCfirst_=step;
  // calculate gamma acceptance
  double MCtrials = std::floor(static_cast<double>(step-MCfirst_) / static_cast<double>(MCstride_))+1.0;
  double accgamma = static_cast<double>(MCaccgamma_) / static_cast<double>(MCsteps_) / MCtrials;
  getPntrToComponent("accgamma")->set(accgamma);

  // do MC at the right time step
  if(step%MCstride_==0&&!getExchangeStep()) doMonteCarlo(igamma, ene, args, bargs);

  // add well-tempered like bias
  if(step%Biaspace_==0) {
    // get updated igamma
    unsigned igamma = static_cast<unsigned>(plumed.passMap[selector_]);
    // add "Gaussian"
    double kbDT = kbt_ * ( biasf_ - 1.0 );
    bias_[igamma] += w0_ * std::exp(-bias_[igamma] / kbDT);
  }

  // print bias
  if(step%Biasstride_==0) print_bias(step);

}


}
}

