/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2020 The plumed team
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
#include "Bias.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "core/ActionRegister.h"
#include "core/ActionWithValue.h"
#include "tools/Communicator.h"
#include "tools/File.h"

// The original implementation of this method was contributed
// by Andrea Cesari (andreacesari90@gmail.com).
// Copyright has been then transferred to PLUMED developers
// (see https://github.com/plumed/plumed2/blob/master/.github/CONTRIBUTING.md)

namespace PLMD {
namespace bias {

//+PLUMEDOC BIAS MAXENT
/*
Add a linear biasing potential on one or more variables \f$f_{i}\left(\boldsymbol{x}\right)\f$ satisfying the maximum entropy principle as proposed in Ref. \cite cesari2016maxent .

\warning
    Notice that syntax is still under revision and might change

The resulting biasing potential is given by:
\f[
  V_{BIAS}(\boldsymbol{x},t)=K_{B}T\sum_{i=1}^{\#arguments}f_{i}(\boldsymbol{x},t)\lambda_{i}(t)
\f]
Lagrangian multipliers \f$ \lambda_{i}\f$ are updated, every PACE steps, according to the following update rule:
\f[
\lambda_{i}=\lambda_{i}+\frac{k_{i}}{1+\frac{t}{\tau_{i}}}\left(f_{exp,i}+\xi_{i}\lambda_{i}-f_{i}(\boldsymbol{x})\right)
\f]
\f$k\f$ set the initial value of the learning rate and its units are \f$[observable]^{-2}ps^{-1}\f$. This can be set with the keyword KAPPA.
The number of components for any KAPPA vector must be equal to the number of arguments of the action.

Variable \f$ \xi_{i}(\lambda) \f$ is related to the chosen prior to model experimental errors. If a GAUSSIAN prior is used then:
\f[
\xi_{i}(\lambda)=-\lambda_{i}\sigma^{2}
\f]
where \f$ \sigma \f$ is the typical expected error on the observable \f$ f_i\f$.
For a LAPLACE prior:
\f[
\xi_{i}(\lambda)=-\frac{\lambda_{i}\sigma^{2}}{1-\frac{\lambda^{2}\sigma^{2}}{2}}

\f]
The value of \f$ \xi(\lambda,t)\f$ is written in output as a component named: argument name followed by the string _error.
Setting \f$ \sigma =0\f$ is equivalent to enforce a pure Maximum Entropy restraint without any noise modelling.
This method can be also used to enforce inequality restraint as shown in following examples.

Notice that a similar method is available as \ref EDS, although with different features and using a different optimization algorithm.

\par Examples

The following input tells plumed to restrain the distance between atoms 7 and 15
and the distance between atoms 2 and 19, at different equilibrium
values, and to print the energy of the restraint.
Lagrangian multiplier will be printed on a file called restraint.LAGMULT with a stride set by the variable PACE to 200ps.
Moreover plumed will compute the average of each Lagrangian multiplier in the window [TSTART,TEND] and use that to continue the simulations with fixed Lagrangian multipliers.
\plumedfile
DISTANCE ATOMS=7,15 LABEL=d1
DISTANCE ATOMS=2,19 LABEL=d2
MAXENT ...
ARG=d1,d2
TYPE=EQUAL
AT=0.2,0.5
KAPPA=35000.0,35000.0
TAU=0.02,0.02
PACE=200
TSTART=100
TEND=500
LABEL=restraint
... MAXENT
PRINT ARG=restraint.bias
\endplumedfile
Lagrangian multipliers will be printed on a file called restraint.bias
The following input tells plumed to restrain the distance between atoms 7 and 15
to be greater than 0.2 and to print the energy of the restraint
\plumedfile
DISTANCE ATOMS=7,15 LABEL=d
MAXENT ARG=d TYPE=INEQUAL> AT=0.02 KAPPA=35000.0 TAU=3 LABEL=restraint
PRINT ARG=restraint.bias
\endplumedfile

(See also \ref DISTANCE and \ref PRINT).

*/
//+ENDPLUMEDOC

class MaxEnt : public Bias {
  std::vector<double> at;
  std::vector<double> kappa;
  std::vector<double> lambda;
  std::vector<double> avgx;
  std::vector<double> work;
  std::vector<double> oldlambda;
  std::vector<double> tau;
  std::vector<double> avglambda;
  std::vector<double> avglambda_restart;
  std::vector<double> expected_eps;
  std::vector<double> apply_weights;
  double sigma;
  double tstart;
  double tend;
  double avgstep; //current number of samples over which to compute the average. Check if could be replaced bu getStep()
  long int pace_;
  long int stride_;
  double totalWork;
  double BetaReweightBias;
  double simtemp;
  std::vector<ActionWithValue*> biases;
  std::string type;
  std::string error_type;
  double alpha;
  double avg_counter;
  int learn_replica;
  Value* valueForce2;
  Value* valueWork;
  OFile lagmultOfile_;
  IFile ifile;
  std::string lagmultfname;
  std::string ifilesnames;
  std::string fmt;
  bool isFirstStep;
  bool reweight;
  bool no_broadcast;
  bool printFirstStep;
  std::vector<bool> done_average;
  int myrep,nrep;
public:
  explicit MaxEnt(const ActionOptions&);
  void calculate() override;
  void update() override;
  void update_lambda();
  static void registerKeywords(Keywords& keys);
  void ReadLagrangians(IFile &ifile);
  void WriteLagrangians(std::vector<double> &lagmult,OFile &file);
  double compute_error(const std::string &err_type,double l);
  double convert_lambda(const std::string &type,double lold);
  void check_lambda_boundaries(const std::string &err_type,double &l);
};
PLUMED_REGISTER_ACTION(MaxEnt,"MAXENT")

void MaxEnt::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  componentsAreNotOptional(keys);
  keys.use("ARG");
  keys.add("compulsory","KAPPA","0.0","specifies the initial value for the learning rate");
  keys.add("compulsory","TAU","Specify the dumping time for the learning rate.");
  keys.add("compulsory","TYPE","specify the restraint type. "
           "EQUAL to restrain the variable at a given equilibrium value "
           "INEQUAL< to restrain the variable to be smaller than a given value "
           "INEQUAL> to restrain the variable to be greater than a given value");
  keys.add("optional","ERROR_TYPE","specify the prior on the error to use."
           "GAUSSIAN: use a Gaussian prior "
           "LAPLACE: use a Laplace prior");
  keys.add("optional","TSTART","time from where to start averaging the Lagrangian multiplier. By default no average is computed, hence lambda is updated every PACE steps");
  keys.add("optional","TEND","time in ps where to stop to compute the average of Lagrangian multiplier. From this time until the end of the simulation Lagrangian multipliers are kept fix to the average computed between TSTART and TEND;");
  keys.add("optional","ALPHA","default=1.0; To be used with LAPLACE KEYWORD, allows to choose a prior function proportional to a Gaussian times an exponential function. ALPHA=1 correspond to the LAPLACE prior.");
  keys.add("compulsory","AT","the position of the restraint");
  keys.add("optional","SIGMA","The typical errors expected on observable");
  keys.add("optional","FILE","Lagrangian multipliers output file. The default name is: label name followed by the string .LAGMULT ");
  keys.add("optional","LEARN_REPLICA","In a multiple replica environment specify which is the reference replica. By default replica 0 will be used.");
  keys.add("optional","APPLY_WEIGHTS","Vector of weights containing 1 in correspondence of each replica that will receive the Lagrangian multiplier from the current one.");
  keys.add("optional","PACE","the frequency for Lagrangian multipliers update");
  keys.add("optional","PRINT_STRIDE","stride of Lagrangian multipliers output file. If no STRIDE is passed they are written every time they are updated (PACE).");
  keys.add("optional","FMT","specify format for Lagrangian multipliers files (useful to decrease the number of digits in regtests)");
  keys.addFlag("REWEIGHT",false,"to be used with plumed driver in order to reweight a trajectory a posteriori");
  keys.addFlag("NO_BROADCAST",false,"If active will avoid Lagrangian multipliers to be communicated to other replicas.");
  keys.add("optional","TEMP","the system temperature.  This is required if you are reweighting.");
  keys.addOutputComponent("force2","default","the instantaneous value of the squared force due to this bias potential");
  keys.addOutputComponent("work","default","the instantaneous value of the work done by the biasing force");
  keys.addOutputComponent("_work","default","the instantaneous value of the work done by the biasing force for each argument. "
                          "These quantities will named with the arguments of the bias followed by "
                          "the character string _work.");
  keys.addOutputComponent("_error","default","Instantaneous values of the discrepancy between the observable and the restraint center");
  keys.addOutputComponent("_coupling","default","Instantaneous values of Lagrangian multipliers. They are also written by default in a separate output file.");
  keys.use("RESTART");
}
MaxEnt::MaxEnt(const ActionOptions&ao):
  PLUMED_BIAS_INIT(ao),
  at(getNumberOfArguments()),
  kappa(getNumberOfArguments(),0.0),
  lambda(getNumberOfArguments(),0.0),
  avgx(getNumberOfArguments(),0.0),
  oldlambda(getNumberOfArguments(),0.0),
  tau(getNumberOfArguments(),getTimeStep()),
  avglambda(getNumberOfArguments(),0.0),
  avglambda_restart(getNumberOfArguments(),0.0),
  expected_eps(getNumberOfArguments(),0.0),
  sigma(0.0),
  pace_(100),
  stride_(100),
  alpha(1.0),
  avg_counter(0.0),
  isFirstStep(true),
  reweight(false),
  no_broadcast(false),
  printFirstStep(true),
  done_average(getNumberOfArguments(),false)
{
  if(comm.Get_rank()==0) nrep=multi_sim_comm.Get_size();
  if(comm.Get_rank()==0) myrep=multi_sim_comm.Get_rank();
  comm.Bcast(nrep,0);
  comm.Bcast(myrep,0);
  parseFlag("NO_BROADCAST",no_broadcast);
  //if(no_broadcast){
  //for(int irep=0;irep<nrep;irep++){
  //  if(irep!=myrep)
  //    apply_weights[irep]=0.0;}
  //}
  avgstep=1.0;
  tstart=-1.0;
  tend=-1.0;
  totalWork=0.0;
  learn_replica=0;

  parse("LEARN_REPLICA",learn_replica);
  parseVector("APPLY_WEIGHTS",apply_weights);
  if(apply_weights.size()==0) apply_weights.resize(nrep,1.0);
  parseVector("KAPPA",kappa);
  parseVector("AT",at);
  parseVector("TAU",tau);
  parse("TYPE",type);
  error_type="GAUSSIAN";
  parse("ERROR_TYPE",error_type);
  parse("ALPHA",alpha);
  parse("SIGMA",sigma);
  parse("TSTART",tstart);
  if(tstart <0 && tstart != -1.0) error("TSTART should be a positive number");
  parse("TEND",tend);
  if(tend<0 && tend != -1.0) error("TSTART should be a positive number");
  if(tend<tstart) error("TEND should be >= TSTART");
  lagmultfname=getLabel()+".LAGMULT";
  parse("FILE",lagmultfname);
  parse("FMT",fmt);
  parse("PACE",pace_);
  if(pace_<=0 ) error("frequency for Lagrangian multipliers update (PACE) is nonsensical");
  stride_=pace_;  //if no STRIDE is passed, then Lagrangian multipliers willbe printed at each update
  parse("PRINT_STRIDE",stride_);
  if(stride_<=0 ) error("frequency for Lagrangian multipliers printing (STRIDE) is nonsensical");
  simtemp=0.;
  parse("TEMP",simtemp);
  if(simtemp>0) simtemp*=plumed.getAtoms().getKBoltzmann();
  else simtemp=plumed.getAtoms().getKbT();
  parseFlag("REWEIGHT",reweight);
  if(simtemp<=0 && reweight) error("Set the temperature (TEMP) if you want to do reweighting.");

  checkRead();

  log.printf("  at");
  for(unsigned i=0; i<at.size(); i++) log.printf(" %f",at[i]);
  log.printf("\n");
  log.printf("  with initial learning rate for optimization of");
  for(unsigned i=0; i<kappa.size(); i++) log.printf(" %f",kappa[i]);
  log.printf("\n");
  log.printf("Dumping times for the learning rates are (ps): ");
  for(unsigned i=0; i<tau.size(); i++) log.printf(" %f",tau[i]);
  log.printf("\n");
  log.printf("Lagrangian multipliers are updated every %ld steps (PACE)\n",pace_);
  log.printf("Lagrangian multipliers output file %s\n",lagmultfname.c_str());
  log.printf("Lagrangian multipliers are written every %ld steps (PRINT_STRIDE)\n",stride_);
  if(fmt.length()>0)
    log.printf("The format for real number in Lagrangian multipliers file is %s\n",fmt.c_str());
  if(tstart >-1.0 && tend>-1.0)
    log.printf("Lagrangian multipliers are averaged from %lf ps to %lf ps\n",tstart,tend);
  if(no_broadcast)
    log.printf("Using NO_BROADCAST options. Lagrangian multipliers will not be comunicated to other replicas.\n");
  //for(int irep=0;irep<nrep;irep++){
  //  if(apply_weights[irep]!=0)
  //    log.printf("%d",irep);
  //  }
  addComponent("force2"); componentIsNotPeriodic("force2");
  addComponent("work"); componentIsNotPeriodic("work");
  valueForce2=getPntrToComponent("force2");
  valueWork=getPntrToComponent("work");

  std::string comp;
  for(unsigned i=0; i< getNumberOfArguments() ; i++) {
    comp=getPntrToArgument(i)->getName()+"_coupling";
    addComponent(comp); componentIsNotPeriodic(comp);
    comp=getPntrToArgument(i)->getName()+"_work";
    addComponent(comp); componentIsNotPeriodic(comp);
    work.push_back(0.); // initialize the work value
    comp=getPntrToArgument(i)->getName()+"_error";
    addComponent(comp); componentIsNotPeriodic(comp);
  }
  std::string fname;
  fname=lagmultfname;
  ifile.link(*this);
  if(ifile.FileExist(fname)) {
    ifile.open(fname);
    if(getRestart()) {
      log.printf("  Restarting from: %s\n",fname.c_str());
      ReadLagrangians(ifile);
      printFirstStep=false;
    }
    ifile.reset(false);
  }

  lagmultOfile_.link(*this);
  lagmultOfile_.open(fname);
  if(fmt.length()>0) {fmt=" "+fmt; lagmultOfile_.fmtField(fmt);}
}
////MEMBER FUNCTIONS
void MaxEnt::ReadLagrangians(IFile &ifile)
{
  double dummy;
  while(ifile.scanField("time",dummy)) {
    for(unsigned j=0; j<getNumberOfArguments(); ++j) {
      ifile.scanField(getPntrToArgument(j)->getName()+"_coupling",lambda[j]);
      if(dummy>=tstart && dummy <=tend)
        avglambda[j]+=lambda[j];
      if(dummy>=tend) {
        avglambda[j]=lambda[j];
        done_average[j]=true;
      }
    }
    if(dummy>=tstart && dummy <=tend)
      avg_counter++;
    ifile.scanField();
  }
}
void MaxEnt::WriteLagrangians(std::vector<double> &lagmult,OFile &file) {
  if(printFirstStep) {
    unsigned ncv=getNumberOfArguments();
    file.printField("time",getTimeStep()*getStep());
    for(unsigned i=0; i<ncv; ++i)
      file.printField(getPntrToArgument(i)->getName()+"_coupling",lagmult[i]);
    file.printField();
  } else {
    if(!isFirstStep) {
      unsigned ncv=getNumberOfArguments();
      file.printField("time",getTimeStep()*getStep());
      for(unsigned i=0; i<ncv; ++i)
        file.printField(getPntrToArgument(i)->getName()+"_coupling",lagmult[i]);
      file.printField();
    }
  }
}
double MaxEnt::compute_error(const std::string &err_type,double l) {
  double sigma2=std::pow(sigma,2.0);
  double l2=convert_lambda(type,l);
  double return_error=0;
  if(err_type=="GAUSSIAN" && sigma!=0.0)
    return_error=-l2*sigma2;
  else {
    if(err_type=="LAPLACE" && sigma!=0) {
      return_error=-l2*sigma2/(1.0-l2*l2*sigma2/(alpha+1));
    }
  }
  return return_error;
}
double MaxEnt::convert_lambda(const std::string &type,double lold) {
  double return_lambda=0;
  if(type=="EQUAL")
    return_lambda=lold;
  else {
    if(type=="INEQUAL>") {
      if(lold>0.0)
        return_lambda=0.0;
      else
        return_lambda=lold;
    }
    else {
      if(type=="INEQUAL<") {
        if(lold<0.0)
          return_lambda=0.0;
        else
          return_lambda=lold;
      }
    }
  }
  return return_lambda;
}
void MaxEnt::check_lambda_boundaries(const std::string &err_type,double &l) {
  if(err_type=="LAPLACE" && sigma !=0 ) {
    double l2=convert_lambda(err_type,l);
    if(l2 <-(std::sqrt(alpha+1)/sigma-0.01)) {
      l=-(std::sqrt(alpha+1)/sigma-0.01);
      log.printf("Lambda exceeded the allowed range\n");
    }
    if(l2>(std::sqrt(alpha+1)/sigma-0.01)) {
      l=std::sqrt(alpha+1)/sigma-0.01;
      log.printf("Lambda exceeded the allowed range\n");
    }
  }
}

void MaxEnt::update_lambda() {

  double totalWork_=0.0;
  const double time=getTime();
  const double step=getStep();
  double KbT=simtemp;
  double learning_rate;
  if(reweight)
    BetaReweightBias=plumed.getBias()/KbT;
  else
    BetaReweightBias=0.0;

  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    const double k=kappa[i];
    double cv=(getArgument(i)+compute_error(error_type,lambda[i])-at[i]);
    if(reweight)
      learning_rate=1.0*k/(1+step/tau[i]);
    else
      learning_rate=1.0*k/(1+time/tau[i]);
    lambda[i]+=learning_rate*cv*std::exp(-BetaReweightBias); //update Lagrangian multipliers and reweight them if REWEIGHT is set
    check_lambda_boundaries(error_type,lambda[i]);      //check that Lagrangians multipliers not exceed the allowed range
    if(time>=tstart && time <=tend && !done_average[i]) {
      avglambda[i]+=convert_lambda(type,lambda[i]); //compute the average of Lagrangian multipliers over the required time window
    }
    if(time>=tend && tend >=0) { //setting tend<0 will disable this feature
      if(!done_average[i]) {
        avglambda[i]=avglambda[i]/avg_counter;
        done_average[i]=true;
        lambda[i]=avglambda[i];
      }
      else
        lambda[i]=avglambda[i]; //keep Lagrangian multipliers fixed to the previously computed average.
    }
    work[i]+=(convert_lambda(type,lambda[i])-oldlambda[i])*getArgument(i); //compute the work performed in updating lambda
    totalWork_+=work[i];
    totalWork=totalWork_;
    oldlambda[i]=convert_lambda(type,lambda[i]);
  };
  if(time>=tstart && time <=tend)
    avg_counter++;
}

void MaxEnt::calculate() {
  double totf2=0.0;
  double ene=0.0;
  double KbT=simtemp;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    getPntrToComponent(getPntrToArgument(i)->getName()+"_error")->set(expected_eps[i]);
    getPntrToComponent(getPntrToArgument(i)->getName()+"_work")->set(work[i]);
    valueWork->set(totalWork);
    getPntrToComponent(getPntrToArgument(i)->getName()+"_coupling")->set(lambda[i]);
    const double f=-KbT*convert_lambda(type,lambda[i])*apply_weights[myrep];
    totf2+=f*f;
    ene+=KbT*convert_lambda(type,lambda[i])*getArgument(i)*apply_weights[myrep];
    setOutputForce(i,f);
  }
  setBias(ene);
  valueForce2->set(totf2);
}

void MaxEnt::update() {

  if(getStep()%stride_ == 0)
    WriteLagrangians(lambda,lagmultOfile_);
  if(getStep()%pace_ == 0) {
    update_lambda();
    if(!no_broadcast) {
      if(comm.Get_rank()==0) //Comunicate Lagrangian multipliers from reference replica to higher ones
        multi_sim_comm.Bcast(lambda,learn_replica);
    }
    comm.Bcast(lambda,0);
  }
  isFirstStep=false;
}

}

}
