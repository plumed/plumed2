/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2016 The plumed team
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
#include "Function.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/Communicator.h"

using namespace std;

namespace PLMD{
namespace function{

//+PLUMEDOC FUNCTION ENSEMBLE
/*
Calculates the replica averaging of a collective variable over multiple replicas.

Each collective variable is averaged separately and stored in a component labelled <em>label</em>.cvlabel.  

Note that in case of variables such as \ref CS2BACKBONE, \ref CH3SHIFTS, \ref NOE and \ref RDC it is possible
to perform the replica-averaging inside the variable, in fact in those cases are the single experimental
values that averaged before calculating the collective variable.

\par Examples
The following input tells plumed to calculate the distance between atoms 3 and 5
and the average it over the available replicas.
\verbatim
dist: DISTANCE ATOMS=3,5 
ens: ENSEMBLE ARG=dist
PRINT ARG=dist,ens.dist
\endverbatim
(See also \ref PRINT and \ref DISTANCE).

*/
//+ENDPLUMEDOC


class Ensemble :
  public Function
{
  unsigned ens_dim;
  unsigned my_repl;
  unsigned narg;
  std::vector<double> bias;
  std::vector<double> mean;
  std::vector<double> var;
  std::vector<double> var_tmp;
  bool     do_reweight;
  bool     master;
  bool     do_variance;
  double   kbt;
public:
  explicit Ensemble(const ActionOptions&);
  void     calculate();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(Ensemble,"ENSEMBLE")

void Ensemble::registerKeywords(Keywords& keys){
  Function::registerKeywords(keys);
  keys.use("ARG");
  keys.addFlag("REWEIGHT",false,"do reweight"); 
  keys.addFlag("VARIANCE",false,"calculate the variance in addition to the mean"); 
  keys.add("optional","TEMP","the system temperature - this is only needed if you are reweighting");
  ActionWithValue::useCustomisableComponents(keys);
}

Ensemble::Ensemble(const ActionOptions&ao):
Action(ao),
Function(ao),
do_reweight(false),
do_variance(false),
kbt(-1.0)
{
  parseFlag("REWEIGHT", do_reweight); 
  parseFlag("VARIANCE", do_variance); 
  double temp=0.0;
  parse("TEMP",temp);
  if(do_reweight) {
    if(temp>0.0) kbt=plumed.getAtoms().getKBoltzmann()*temp;
    else kbt=plumed.getAtoms().getKbT();
    if(kbt==0.0) error("Unless the MD engine passes the temperature to plumed, with REWEIGHT you must specify TEMP");
  }

  master = (comm.Get_rank()==0);
  if(master) {
    if(multi_sim_comm.Get_size()<2) log.printf("WARNING: ENSEMBLE with one replica is not doing any averaging!\n");
    ens_dim=multi_sim_comm.Get_size();
    my_repl=multi_sim_comm.Get_rank();
  } else { 
    ens_dim=0; 
    my_repl=0;
  }
  comm.Sum(&ens_dim, 1);
  comm.Sum(&my_repl, 1);
  
  // prepare output components, the number depending on reweighing or not
  narg = getNumberOfArguments();
  if(do_reweight) narg--;
  
  // these are the averages
  for(unsigned i=0;i<narg;i++) {
     std::string s=getPntrToArgument(i)->getName();
     addComponentWithDerivatives(s); 
     getPntrToComponent(i)->setNotPeriodic();
  }
  // these are the variances
  if(do_variance) {
    for(unsigned i=0;i<narg;i++) {
      std::string s=getPntrToArgument(i)->getName()+"_v";
      addComponentWithDerivatives(s); 
      getPntrToComponent(i+narg)->setNotPeriodic();
    }
  }
  log.printf("  using %u replicas.\n", ens_dim);
  if(do_reweight) log.printf("  doing simple REWEIGHT using the latest ARGUMENT as energy.\n");
  checkRead();

  // prepare vector for biases, means, variances, and derivatives
  bias.resize(ens_dim);
  mean.resize(narg);
  var.resize(narg);
  var_tmp.resize(narg);
}

void Ensemble::calculate(){

  double norm = 0.0;

  // in case of reweight
  if(do_reweight){
    // put bias to zero
    for(unsigned i=0; i<ens_dim; ++i) bias[i] = 0.0;
    if(master){
      // first we share the bias - the narg
      bias[my_repl] = getArgument(narg); 
      if(ens_dim>1) multi_sim_comm.Sum(&bias[0], ens_dim);  
    }
    // inside each replica
    comm.Sum(&bias[0], ens_dim);

    // find maximum value of bias
    double maxbias = *(std::max_element(bias.begin(), bias.end()));

    // calculate weights
    for(unsigned i=0; i<ens_dim; ++i){
       bias[i] = exp((bias[i]-maxbias)/kbt); 
       norm += bias[i];
    }
  // otherwise set weight to 1 and norm to ens_dim  
  } else {
    for(unsigned i=0; i<ens_dim; ++i){
       bias[i] = 1.0; 
       norm += bias[i];
    }   
  }

  double  w = bias[my_repl];
  double  fact = w/norm;
  double  fact_kbt = fact/kbt;

  // 1) calculate means
  // cycle on number of arguments - bias excluded
  if(master) for(unsigned i=0;i<narg;++i) mean[i] = getArgument(i) * fact;
  else       for(unsigned i=0;i<narg;++i) mean[i] = 0.;
  // among replicas
  if(master&&ens_dim>1) multi_sim_comm.Sum(&mean[0], narg);
  // inside each replica
  comm.Sum(&mean[0], narg);
  
  // 2) calculate variances
  // cycle on number of arguments - bias excluded
  if(do_variance) {
    if(master) { 
      for(unsigned i=0;i<narg;++i){
        var_tmp[i] = ( getArgument(i) - mean[i] ) * fact;
        var[i]     = ( getArgument(i) - mean[i] ) * var_tmp[i];
      }
    } else {
      for(unsigned i=0;i<narg;++i){
        var_tmp[i] = 0.; 
        var[i]     = 0.;
      }
    }
    // among replicas
    if(master&&ens_dim>1){
      multi_sim_comm.Sum(&var[0], narg);
      multi_sim_comm.Sum(&var_tmp[0], narg);
    }
    // inside each replica
    comm.Sum(&var[0], narg);
    comm.Sum(&var_tmp[0], narg);
  }
  
  // 3) set components
  for(unsigned i=0;i<narg;++i){
    // set mean
    Value* v=getPntrToComponent(i);
    v->set(mean[i]);
    setDerivative(v, i, fact);
    // useful tmp quantity
    const double der_tmp = getArgument(i) - mean[i];

    if(do_variance) {
      // set variance
      Value* u=getPntrToComponent(i+narg);
      u->set(var[i]);
      const double der = 2.0 * fact * ( der_tmp - var_tmp[i] );
      setDerivative(u, i, der);
      // if reweighing, derivative also wrt to bias
      if(do_reweight){
        const double der = ( der_tmp * der_tmp -2.0 * der_tmp * var_tmp[i] - var[i] ) * fact_kbt;
        setDerivative(u, narg, der);
      }
    }
    // if reweighing, derivative also wrt to bias
    if(do_reweight){
      // of the mean
      setDerivative(v, narg, der_tmp * fact_kbt);
    }
  };
}

}
}


