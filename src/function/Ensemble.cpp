/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
  std::vector<double> cv;
  bool     do_reweight;
  bool     master;
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
  keys.add("optional","TEMP","the system temperature - this is only needed if you are reweighting");
  ActionWithValue::useCustomisableComponents(keys);
}

Ensemble::Ensemble(const ActionOptions&ao):
Action(ao),
Function(ao),
do_reweight(false),
kbt(-1.0)
{
  parseFlag("REWEIGHT", do_reweight); 
  double temp=0.0;
  parse("TEMP",temp);
  if(do_reweight) {
    if(temp>0.0) kbt=plumed.getAtoms().getKBoltzmann()*temp;
    else kbt=plumed.getAtoms().getKbT();
    if(kbt==0.0) error("Unless the MD engine passes the temperature to plumed, with REWEIGHT you must specify TEMP");
  }

  master = (comm.Get_rank()==0);
  if(master) {
    if(multi_sim_comm.Get_size()<2) error("You CANNOT run Replica-Averaged simulations without running multiple replicas!\n");
    else {ens_dim=multi_sim_comm.Get_size(); my_repl=multi_sim_comm.Get_rank();}
  } else {ens_dim=0; my_repl=0;}
  comm.Sum(&ens_dim, 1);
  comm.Sum(&my_repl, 1);
  
  // prepare output components, the number depending on reweighing or not
  narg = getNumberOfArguments();
  if(do_reweight) narg--;
  
  for(unsigned i=0;i<narg;i++) {
     std::string s=getPntrToArgument(i)->getName();
     addComponentWithDerivatives(s); 
     getPntrToComponent(i)->setNotPeriodic();
  }
  log.printf("  using %u replicas.\n", ens_dim);

  checkRead();

  // prepare vector for biases and cvs
  for(unsigned i=0; i<ens_dim; ++i) bias.push_back(0.0);
  for(unsigned i=0; i<narg; ++i) cv.push_back(0.0);
}

void Ensemble::calculate(){

  // put bias to zero
  for(unsigned i=0; i<ens_dim; ++i) bias[i] = 0.0;
  double norm = 0.0;

  // in case of reweight
  if(do_reweight){
    // first we share the bias - the narg
    double b = getArgument(narg);
    if(master){
      bias[my_repl] = b;
      multi_sim_comm.Sum(&bias[0], ens_dim);  
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
  // cycle on number of arguments - bias excluded
  for(unsigned i=0;i<narg;++i){
    if(master) cv[i] = getArgument(i) * fact;
    else cv[i] = 0.;
  }
  // among replicas
  if(master) multi_sim_comm.Sum(&cv[0], narg);
  // inside each replica
  comm.Sum(&cv[0], narg);
  for(unsigned i=0;i<narg;++i){
    Value* v=getPntrToComponent(i);
    v->set(cv[i]);
    setDerivative(v, i, fact);
    // if reweighing, derivative also wrt to bias
    if(do_reweight){
      double der = (getArgument(i) - cv[i]) * fact_kbt;
      setDerivative(v, narg, der);
    }
  };
}

}
}


