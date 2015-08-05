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
  unsigned  ens_dim;
  double    fact;
public:
  explicit Ensemble(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(Ensemble,"ENSEMBLE")

void Ensemble::registerKeywords(Keywords& keys){
  Function::registerKeywords(keys);
  keys.use("ARG");
  ActionWithValue::useCustomisableComponents(keys);
}

Ensemble::Ensemble(const ActionOptions&ao):
Action(ao),
Function(ao)
{
  if(comm.Get_rank()==0) {
    if(multi_sim_comm.Get_size()<2) error("You CANNOT run Replica-Averaged simulations without running multiple replicas!\n");
    else ens_dim=multi_sim_comm.Get_size(); 
  } else ens_dim=0; 
  comm.Sum(&ens_dim, 1);
  fact = 1./((double) ens_dim);
  for(unsigned i=0;i<getNumberOfArguments();i++) {
     std::string s=getPntrToArgument(i)->getName();
     addComponentWithDerivatives(s); 
     getPntrToComponent(i)->setNotPeriodic();
  }
  log.printf("  using %u replicas.\n", ens_dim);
  checkRead();
}

void Ensemble::calculate(){
  for(unsigned i=0;i<getNumberOfArguments();++i){
    double cv=0.;
    if(comm.Get_rank()==0) { // I am the master of my replica
      // among replicas
      cv = getArgument(i);
      multi_sim_comm.Sum(&cv, 1);
      cv *= fact; 
    }
    // inside each replica
    comm.Sum(&cv, 1);
    Value* v=getPntrToComponent(i);
    v->set(cv);
    setDerivative(v,i,fact);
  };
}

}
}


