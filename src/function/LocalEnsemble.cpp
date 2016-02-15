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

using namespace std;

namespace PLMD{
namespace function{

//+PLUMEDOC FUNCTION LOCALENSEMBLE
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


class LocalEnsemble :
  public Function
{
  unsigned ens_dim;
  unsigned narg;
public:
  explicit LocalEnsemble(const ActionOptions&);
  void     calculate();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(LocalEnsemble,"LOCALENSEMBLE")

void LocalEnsemble::registerKeywords(Keywords& keys){
  Function::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","NUM","the number of local replicas");
  ActionWithValue::useCustomisableComponents(keys);
}

LocalEnsemble::LocalEnsemble(const ActionOptions&ao):
Action(ao),
Function(ao),
ens_dim(0)
{
  parse("NUM",ens_dim);
  if(ens_dim==0) error("NUM should be greater or equal to 1");

  vector<Value*> arg;
  int oldsize=-1;
  for(unsigned i=1;i<=ens_dim;++i ){
    vector<Value*> larg;
    if(!parseArgumentList("ARG",i,larg)) break;
    for(unsigned j=0;j<larg.size();j++) arg.push_back(larg[j]);
    if(oldsize!=-1&&oldsize!=larg.size()) error("In LOCALENSEMBLE you should have the same number of arguments for each ARG keyword");
    oldsize = larg.size();
    if(!larg.empty()){
      log.printf("  with arguments %u: ", i);
      for(unsigned j=0;j<larg.size();j++) log.printf(" %s",larg[j]->getName().c_str());
      log.printf("\n");
    }
  }
  requestArguments(arg);
  narg = arg.size()/ens_dim;

  // these are the averages
  for(unsigned i=0;i<narg;i++) {
    std::string s=getPntrToArgument(i)->getName();
    addComponentWithDerivatives(s); 
    getPntrToComponent(i)->setNotPeriodic();
  }

  log.printf("  averaging over %u replicas.\n", ens_dim);
}

void LocalEnsemble::calculate(){
  const double norm = static_cast<double>(ens_dim); 
  const double fact = 1.0/norm; 

  vector<double> mean(narg);
  vector<double> dmean(narg,fact);
  // calculate the mean 
  for(unsigned i=0;i<narg;++i) for(unsigned j=0;j<ens_dim;++j) mean[i] += fact*getArgument(j*narg+i); 

  // set components
  for(unsigned i=0;i<narg;++i){
    // set mean
    Value* v=getPntrToComponent(i);
    v->set(mean[i]);
    setDerivative(v, i, dmean[i]);
  } 
}

}
}


