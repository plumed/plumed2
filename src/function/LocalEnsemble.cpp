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
#include "Function.h"
#include "ActionRegister.h"
#include "tools/OpenMP.h"

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION LOCALENSEMBLE
/*
Calculates the average over multiple arguments.

If more than one collective variable is given for each argument then they
are averaged separately. The average is stored in a component labelled <em>label</em>.cvlabel.

\par Examples

The following input tells plumed to calculate the chemical shifts for four
different proteins in the same simulation box then average them, calculated
the sum of the squared deviation with respect to the experimental values and
applies a linear restraint.
\plumedfile
MOLINFO STRUCTURE=data/template.pdb

chaina: GROUP ATOMS=1-1640
chainb: GROUP ATOMS=1641-3280
chainc: GROUP ATOMS=3281-4920
chaind: GROUP ATOMS=4921-6560

WHOLEMOLECULES ENTITY0=chaina ENTITY1=chainb ENTITY2=chainc ENTITY3=chaind

csa: CS2BACKBONE ATOMS=chaina NRES=100 DATA=data/ TEMPLATE=chaina.pdb NOPBC
csb: CS2BACKBONE ATOMS=chainb NRES=100 DATA=data/ TEMPLATE=chainb.pdb NOPBC
csc: CS2BACKBONE ATOMS=chainc NRES=100 DATA=data/ TEMPLATE=chainc.pdb NOPBC
csd: CS2BACKBONE ATOMS=chaind NRES=100 DATA=data/ TEMPLATE=chaind.pdb NOPBC

ensca: LOCALENSEMBLE NUM=4 ARG1=(csa\.ca_.*) ARG2=(csb\.ca_.*) ARG3=(csc\.ca_.*) ARG4=(csd\.ca_.*)
enscb: LOCALENSEMBLE NUM=4 ARG1=(csa\.cb_.*) ARG2=(csb\.cb_.*) ARG3=(csc\.cb_.*) ARG4=(csd\.cb_.*)
ensco: LOCALENSEMBLE NUM=4 ARG1=(csa\.co_.*) ARG2=(csb\.co_.*) ARG3=(csc\.co_.*) ARG4=(csd\.co_.*)
enshn: LOCALENSEMBLE NUM=4 ARG1=(csa\.hn_.*) ARG2=(csb\.hn_.*) ARG3=(csc\.hn_.*) ARG4=(csd\.hn_.*)
ensnh: LOCALENSEMBLE NUM=4 ARG1=(csa\.nh_.*) ARG2=(csb\.nh_.*) ARG3=(csc\.nh_.*) ARG4=(csd\.nh_.*)

stca: STATS ARG=(ensca\.csa\.ca_.*) PARARG=(csa\.expca_.*) SQDEVSUM
stcb: STATS ARG=(enscb\.csa\.cb_.*) PARARG=(csa\.expcb_.*) SQDEVSUM
stco: STATS ARG=(ensco\.csa\.co_.*) PARARG=(csa\.expco_.*) SQDEVSUM
sthn: STATS ARG=(enshn\.csa\.hn_.*) PARARG=(csa\.exphn_.*) SQDEVSUM
stnh: STATS ARG=(ensnh\.csa\.nh_.*) PARARG=(csa\.expnh_.*) SQDEVSUM

res: RESTRAINT ARG=stca.*,stcb.*,stco.*,sthn.*,stnh.* AT=0.,0.,0.,0.,0. KAPPA=0.,0.,0.,0.,0 SLOPE=16.,16.,12.,24.,0.5
\endplumedfile

*/
//+ENDPLUMEDOC


class LocalEnsemble :
  public Function
{
  unsigned ens_dim;
  unsigned narg;
public:
  explicit LocalEnsemble(const ActionOptions&);
  void     calculate() override;
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(LocalEnsemble,"LOCALENSEMBLE")

void LocalEnsemble::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","NUM","the number of local replicas");
  useCustomisableComponents(keys);
}

LocalEnsemble::LocalEnsemble(const ActionOptions&ao):
  Action(ao),
  Function(ao),
  ens_dim(0)
{
  parse("NUM",ens_dim);
  if(ens_dim==0) error("NUM should be greater or equal to 1");

  std::vector<Value*> arg;
  int oldsize=-1;
  for(unsigned i=1; i<=ens_dim; ++i ) {
    std::vector<Value*> larg;
    if(!parseArgumentList("ARG",i,larg)) break;
    for(unsigned j=0; j<larg.size(); j++) arg.push_back(larg[j]);
    if(oldsize!=-1&&oldsize!=static_cast<int>(larg.size())) error("In LOCALENSEMBLE you should have the same number of arguments for each ARG keyword");
    oldsize = larg.size();
    if(!larg.empty()) {
      log.printf("  with arguments %u: ", i);
      for(unsigned j=0; j<larg.size(); j++) log.printf(" %s",larg[j]->getName().c_str());
      log.printf("\n");
    }
  }
  requestArguments(arg);
  narg = arg.size()/ens_dim;

  // these are the averages
  for(unsigned i=0; i<narg; i++) {
    std::string s=getPntrToArgument(i)->getName();
    addComponentWithDerivatives(s);
    getPntrToComponent(i)->setNotPeriodic();
  }

  log.printf("  averaging over %u replicas.\n", ens_dim);
}

void LocalEnsemble::calculate()
{
  const double fact = 1.0/static_cast<double>(ens_dim);
  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for(unsigned i=0; i<narg; ++i) {
    double mean = 0.;
    Value* v=getPntrToComponent(i);
    for(unsigned j=0; j<ens_dim; ++j) {
      const unsigned index = j*narg+i;
      setDerivative(v, index, fact);
      mean += fact*getArgument(index);
    }
    v->set(mean);
  }
}

}
}


