/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2020 The plumed team
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

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION ENSEMBLE
/*
Calculates the replica averaging of a collective variable over multiple replicas.

Each collective variable is averaged separately and stored in a component labelled <em>label</em>.cvlabel.

\par Examples

The following input tells plumed to calculate the distance between atoms 3 and 5
and the average it over the available replicas.
\plumedfile
dist: DISTANCE ATOMS=3,5
ens: ENSEMBLE ARG=dist
PRINT ARG=dist,ens.dist
\endplumedfile

*/
//+ENDPLUMEDOC


class Ensemble :
  public Function
{
  unsigned ens_dim;
  unsigned my_repl;
  unsigned narg;
  bool     master;
  bool     do_reweight;
  bool     do_moments;
  bool     do_central;
  bool     do_powers;
  double   kbt;
  double   moment;
  double   power;
public:
  explicit Ensemble(const ActionOptions&);
  void     calculate() override;
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(Ensemble,"ENSEMBLE")

void Ensemble::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys);
  keys.use("ARG");
  keys.addFlag("REWEIGHT",false,"simple REWEIGHT using the latest ARG as energy");
  keys.addFlag("CENTRAL",false,"calculate a central moment instead of a standard moment");
  keys.add("optional","TEMP","the system temperature - this is only needed if you are reweighting");
  keys.add("optional","MOMENT","the moment you want to calculate in alternative to the mean or the variance");
  keys.add("optional","POWER","the power of the mean (and moment)");
  ActionWithValue::useCustomisableComponents(keys);
}

Ensemble::Ensemble(const ActionOptions&ao):
  Action(ao),
  Function(ao),
  do_reweight(false),
  do_moments(false),
  do_central(false),
  do_powers(false),
  kbt(-1.0),
  moment(0),
  power(0)
{
  parseFlag("REWEIGHT", do_reweight);
  double temp=0.0;
  parse("TEMP",temp);
  if(do_reweight) {
    if(temp>0.0) kbt=plumed.getAtoms().getKBoltzmann()*temp;
    else kbt=plumed.getAtoms().getKbT();
    if(kbt==0.0) error("Unless the MD engine passes the temperature to plumed, with REWEIGHT you must specify TEMP");
  }

  parse("MOMENT",moment);
  if(moment==1) error("MOMENT can be any number but for 0 and 1");
  if(moment!=0) do_moments=true;
  parseFlag("CENTRAL", do_central);
  if(!do_moments&&do_central) error("To calculate a CENTRAL moment you need to define for which MOMENT");

  parse("POWER",power);
  if(power==1) error("POWER can be any number but for 0 and 1");
  if(power!=0) do_powers=true;

  checkRead();

  master = (comm.Get_rank()==0);
  ens_dim=0;
  my_repl=0;
  if(master) {
    ens_dim=multi_sim_comm.Get_size();
    my_repl=multi_sim_comm.Get_rank();
  }
  comm.Bcast(ens_dim,0);
  comm.Bcast(my_repl,0);
  if(ens_dim<2) log.printf("WARNING: ENSEMBLE with one replica is not doing any averaging!\n");

  // prepare output components, the number depending on reweighing or not
  narg = getNumberOfArguments();
  if(do_reweight) narg--;

  // these are the averages
  for(unsigned i=0; i<narg; i++) {
    std::string s=getPntrToArgument(i)->getName();
    addComponentWithDerivatives(s);
    getPntrToComponent(i)->setNotPeriodic();
  }
  // these are the moments
  if(do_moments) {
    for(unsigned i=0; i<narg; i++) {
      std::string s=getPntrToArgument(i)->getName()+"_m";
      addComponentWithDerivatives(s);
      getPntrToComponent(i+narg)->setNotPeriodic();
    }
  }

  log.printf("  averaging over %u replicas.\n", ens_dim);
  if(do_reweight) log.printf("  doing simple REWEIGHT using the latest ARGUMENT as energy.\n");
  if(do_moments&&!do_central)  log.printf("  calculating also the %lf standard moment\n", moment);
  if(do_moments&&do_central)   log.printf("  calculating also the %lf central moment\n", moment);
  if(do_powers)                log.printf("  calculating the %lf power of the mean (and moment)\n", power);
}

void Ensemble::calculate() {
  double norm = 0.0;
  double fact = 0.0;

  // calculate the weights either from BIAS
  if(do_reweight) {
    std::vector<double> bias;
    bias.resize(ens_dim);
    if(master) {
      bias[my_repl] = getArgument(narg);
      if(ens_dim>1) multi_sim_comm.Sum(&bias[0], ens_dim);
    }
    comm.Sum(&bias[0], ens_dim);
    const double maxbias = *(std::max_element(bias.begin(), bias.end()));
    for(unsigned i=0; i<ens_dim; ++i) {
      bias[i] = exp((bias[i]-maxbias)/kbt);
      norm += bias[i];
    }
    fact = bias[my_repl]/norm;
    // or arithmetic ones
  } else {
    norm = static_cast<double>(ens_dim);
    fact = 1.0/norm;
  }

  const double fact_kbt = fact/kbt;

  std::vector<double> mean(narg);
  std::vector<double> dmean(narg,fact);
  // calculate the mean
  if(master) {
    for(unsigned i=0; i<narg; ++i) mean[i] = fact*getArgument(i);
    if(ens_dim>1) multi_sim_comm.Sum(&mean[0], narg);
  }
  comm.Sum(&mean[0], narg);

  std::vector<double> v_moment, dv_moment;
  // calculate other moments
  if(do_moments) {
    v_moment.resize(narg);
    dv_moment.resize(narg);
    // standard moment
    if(!do_central) {
      if(master) {
        for(unsigned i=0; i<narg; ++i) {
          const double tmp = fact*std::pow(getArgument(i),moment-1);
          v_moment[i]      = tmp*getArgument(i);
          dv_moment[i]     = moment*tmp;
        }
        if(ens_dim>1) multi_sim_comm.Sum(&v_moment[0], narg);
      } else {
        for(unsigned i=0; i<narg; ++i) {
          const double tmp = fact*std::pow(getArgument(i),moment-1);
          dv_moment[i]     = moment*tmp;
        }
      }
      // central moment
    } else {
      if(master) {
        for(unsigned i=0; i<narg; ++i) {
          const double tmp = std::pow(getArgument(i)-mean[i],moment-1);
          v_moment[i]      = fact*tmp*(getArgument(i)-mean[i]);
          dv_moment[i]     = moment*tmp*(fact-fact/norm);
        }
        if(ens_dim>1) multi_sim_comm.Sum(&v_moment[0], narg);
      } else {
        for(unsigned i=0; i<narg; ++i) {
          const double tmp = std::pow(getArgument(i)-mean[i],moment-1);
          dv_moment[i]     = moment*tmp*(fact-fact/norm);
        }
      }
    }
    comm.Sum(&v_moment[0], narg);
  }

  // calculate powers of moments
  if(do_powers) {
    for(unsigned i=0; i<narg; ++i) {
      const double tmp1 = std::pow(mean[i],power-1);
      mean[i]          *= tmp1;
      dmean[i]         *= power*tmp1;
      if(do_moments) {
        const double tmp2 = std::pow(v_moment[i],power-1);
        v_moment[i]      *= tmp2;
        dv_moment[i]     *= power*tmp2;
      }
    }
  }

  // set components
  for(unsigned i=0; i<narg; ++i) {
    // set mean
    Value* v=getPntrToComponent(i);
    v->set(mean[i]);
    setDerivative(v, i, dmean[i]);
    if(do_reweight) {
      const double w_tmp = fact_kbt*(getArgument(i) - mean[i]);
      setDerivative(v, narg, w_tmp);
    }
    if(do_moments) {
      // set moments
      Value* u=getPntrToComponent(i+narg);
      u->set(v_moment[i]);
      setDerivative(u, i, dv_moment[i]);
      if(do_reweight) {
        const double w_tmp = fact_kbt*(pow(getArgument(i),moment) - v_moment[i]);
        setDerivative(u, narg, w_tmp);
      }
    }
  }
}

}
}
