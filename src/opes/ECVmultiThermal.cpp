/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2020-2021 of Michele Invernizzi.

   This file is part of the OPES plumed module.

   The OPES plumed module is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The OPES plumed module is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "ExpansionCVs.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace opes {

//+PLUMEDOC OPES_EXPANSION_CV ECV_MULTITHERMAL
/*
Expand a simulation to sample multiple temperatures simultaneously.

The internal energy $U$ of of the system should be used as ARG.

$$
  \Delta u_{\beta'}=(\beta'-\beta) U\, ,
$$

where $\beta'$ are the temperatures to be sampled and $\beta$ is the temperature at which the simulation is conducted.
In case of fixed volume, the internal energy is simply the potential energy given by the [ENERGY](ENERGY.md) colvar $U=E$, and you will run a multicanonical simulation.
If instead the simulation is at fixed pressure $p$, the contribution of the volume must be added $U=E+pV$ (see example below).

By defauly the needed steps in temperatures are automatically guessed from few initial unbiased MD steps, as descibed in the paper cited below.
Otherwise you can manually set this number with TEMP_STEPS.
In both cases the steps will be geometrically spaced in temperature.
Use instead the keyword NO_GEOM_SPACING for a linear spacing in the inverse temperature (beta), that typically increases the focus on lower temperatures.
Finally, you can use instead the keyword TEMP_SET_ALL and explicitly provide each temperature.

You can reweight the resulting simulation at any temperature in the chosen range, using e.g. [REWEIGHT_TEMP_PRESS](REWEIGHT_TEMP_PRESS.md).
A similar target distribution can be sampled using [TD_MULTICANONICAL](TD_MULTICANONICAL.md).

## Examples

Fixed volume, multicanonical simulation:

```plumed
ene: ENERGY
ecv: ECV_MULTITHERMAL ARG=ene TEMP=300 TEMP_MIN=300 TEMP_MAX=800
opes: OPES_EXPANDED ARG=ecv.ene PACE=500
```

which, if your MD code passes the temperature to PLUMED, is equivalent to:

```plumed
ene: ENERGY
ecv: ECV_MULTITHERMAL ARG=ene TEMP_MAX=800
opes: OPES_EXPANDED ARG=ecv.ene PACE=500
```

If instead the pressure is fixed and the volume changes, you shuld calculate the internal energy first, $U=E+pV$

```plumed
ene: ENERGY
vol: VOLUME
intEne: CUSTOM PERIODIC=NO ARG=ene,vol FUNC=x+0.06022140857*y
ecv: ECV_MULTITHERMAL ARG=intEne TEMP_MAX=800
opes: OPES_EXPANDED ARG=ecv.intEne PACE=500
```

Notice that $p=0.06022140857$ corresponds to 1 bar when using the default PLUMED units.

*/
//+ENDPLUMEDOC

class ECVmultiThermal :
  public ExpansionCVs {
private:
  bool todoAutomatic_;
  bool geom_spacing_;
  std::vector<double> ECVs_;
  std::vector<double> derECVs_; //(beta_k-beta0) or (temp0/temp_k-1)/kbt
  void initECVs();

public:
  explicit ECVmultiThermal(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  void calculateECVs(const double *) override;
  const double * getPntrToECVs(unsigned) override;
  const double * getPntrToDerECVs(unsigned) override;
  std::vector<std::string> getLambdas() const override;
  void initECVs_observ(const std::vector<double>&,const unsigned,const unsigned) override;
  void initECVs_restart(const std::vector<std::string>&) override;
};

PLUMED_REGISTER_ACTION(ECVmultiThermal,"ECV_MULTITHERMAL")

void ECVmultiThermal::registerKeywords(Keywords& keys) {
  ExpansionCVs::registerKeywords(keys);
  keys.addInputKeyword("compulsory","ARG","scalar","the label of the internal energy of the system. If volume is fixed it is calculated by the ENERGY colvar");
  keys.add("optional","TEMP_MIN","the minimum of the temperature range");
  keys.add("optional","TEMP_MAX","the maximum of the temperature range");
  keys.add("optional","TEMP_STEPS","the number of steps in temperature");
  keys.add("optional","TEMP_SET_ALL","manually set all the temperatures");
  keys.addFlag("NO_GEOM_SPACING",false,"do not use geometrical spacing in temperature, but instead linear spacing in inverse temperature");
  keys.addDOI("10.1103/PhysRevX.10.041034");
}

ECVmultiThermal::ECVmultiThermal(const ActionOptions&ao)
  : Action(ao)
  , ExpansionCVs(ao)
  , todoAutomatic_(false) {
  plumed_massert(getNumberOfArguments()==1,"only the internal energy should be given as ARG");

//set temp0
  const double temp0=kbt_/getKBoltzmann();

//parse temp range
  double temp_min=-1;
  double temp_max=-1;
  parse("TEMP_MIN",temp_min);
  parse("TEMP_MAX",temp_max);
  unsigned temp_steps=0;
  parse("TEMP_STEPS",temp_steps);
  std::vector<double> temps;
  parseVector("TEMP_SET_ALL",temps);
  parseFlag("NO_GEOM_SPACING",geom_spacing_);
  geom_spacing_=!geom_spacing_;

  checkRead();

//set the intermediate temperatures
  if(temps.size()>0) {
    plumed_massert(temp_steps==0,"cannot set both TEMP_STEPS and TEMP_SET_ALL");
    plumed_massert(temp_min==-1 && temp_max==-1,"cannot set both TEMP_SET_ALL and TEMP_MIN/MAX");
    plumed_massert(temps.size()>=2,"set at least 2 temperatures");
    temp_min=temps[0];
    temp_max=temps[temps.size()-1];
    derECVs_.resize(temps.size());
    for(unsigned k=0; k<derECVs_.size(); k++) {
      derECVs_[k]=(temp0/temps[k]-1.)/kbt_;
      if(k<derECVs_.size()-1) {
        plumed_massert(temps[k]<=temps[k+1],"TEMP_SET_ALL must be properly ordered");
      }
    }
  } else {
    //get TEMP_MIN and TEMP_MAX
    plumed_massert(temp_min!=-1 || temp_max!=-1,"TEMP_MIN, TEMP_MAX or both, should be set");
    if(temp_min==-1) {
      temp_min=temp0;
      log.printf("  no TEMP_MIN provided, using TEMP_MIN=TEMP\n");
    }
    if(temp_max==-1) {
      temp_max=temp0;
      log.printf("  no TEMP_MAX provided, using TEMP_MAX=TEMP\n");
    }
    plumed_massert(temp_max>=temp_min,"TEMP_MAX should be bigger than TEMP_MIN");
    derECVs_.resize(2);
    derECVs_[0]=(temp0/temp_min-1.)/kbt_;
    derECVs_[1]=(temp0/temp_max-1.)/kbt_;
    if(temp_min==temp_max && temp_steps==0) {
      temp_steps=1;
    }
    if(temp_steps>0) {
      derECVs_=getSteps(derECVs_[0],derECVs_[1],temp_steps,"TEMP",geom_spacing_,1./kbt_);
    } else {
      todoAutomatic_=true;
    }
  }
  const double tol=1e-3; //if temp is taken from MD engine it might be numerically slightly different
  if(temp0<(1-tol)*temp_min || temp0>(1+tol)*temp_max) {
    log.printf(" +++ WARNING +++ running at TEMP=%g which is outside the chosen temperature range\n",temp0);
  }

//print some info
  log.printf("  targeting a temperature range from TEMP_MIN=%g to TEMP_MAX=%g\n",temp_min,temp_max);
  if(!geom_spacing_) {
    log.printf(" -- NO_GEOM_SPACING: inverse temperatures will be linearly spaced\n");
  }
}

void ECVmultiThermal::calculateECVs(const double * ene) {
  for(unsigned k=0; k<derECVs_.size(); k++) {
    ECVs_[k]=derECVs_[k]*ene[0];
  }
// derivatives never change: derECVs_k=(beta_k-beta0)
}

const double * ECVmultiThermal::getPntrToECVs(unsigned j) {
  plumed_massert(isReady_,"cannot access ECVs before initialization");
  plumed_massert(j==0,getName()+" has only one CV, the ENERGY");
  return &ECVs_[0];
}

const double * ECVmultiThermal::getPntrToDerECVs(unsigned j) {
  plumed_massert(isReady_,"cannot access ECVs before initialization");
  plumed_massert(j==0,getName()+" has only one CV, the ENERGY");
  return &derECVs_[0];
}

std::vector<std::string> ECVmultiThermal::getLambdas() const {
  plumed_massert(!todoAutomatic_,"cannot access lambdas before initializing them");
  const double temp0=kbt_/getKBoltzmann();
  std::vector<std::string> lambdas(derECVs_.size());
  for(unsigned k=0; k<derECVs_.size(); k++) {
    std::ostringstream subs;
    subs<<temp0/(derECVs_[k]*kbt_+1); //temp_k
    lambdas[k]=subs.str();
  }
  return lambdas;
}

void ECVmultiThermal::initECVs() {
  plumed_massert(!isReady_,"initialization should not be called twice");
  plumed_massert(!todoAutomatic_,"this should not happen");
  totNumECVs_=derECVs_.size();
  ECVs_.resize(derECVs_.size());
  isReady_=true;
  log.printf("  *%4lu temperatures for %s\n",derECVs_.size(),getName().c_str());
}

void ECVmultiThermal::initECVs_observ(const std::vector<double>& all_obs_cvs,const unsigned ncv,const unsigned index_j) {
  if(todoAutomatic_) { //estimate the steps in beta from observations
    plumed_massert(all_obs_cvs.size()%ncv==0 && index_j<ncv,"initECVs_observ parameters are inconsistent");
    std::vector<double> obs_ene(all_obs_cvs.size()/ncv); //copy only useful observation (would be better not to copy...)
    for(unsigned t=0; t<obs_ene.size(); t++) {
      obs_ene[t]=all_obs_cvs[t*ncv+index_j];
    }
    const unsigned temp_steps=estimateNumSteps(derECVs_[0],derECVs_[1],obs_ene,"TEMP");
    log.printf("    (spacing is in beta, not in temperature)\n");
    derECVs_=getSteps(derECVs_[0],derECVs_[1],temp_steps,"TEMP",geom_spacing_,1./kbt_);
    todoAutomatic_=false;
  }
  initECVs();
  calculateECVs(&all_obs_cvs[index_j]);
}

void ECVmultiThermal::initECVs_restart(const std::vector<std::string>& lambdas) {
  std::size_t pos=lambdas[0].find("_");
  plumed_massert(pos==std::string::npos,"this should not happen, only one CV is used in "+getName());
  if(todoAutomatic_) {
    derECVs_=getSteps(derECVs_[0],derECVs_[1],lambdas.size(),"TEMP",geom_spacing_,1./kbt_);
    todoAutomatic_=false;
  }
  std::vector<std::string> myLambdas=getLambdas();
  plumed_massert(myLambdas.size()==lambdas.size(),"RESTART - mismatch in number of "+getName());
  plumed_massert(std::equal(myLambdas.begin(),myLambdas.end(),lambdas.begin()),"RESTART - mismatch in lambda values of "+getName());

  initECVs();
}

}
}
