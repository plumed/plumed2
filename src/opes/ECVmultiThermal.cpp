/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2020 of Michele Invernizzi.

The opes module is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The opes module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "ExpansionCVs.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

namespace PLMD {
namespace opes {

//+PLUMEDOC EXPANSION_CV ECV_MULTITHERMAL
/*
Expand a simulation to sample multiple temperatures simultaneously.

The internal energy \f$U\f$ of of the system should be used as ARG.
\f[
  \Delta u_{\beta'}=(\beta-\beta') U\, .
\f]
In case of fixed volume, the internal energy is simply the potential energy given by the \ref ENERGY colvar\f$U=E\f$, and you will run a multicanonical simulation.
If instead the simulation is at fixed pressure \f$p\f$, the contribution of the volume must be added \f$U=E+pV\f$ (see example below).

By defauly the needed steps in temperatures are automatically guessed from few initial unbiased MD steps.
Otherwise you can manually set this number with STEPS_TEMP.
In both cases the steps will have a geometric distriution in the temperature, thus uniform in beta.
For more fine controll you can use the keywork SET_ALL_TEMPS and explicitly provide each temperature.

You can reweight the resulting simulation at any temperature in the chosen range, using e.g. \ref REWEIGHT_TEMP_PRESS.
A similar target distribution can be sampled using \ref TD_MULTICANONICAL.


\par Examples

Fixed volume, multicanonical simulation:

\plumedfile
ene: ENERGY
ecv: ECV_MULTITHERMAL ARG=ene TEMP=300 MIN_TEMP=300 MAX_TEMP=800
opes: OPES_EXPANDED ARG=ecv.ene PACE=500
\endplumedfile

which, if your MD code passes the temperature to PLUMED, is equivalent to:

\plumedfile
ene: ENERGY
ecv: ECV_MULTITHERMAL ARG=ene MAX_TEMP=800
opes: OPES_EXPANDED ARG=ecv.ene PACE=500
\endplumedfile

If instead the pressure is fixed and the volume changes, you must calculate the internal energy first, \f$U=E+pV\f$

\plumedfile
ene: ENERGY
vol: VOLUME
intEne: CUSTOM PERIODIC=NO ARG=ene,vol FUNC=x+0.06022140857*y
ecv: ECV_MULTITHERMAL ARG=intEne MAX_TEMP=800
opes: OPES_EXPANDED ARG=ecv.intEne PACE=500
\endplumedfile

Notice that \f$p=0.06022140857\f$ corresponds to 1 bar when using the default PLUMED units.

*/
//+ENDPLUMEDOC

class ECVmultiCanonical :
  public ExpansionCVs
{
private:
  bool todoAutomatic_;
  double beta0_;
  std::vector<double> beta_;
  std::vector<double> ECVs_;
  std::vector<double> derECVs_;
  void initECVs();

public:
  explicit ECVmultiCanonical(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  void calculateECVs(const double *) override;
  const double * getPntrToECVs(unsigned) override;
  const double * getPntrToDerECVs(unsigned) override;
  std::vector<std::string> getLambdas() const override;
  void initECVs_observ(const std::vector<double>&,const unsigned,const unsigned) override;
  void initECVs_restart(const std::vector<std::string>&) override;
};

PLUMED_REGISTER_ACTION(ECVmultiCanonical,"ECV_MULTITHERMAL")

void ECVmultiCanonical::registerKeywords(Keywords& keys)
{
  ExpansionCVs::registerKeywords(keys);
  keys.remove("ARG");
  keys.add("compulsory","ARG","the label of the internal energy of the system. If volume is fixed it is calculated by the \\ref ENERGY colvar");
  keys.add("optional","MIN_TEMP","the minimum of the temperature range");
  keys.add("optional","MAX_TEMP","the maximum of the temperature range");
  keys.add("optional","STEPS_TEMP","the number of steps in temperature");
  keys.add("optional","SET_ALL_TEMPS","manually set all the temperatures");
}

ECVmultiCanonical::ECVmultiCanonical(const ActionOptions&ao)
  : Action(ao)
  , ExpansionCVs(ao)
  , todoAutomatic_(false)
  , beta0_(1./kbt_)
{
  plumed_massert(getNumberOfArguments()==1,"only the internal energy should be given as ARG");

//set temp0 and beta0_
  const double Kb=plumed.getAtoms().getKBoltzmann();
  double temp0=kbt_/Kb;

//parse temp range
  double min_temp=-1;
  double max_temp=-1;
  parse("MIN_TEMP",min_temp);
  parse("MAX_TEMP",max_temp);
  unsigned steps_temp=0;
  parse("STEPS_TEMP",steps_temp);
  std::vector<double> temps;
  parseVector("SET_ALL_TEMPS",temps);

  checkRead();

//set the intermediate temperatures
  if(temps.size()>0)
  {
    plumed_massert(steps_temp==0,"cannot set both STEPS_TEMP and SET_ALL_TEMPS");
    plumed_massert(min_temp==-1 && max_temp==-1,"cannot set both SET_ALL_TEMPS and MIN/MAX_TEMP");
    plumed_massert(temps.size()>=2,"set at least 2 temperatures");
    min_temp=temps[0];
    max_temp=temps[temps.size()-1];
    beta_.resize(temps.size());
    for(unsigned k=0; k<beta_.size(); k++)
    {
      beta_[k]=1./(Kb*temps[k]);
      if(k<beta_.size()-1)
        plumed_massert(temps[k]<=temps[k+1],"SET_ALL_TEMPS must be properly ordered");
    }
  }
  else
  { //get MIN_TEMP and MAX_TEMP
    plumed_massert(min_temp!=-1 || max_temp!=-1,"MIN_TEMP, MAX_TEMP or both, should be set");
    if(min_temp==-1)
    {
      min_temp=temp0;
      log.printf("  no MIN_TEMP provided, using MIN_TEMP=TEMP\n");
    }
    if(max_temp==-1)
    {
      max_temp=temp0;
      log.printf("  no MAX_TEMP provided, using MAX_TEMP=TEMP\n");
    }
    plumed_massert(max_temp>=min_temp,"MAX_TEMP should be bigger than MIN_TEMP");
    beta_.resize(2);
    beta_[0]=1./(Kb*min_temp); //ordered temp, inverted beta
    beta_[1]=1./(Kb*max_temp);
    if(min_temp==max_temp && steps_temp==0)
      steps_temp=1;
    if(steps_temp>0)
      setSteps(beta_,steps_temp,"TEMP");
    else
      todoAutomatic_=true;
  }
  const double tol=1e-3; //if temp is taken from MD engine it might be numerically slightly different
  if(temp0<(1-tol)*min_temp || temp0>(1+tol)*max_temp)
    log.printf(" +++ WARNING +++ running at TEMP=%g which is outside the chosen temperature range\n",temp0);

//print some info
  log.printf("  targeting a temperature range from MIN_TEMP=%g to MAX_TEMP=%g\n",min_temp,max_temp);
}

void ECVmultiCanonical::calculateECVs(const double * ene)
{
  for(unsigned k=0; k<beta_.size(); k++)
  {
    const double diff_k=(beta_[k]-beta0_);
    ECVs_[k]=diff_k*ene[0];
    derECVs_[k]=diff_k;
  }
}

const double * ECVmultiCanonical::getPntrToECVs(unsigned j)
{
  plumed_massert(isReady_,"cannot access ECVs before initialization");
  plumed_massert(j==0,getName()+" has only one CV, the ENERGY");
  return &ECVs_[0];
}

const double * ECVmultiCanonical::getPntrToDerECVs(unsigned j)
{
  plumed_massert(isReady_,"cannot access ECVs before initialization");
  plumed_massert(j==0,getName()+" has only one CV, the ENERGY");
  return &derECVs_[0];
}

std::vector<std::string> ECVmultiCanonical::getLambdas() const
{
  plumed_massert(!todoAutomatic_,"cannot access lambdas before initializing them");
  std::vector<std::string> lambdas(beta_.size());
  const double Kb=plumed.getAtoms().getKBoltzmann();
  for(unsigned k=0; k<beta_.size(); k++)
  {
    std::ostringstream subs;
    subs<<1./(Kb*beta_[k]);
    lambdas[k]=subs.str();
  }
  return lambdas;
}

void ECVmultiCanonical::initECVs()
{
  plumed_massert(!isReady_,"initialization should not be called twice");
  plumed_massert(!todoAutomatic_,"this should not happen");
  totNumECVs_=beta_.size();
  ECVs_.resize(beta_.size());
  derECVs_.resize(beta_.size());
  isReady_=true;
  log.printf("  *%4lu temperatures for %s\n",beta_.size(),getName().c_str());
}

void ECVmultiCanonical::initECVs_observ(const std::vector<double>& all_obs_cvs,const unsigned ncv,const unsigned index_j)
{
  if(todoAutomatic_) //estimate the steps in beta from observations
  {
    plumed_massert(all_obs_cvs.size()%ncv==0 && index_j<ncv,"initECVs_observ parameters are inconsistent");
    std::vector<double> obs_ene(all_obs_cvs.size()/ncv); //copy only useful observation (would be better not to copy...)
    for(unsigned t=0; t<obs_ene.size(); t++)
      obs_ene[t]=all_obs_cvs[t*ncv+index_j];
    const unsigned steps_temp=estimateSteps(beta_[0]-beta0_,beta_[1]-beta0_,obs_ene,"TEMP");
    log.printf("    (spacing is in beta, not in temperature)\n");
    setSteps(beta_,steps_temp,"TEMP");
    todoAutomatic_=false;
  }
  initECVs();
  calculateECVs(&all_obs_cvs[index_j]);
}

void ECVmultiCanonical::initECVs_restart(const std::vector<std::string>& lambdas)
{
  std::size_t pos=lambdas[0].find("_");
  plumed_massert(pos==std::string::npos,"this should not happen, only one CV is used in "+getName());
  if(todoAutomatic_)
  {
    setSteps(beta_,lambdas.size(),"TEMP");
    todoAutomatic_=false;
  }
  std::vector<std::string> myLambdas=getLambdas();
  plumed_massert(myLambdas.size()==lambdas.size(),"RESTART - mismatch in number of "+getName());
  plumed_massert(std::equal(myLambdas.begin(),myLambdas.end(),lambdas.begin()),"RESTART - mismatch in lambda values of "+getName());

  initECVs();
}

}
}
