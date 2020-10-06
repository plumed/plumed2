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

//+PLUMEDOC EXPANSION_CV ECV_MULTICANONICAL
/*
Expand a canonical simulation to sample multiple temperatures.
If instead of fixed volume NVT you are running with fixed pressure NPT, you must use \ref ECV_MULTITHERMAL_MULTIBARIC and add the volume contribution.
The \ref ENERGY of the system should be used as ARG.

\par Examples

mc: ECV_MULTICANONICAL ARG=ene TEMP=300 MIN_TEMP=300 MAX_TEMP=500

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
  void setBetaSteps(double,double,unsigned);
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

PLUMED_REGISTER_ACTION(ECVmultiCanonical,"ECV_MULTICANONICAL")

void ECVmultiCanonical::registerKeywords(Keywords& keys) {
  ExpansionCVs::registerKeywords(keys);
  keys.remove("ARG");
  keys.add("compulsory","ARG","provide the label of the potential energy of the system, you can calculate it with the \\ref ENERGY colvar");
  keys.add("compulsory","TEMP","-1","temperature. If not specified tries to get it from MD engine");
  keys.add("optional","MIN_TEMP","the minimum of the temperature range");
  keys.add("optional","MAX_TEMP","the maximum of the temperature range");
  keys.add("optional","STEPS_TEMP","the number of steps in temperature");
  keys.add("optional","SET_ALL_TEMPS","manually set all the temperatures");
//  keys.add("optional","BORDER_WEIGHT","set it greater than 1 to obtain better sampling of the max and min thermodynamics conditions");
}

ECVmultiCanonical::ECVmultiCanonical(const ActionOptions&ao):
  Action(ao),
  ExpansionCVs(ao)
{
  plumed_massert(getNumberOfArguments()==1,"only ENERGY should be given as ARG");

//set Kb and beta0_
  const double Kb=plumed.getAtoms().getKBoltzmann();
  double KbT=plumed.getAtoms().getKbT();
  double temp=-1;
  parse("TEMP",temp);
  if(temp>0)
  {
    if(KbT>0 && std::abs(KbT-Kb*temp)>1e-4)
      log.printf(" +++ WARNING +++ using TEMP=%g while MD engine uses %g\n",temp,KbT/Kb);
    KbT=Kb*temp;
  }
  plumed_massert(KbT>0,"your MD engine does not pass the temperature to plumed, you must specify it using TEMP");
  beta0_=1./KbT;
  temp=KbT/Kb;

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
    todoAutomatic_=false;
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
        plumed_massert(temps[k]<temps[k+1],"SET_ALL_TEMPS must be properly ordered");
    }
  }
  else
  { //get MIN_TEMP and MAX_TEMP
    plumed_massert(min_temp!=-1 || max_temp!=-1,"MIN_TEMP, MAX_TEMP or both, should be set");
    if(min_temp==-1)
    {
      min_temp=temp;
      log.printf("  no MIN_TEMP provided, using MIN_TEMP=TEMP\n");
    }
    if(max_temp==-1)
    {
      max_temp=temp;
      log.printf("  no MAX_TEMP provided, using MAX_TEMP=TEMP\n");
    }
    plumed_massert(max_temp>=min_temp,"MAX_TEMP should be bigger than MIN_TEMP");
    if(min_temp==max_temp && steps_temp==0)
      steps_temp=1;
    const double min_beta=1./(Kb*max_temp);
    const double max_beta=1./(Kb*min_temp);
    if(steps_temp>0)
      setBetaSteps(min_beta,max_beta,steps_temp);
    else
    {
      todoAutomatic_=true;
      beta_.resize(2);
      beta_[0]=min_beta;
      beta_[1]=max_beta;
    }
  }
  const double tol=1e-3; //if temp is taken from MD engine it might be numerically slightly different
  if(temp<(1-tol)*min_temp || temp>(1+tol)*max_temp)
    log.printf(" +++ WARNING +++ running at TEMP=%g which is outside the chosen temperature range\n",temp);

//print some info
  log.printf("  running at TEMP=%g\n",temp);
  log.printf("  targeting a temperature range from MIN_TEMP=%g to MAX_TEMP=%g\n",min_temp,max_temp);
}

void ECVmultiCanonical::calculateECVs(const double * ene) {
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
  plumed_massert(j==0,"ECV_MULTICANONICAL has only one CV, the ENERGY");
  return &ECVs_[0];
}

const double * ECVmultiCanonical::getPntrToDerECVs(unsigned j)
{
  plumed_massert(isReady_,"cannot access ECVs before initialization");
  plumed_massert(j==0,"ECV_MULTICANONICAL has only one CV, the ENERGY");
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

void ECVmultiCanonical::setBetaSteps(double min_beta,double max_beta,unsigned steps_beta)
{
  plumed_massert(beta_.size()==0 || beta_.size()==2,"you should not set the beta steps twice...");
  plumed_massert(min_beta<=max_beta,"this should not happen");
  plumed_massert(!(min_beta==max_beta && steps_beta>1),"cannot have multiple STEPS_TEMP if MIN_TEMP==MAX_TEMP");
  beta_.resize(steps_beta);
  if(steps_beta==1)
  {
    beta_[0]=(min_beta+max_beta)/2.;
    log.printf(" +++ WARNING +++ using one single temperature as target, corresponding to beta = %g\n",beta_[0]);
  }
  else
    for(unsigned k=0; k<beta_.size(); k++)//betas are stored in reversed order, so temps are in correct order
      beta_[k]=max_beta-k*(max_beta-min_beta)/(steps_beta-1);
  todoAutomatic_=false;
}

void ECVmultiCanonical::initECVs()
{
  plumed_massert(!isReady_,"initialization should not be called twice");
  plumed_massert(!todoAutomatic_,"this should not happen");
  totNumECVs_=beta_.size();
  ECVs_.resize(beta_.size());
  derECVs_.resize(beta_.size());
  isReady_=true;
  log.printf("  *%4u temperatures for ECV_MULTICANONICAL\n",beta_.size());
}

void ECVmultiCanonical::initECVs_observ(const std::vector<double>& all_obs_cvs,const unsigned ncv,const unsigned index_j)
{
  if(todoAutomatic_)//estimate the steps in beta from observations
  {
    plumed_massert(all_obs_cvs.size()%ncv==0 && index_j<ncv,"initECVs_observ parameters are inconsistent");
    std::vector<double> obs_ene(all_obs_cvs.size()/ncv);//copy only useful observation //TODO we should avoid this...
    for(unsigned t=0; t<obs_ene.size(); t++)
      obs_ene[t]=all_obs_cvs[t*ncv+index_j];
    const double min_beta=beta_[0];
    const double max_beta=beta_[1];
    const unsigned steps_temp=estimate_steps(max_beta-beta0_,min_beta-beta0_,obs_ene,"TEMP");
    log.printf("    (spacing is in beta, not in temperature)\n");
    setBetaSteps(min_beta,max_beta,steps_temp);
  }
  initECVs();
}

void ECVmultiCanonical::initECVs_restart(const std::vector<std::string>& lambdas)
{
  std::size_t pos=lambdas[0].find("_");
  plumed_massert(pos==std::string::npos,"this should not happen, only one CV is used in ECV_MULTICANONICAL");
  if(todoAutomatic_)
    setBetaSteps(beta_[0],beta_[1],lambdas.size());
  std::vector<std::string> myLambdas=getLambdas();
  plumed_massert(myLambdas.size()==lambdas.size(),"RESTART - mismatch in number of ECV_MULTICANONICAL");
  plumed_massert(std::equal(myLambdas.begin(),myLambdas.end(),lambdas.begin()),"RESTART - mismatch in lambda values of ECV_MULTICANONICAL");

  initECVs();
}

}
}
