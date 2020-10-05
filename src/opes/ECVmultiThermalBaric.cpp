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

//+PLUMEDOC EXPANSION_CV ECV_MULTITHERMAL_MULTIBARIC
/*
Expand a canonical simulation to sample multiple temperatures.
If instead of fixed volume NVT you are running with fixed pressure NPT, you must use \ref ECV_MULTITHERMAL_MULTIBARIC and add the volume contribution.
The \ref ENERGY of the system should be used as ARG.

\par Examples

mc: ECV_MULTITHERMAL_MULTIBARIC ARG=ene,vol TEMP=300 MAX_TEMP=500 PRESSURE=0.06022140857

*/
//+ENDPLUMEDOC

class ECVmultiThermalBaric :
  public ExpansionCVs
{
private:
  bool todoAutomatic_beta_;
  bool todoAutomatic_pres_;
  double beta0_;
  double pres0_;
  std::vector<double> beta_;
  std::vector<double> pres_;
  std::vector<double> ECVs_beta_;
  std::vector<double> ECVs_pres_;
  std::vector<double> derECVs_beta_;
  std::vector<double> derECVs_pres_;
  void setBetaSteps(double,double,unsigned);
  void setPresSteps(double,double,unsigned);
  void initECVs();

  std::vector<double> cut_corner_;

public:
  explicit ECVmultiThermalBaric(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  void calculateECVs(const double *) override;
  const double * getPntrToECVs(unsigned) override;
  const double * getPntrToDerECVs(unsigned) override;
  std::vector< std::vector<unsigned> > getIndex_k() const override;
  std::vector<std::string> getLambdas() const override;
  void initECVs_observ(const std::vector<double>&,const unsigned,const unsigned) override;
  void initECVs_restart(const std::vector<std::string>&) override;
};

PLUMED_REGISTER_ACTION(ECVmultiThermalBaric,"ECV_MULTITHERMAL_MULTIBARIC")

void ECVmultiThermalBaric::registerKeywords(Keywords& keys) {
  ExpansionCVs::registerKeywords(keys);
  keys.remove("ARG");
  keys.add("compulsory","ARG","provide the labels of the potential energy and of the volume of the system, you can calculate them with \\ref ENERGY and \\ref VOLUME respectively");
//temperature
  keys.add("compulsory","TEMP","-1","temperature. If not specified tries to get it from MD engine");
  keys.add("optional","MIN_TEMP","the minimum of the temperature range");
  keys.add("optional","MAX_TEMP","the maximum of the temperature range");
  keys.add("optional","STEPS_TEMP","the number of steps in temperature");
  keys.add("optional","SET_ALL_TEMPS","manually set all the temperatures");
//pressure
  keys.add("compulsory","PRESSURE","pressure. Use the proper units");
  keys.add("optional","MIN_PRESSURE","the minimum of the pressure range");
  keys.add("optional","MAX_PRESSURE","the maximum of the pressure range");
  keys.add("optional","STEPS_PRESSURE","the number of steps in pressure");
  keys.add("optional","SET_ALL_PRESSURES","manually set all the pressures");
//other
  keys.add("optional","CUT_CORNER","avoid region of high temperature and low pressure. This excludes all points below a line in the temperature-pressure plane, defined by two points: T_low,P_low,T_high,P_high");
//  keys.add("optional","BORDER_WEIGHT","set it greater than 1 to obtain better sampling of the max and min thermodynamics conditions");
}

ECVmultiThermalBaric::ECVmultiThermalBaric(const ActionOptions&ao):
  Action(ao),
  ExpansionCVs(ao)
{
  plumed_massert(getNumberOfArguments()==2,"ENERGY and VOLUME should be given as ARG");

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

//parse pressures
  parse("PRESSURE",pres0_);
  double min_pres=std::numeric_limits<double>::quiet_NaN();//-1 might be a meaningful pressure
  double max_pres=std::numeric_limits<double>::quiet_NaN();
  parse("MIN_PRESSURE",min_pres);
  parse("MAX_PRESSURE",max_pres);
  unsigned steps_pres=0;
  parse("STEPS_PRESSURE",steps_pres);
  parseVector("SET_ALL_PRESSURES",pres_);

//other
  std::vector<double> cut_corner_;
  parseVector("CUT_CORNER",cut_corner_);

  checkRead();
  log.printf("  running at TEMP=%g and PRESSURE=%g\n",temp,pres0_);

//set the intermediate temperatures
  if(temps.size()>0)
  {
    todoAutomatic_beta_=false;
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
    if(min_temp!=-1 && max_temp!=-1)
      log.printf(" +++ WARNING +++ running at single temperature. If you only need a multibaric simulation it is more efficient to set it up with ECV_LINEAR\n");
    plumed_massert(max_temp>=min_temp,"MAX_TEMP should be bigger than MIN_TEMP");
    const double min_beta=1./(Kb*max_temp);
    const double max_beta=1./(Kb*min_temp);
    if(steps_temp>0)
      setBetaSteps(min_beta,max_beta,steps_temp);
    else
    {
      todoAutomatic_beta_=true;
      beta_.resize(2);
      beta_[0]=min_beta;
      beta_[1]=max_beta;
    }
  }
  const double tol=1e-3; //if temp is taken from MD engine it might be numerically slightly different
  if(temp<(1-tol)*min_temp || temp>(1+tol)*max_temp)
    log.printf(" +++ WARNING +++ running at TEMP=%g which is outside the chosen temperature range\n",temp);

//set the intermediate pressures
  if(pres_.size()>0)
  {
    todoAutomatic_pres_=false;
    plumed_massert(steps_pres==0,"cannot set both STEPS_PRESSURE and SET_ALL_PRESSURES");
    plumed_massert(std::isnan(min_pres) && std::isnan(max_pres),"cannot set both SET_ALL_PRESSURES and MIN/MAX_PRESSURE");
    plumed_massert(pres_.size()>=2,"set at least 2 pressures");
    for(unsigned kk=0; kk<pres_.size()-1; kk++)
      plumed_massert(pres_[kk]<pres_[kk+1],"SET_ALL_PRESSURES must be properly ordered");
    min_pres=pres_[0];
    max_pres=pres_[pres_.size()-1];
  }
  else
  { //get MIN_PRESSURE and MAX_PRESSURE
    if(std::isnan(min_pres))
    {
      min_pres=pres0_;
      log.printf("  no MIN_PRESSURE provided, using MIN_PRESSURE=PRESSURE\n");
    }
    if(std::isnan(max_pres))
    {
      max_pres=pres0_;
      log.printf("  no MAX_PRESSURE provided, using MAX_PRESSURE=PRESSURE\n");
    }
    plumed_massert(max_pres>=min_pres,"MAX_PRESSURE should be bigger than MIN_PRESSURE");
    if(steps_pres>0)
      setPresSteps(min_pres,max_pres,steps_pres);
    else
    {
      todoAutomatic_pres_=true;
      pres_.resize(2);
      pres_[0]=min_pres;
      pres_[1]=max_pres;
    }
  }
  if(pres0_<min_pres || pres0_>max_pres)
    log.printf(" +++ WARNING +++ running at PRESSURE=%g which is outside the chosen pressure range\n",pres0_);

//other
  if(cut_corner_.size()>0)
  {
    std::string cc_usage("CUT_CORNER=low_temp,low_pres,high_temp,high_pres");
    plumed_massert(cut_corner_.size()==4,"expected 4 values for "+cc_usage);
    plumed_massert(cut_corner_[0]<cut_corner_[2],"low_temp="+std::to_string(cut_corner_[0])+" should be smaller than high_temp="+std::to_string(cut_corner_[2])+", "+cc_usage);
    plumed_massert(cut_corner_[0]>=min_temp && cut_corner_[0]<=max_temp,"low_temp="+std::to_string(cut_corner_[0])+" is out of temperature range. "+cc_usage);
    plumed_massert(cut_corner_[2]>=min_temp && cut_corner_[2]<=max_temp,"high_temp="+std::to_string(cut_corner_[2])+" is out of temperature range. "+cc_usage);
    plumed_massert(cut_corner_[1]<cut_corner_[3],"low_pres="+std::to_string(cut_corner_[1])+" should be smaller than high_pres="+std::to_string(cut_corner_[3])+", "+cc_usage);
    plumed_massert(cut_corner_[1]>=min_pres && cut_corner_[1]<=max_pres,"low_pres="+std::to_string(cut_corner_[1])+" is out of pressure range. "+cc_usage);
    plumed_massert(cut_corner_[3]>=min_pres && cut_corner_[3]<=max_pres,"high_pres="+std::to_string(cut_corner_[3])+" is out of pressure range. "+cc_usage);
  }
}

void ECVmultiThermalBaric::calculateECVs(const double * ene_vol) {
  for(unsigned k=0; k<beta_.size(); k++)
  {
    const double diff_k=(beta_[k]-beta0_);
    ECVs_beta_[k]=diff_k*ene_vol[0];
    derECVs_beta_[k]=diff_k;
  }
  for(unsigned kk=0; kk<pres_.size(); kk++)
  {
    const double diff_kk=(pres_[kk]-pres0_);
    ECVs_pres_[kk]=diff_kk*ene_vol[1];
    derECVs_pres_[kk]=diff_kk;
  }
}

const double * ECVmultiThermalBaric::getPntrToECVs(unsigned j)
{
  plumed_massert(isReady_,"cannot access ECVs before initialization");
  plumed_massert(j==0 || j==1,"ECV_MULTITHERMAL_MULTIBARIC has only two CVs, the ENERGY and the VOLUME");
  if(j==0)
    return &ECVs_beta_[0];
  else //if (j==1)
    return &ECVs_pres_[0];
}

const double * ECVmultiThermalBaric::getPntrToDerECVs(unsigned j)
{
  plumed_massert(isReady_,"cannot access ECVs before initialization");
  plumed_massert(j==0 || j==1,"ECV_MULTITHERMAL_MULTIBARIC has only two CVs, the ENERGY and the VOLUME");
  if(j==0)
    return &derECVs_beta_[0];
  else //if (j==1)
    return &derECVs_pres_[0];
}

std::vector< std::vector<unsigned> > ECVmultiThermalBaric::getIndex_k() const
{
  plumed_massert(!todoAutomatic_beta_ && !todoAutomatic_pres_,"cannot access index_k before initializing lambdas");
  std::vector< std::vector<unsigned> > index_k;
  if(cut_corner_.size()==0)
  {
    index_k.resize(beta_.size()*pres_.size(),std::vector<unsigned>(2));
    for(unsigned k=0; k<beta_.size(); k++)
      for(unsigned kk=0; kk<pres_.size(); kk++)
      {
        const unsigned i=k*pres_.size()+kk;
        index_k[i][0]=k;
        index_k[i][1]=kk;
      }
  }
  else
  {
    const double Kb=plumed.getAtoms().getKBoltzmann();
    const double low_temp=cut_corner_[0];
    const double low_pres=cut_corner_[1];
    const double high_temp=cut_corner_[2];
    const double high_pres=cut_corner_[3];
    const double coeff=(high_pres-low_pres)/(high_temp-low_temp);
    for(unsigned k=0; k<beta_.size(); k++)
    {
      const double line_k=(1./(Kb*beta_[k])-low_temp)*coeff+low_pres;
      for(unsigned kk=0; kk<pres_.size(); kk++)
      {
        if(pres_[kk]>=line_k) //important to be inclusive, thus >=, not just >
          index_k.emplace_back(std::vector<unsigned>{k,kk});
      }
    }
    plumed_massert(beta_.size()*pres_.size()>=index_k.size(),"this should not happen, something is wrong with CUT_CORNER");
  }
  return index_k;
}

std::vector<std::string> ECVmultiThermalBaric::getLambdas() const
{
  plumed_massert(!todoAutomatic_beta_ && !todoAutomatic_pres_,"cannot access lambdas before initializing them");
  std::vector< std::vector<unsigned> > index_k=getIndex_k();//slow, but not in a critical place...
  std::vector<std::string> lambdas(index_k.size());
  const double Kb=plumed.getAtoms().getKBoltzmann();
  for(unsigned i=0; i<lambdas.size(); i++)
  {
    std::ostringstream subs;
    subs<<1./(Kb*beta_[index_k[i][0]])<<"_"<<pres_[index_k[i][1]];
    lambdas[i]=subs.str();
  }
  return lambdas;
}

void ECVmultiThermalBaric::setBetaSteps(double min_beta,double max_beta,unsigned steps_beta)
{
  plumed_massert(beta_.size()==0 || beta_.size()==2,"you should not set the beta steps twice...");
  plumed_massert(min_beta==max_beta && steps_beta>1,"cannot have multiple STEPS_TEMP if MIN_TEMP==MAX_TEMP");
  beta_.resize(steps_beta);
  if(steps_beta==1)
  {
    beta_[0]=(min_beta+max_beta)/2.;
    log.printf(" +++ WARNING +++ using one single temperature as target, corresponding to beta = %g\n",beta_[0]);
  }
  else
    for(unsigned k=0; k<beta_.size(); k++)//betas are stored in reversed order, so temps are in correct order
      beta_[k]=max_beta-k*(max_beta-min_beta)/(steps_beta-1);
  todoAutomatic_beta_=false;
}

void ECVmultiThermalBaric::setPresSteps(double min_pres,double max_pres,unsigned steps_pres)
{
  plumed_massert(pres_.size()==0 || pres_.size()==2,"you should not set the pres steps twice...");
  plumed_massert(min_pres==max_pres && steps_pres>1,"cannot have multiple STEPS_PRESSURE if MIN_PRESSURE==MAX_PRESSURE");
  pres_.resize(steps_pres);
  if(steps_pres==1)
  {
    pres_[0]=(min_pres+max_pres)/2.;
    log.printf(" +++ WARNING +++ using one single pressure as target = %g\n",pres_[0]);
  }
  else
    for(unsigned kk=0; kk<pres_.size(); kk++)//betas are stored in reversed order, so temps are in correct order
      pres_[kk]=min_pres+kk*(max_pres-min_pres)/(steps_pres-1);
  todoAutomatic_pres_=false;
}

void ECVmultiThermalBaric::initECVs()
{
  plumed_massert(!isReady_,"initialization should not be called twice");
  plumed_massert(!todoAutomatic_beta_ && !todoAutomatic_pres_,"this should not happen");
  totNumECVs_=getIndex_k().size();//get it again, to be sure
  ECVs_beta_.resize(beta_.size());
  ECVs_pres_.resize(pres_.size());
  derECVs_beta_.resize(beta_.size());
  derECVs_pres_.resize(pres_.size());
  isReady_=true;
  log.printf("  *%4u temperatures for ECV_MULTITHERMAL_MULTIBARIC\n",beta_.size());
  log.printf("  *%4u pressures for ECV_MULTITHERMAL_MULTIBARIC\n",beta_.size());
  if(cut_corner_.size()!=0)
    log.printf("    -- CUT_CORNER: %u temp-pres points were excluded, thus total is %u\n",beta_.size()*pres_.size()-totNumECVs_,totNumECVs_);
}

void ECVmultiThermalBaric::initECVs_observ(const std::vector<double>& all_obs_cvs,const unsigned ncv,const unsigned index_j)
{
  if(todoAutomatic_beta_)//estimate the steps in beta from observations
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
  if(todoAutomatic_pres_)//estimate the steps in pres from observations
  {
    plumed_massert(all_obs_cvs.size()%ncv==0 && index_j+1<ncv,"initECVs_observ parameters are inconsistent");
    std::vector<double> obs_vol(all_obs_cvs.size()/ncv);//copy only useful observation //TODO we should avoid this...
    for(unsigned t=0; t<obs_vol.size(); t++)
      obs_vol[t]=all_obs_cvs[t*ncv+index_j+1];
    const double min_pres=pres_[0];
    const double max_pres=pres_[1];
    const unsigned steps_pres=estimate_steps(beta0_*(min_pres-pres0_),beta0_*(max_pres-pres0_),obs_vol,"PRESSURE");
    log.printf("    (spacing is in beta0 units)\n");
    setPresSteps(min_pres,max_pres,steps_pres);
  }
  initECVs();
}

void ECVmultiThermalBaric::initECVs_restart(const std::vector<std::string>& lambdas)
{
  std::size_t pos=lambdas[0].find("_");
  plumed_massert(pos!=std::string::npos,"this should not happen, two CVs are used in ECV_MULTITHERMAL_MULTIBARIC, not less");
  lambdas[0].find("_",pos+1);
  plumed_massert(pos==std::string::npos,"this should not happen, two CVs are used in ECV_MULTITHERMAL_MULTIBARIC, not more");

  auto getLambdaName=[](std::string name,const unsigned index_j) //slightly different from the one in OPESexpanded.cpp
  {
    std::size_t pos=0;
    for(unsigned j=0; j<index_j; j++)
      pos=name.find("_",pos+1);
    std::size_t pos_end=name.find("_",pos+1);
    pos+=index_j;
    return name.substr(pos,pos_end-pos);
  };
  if(todoAutomatic_pres_)
  {
    unsigned steps_pres=0;
    std::string pres_min=getLambdaName(lambdas[0],1);
    for(unsigned i=1; i<lambdas.size(); i++)//pres is second, thus increas by 1
    {
      if(getLambdaName(lambdas[i],1)==pres_min)
        break;
      steps_pres++;
    }
    setPresSteps(pres_[0],pres_[1],steps_pres);
  }
  if(todoAutomatic_beta_)
  {
    unsigned steps_temp=1;
    std::string pres_max=getLambdaName(lambdas[pres_.size()-1],1);
    for(unsigned i=pres_.size(); i<lambdas.size(); i++)
    { //if CUT_CORNER below, the max pressures are all present, for each temp
      if(getLambdaName(lambdas[i],1)==pres_max)
        steps_temp++;
    }
    setBetaSteps(beta_[0],beta_[1],steps_temp);
  }
  std::vector<std::string> myLambdas=getLambdas();
  plumed_massert(myLambdas.size()==lambdas.size(),"RESTART - mismatch in number of ECV_MULTITHERMAL_MULTIBARIC");
  plumed_massert(std::equal(myLambdas.begin(),myLambdas.end(),lambdas.begin()),"RESTART - mismatch in lambda values of ECV_MULTITHERMAL_MULTIBARIC");

  initECVs();
}

}
}
