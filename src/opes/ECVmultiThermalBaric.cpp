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
Expand a simulation to sample multiple temperatures and pressures.

The potential \ref ENERGY, \f$E\f$, and the \ref VOLUME, \f$V\f$, of the system should be used as ARG.
\f[
  \Delta u_{\beta',p'}=(\beta-\beta') E + (\beta p -\beta' p') V\, .
\f]

If instead you wish to sample multiple temperatures and a single pressure, you should use \ref ECV_MULTITHERMAL with as ARG the internal energy \f$U=E+pV\f$.

The STEPS_TEMP and STEPS_PRESSURE are automatically guessed from the initial unbiased steps (see OBSERVATION_STEPS in \ref OPES_EXPANDED), unless explicitly set.
The temperatures are chosen with geometric distribution (uniform in beta), while the pressures are uniformely spaced.
For more detailed control you can use SET_ALL_TEMPS and SET_ALL_PRESSURES.
The temperatures and pressures are then combined in a 2D grid.

You can use CUT_CORNER to avoid a high-temperature/low-pressure region.
This can be useful e.g. to increase the temperature for greater ergodicity, while avoiding water to vaporize, as in Ref.\cite Invernizzi2020unified.

You can reweight the resulting simulation at any temperature and pressure in chosen target, using e.g. \ref REWEIGHT_TEMP_PRESS.
A similar target distribution can be sampled using \ref TD_MULTITHERMAL_MULTIBARIC.

\par Examples

\plumedfile
ene: ENERGY
vol: VOLUME
ecv: ECV_MULTITHERMAL_MULTIBARIC ...
  ARG=ene,vol
  TEMP=500
  MIN_TEMP=270
  MAX_TEMP=800
  PRESSURE=0.06022140857*2000 #2 kbar
  MIN_PRESSURE=0.06022140857  #1 bar
  MAX_PRESSURE=0.06022140857*4000 #4 kbar
  CUT_CORNER=500,0.06022140857,800,0.06022140857*1000
...
opes: OPES_EXPANDED ARG=ecv.* FILE=DeltaF.data PACE=500 WALKERS_MPI
\endplumedfile

Notice that \f$p=0.06022140857\f$ corresponds to 1 bar only when using the default PLUMED units.
If you modify them via the \ref UNITS command, then the pressure has to be rescaled accordingly.

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
  void initECVs();

//CUT_CORNER stuff
  double coeff_;
  double low_pres_;
  double low_temp_Kb_;

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

void ECVmultiThermalBaric::registerKeywords(Keywords& keys)
{
  ExpansionCVs::registerKeywords(keys);
  keys.remove("ARG");
  keys.add("compulsory","ARG","the labels of the potential energy and of the volume of the system. You can calculate them with \\ref ENERGY and \\ref VOLUME respectively");
//temperature
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
//  keys.add("optional","SET_GRID_TEMPS_PRESSURES","manually set the whole temperature-pressure grid. Comma separated points, with internal underscore, e.g.: temp1_pres1,temp1_pres2,...");
  keys.add("optional","CUT_CORNER","avoid region of high temperature and low pressure. Exclude all points below a line in the temperature-pressure plane, defined by two points: \\f$T_{\\text{low}},P_{\\text{low}},T_{\\text{high}},P_{\\text{high}}\\f$");
}

ECVmultiThermalBaric::ECVmultiThermalBaric(const ActionOptions&ao)
  : Action(ao)
  , ExpansionCVs(ao)
  , todoAutomatic_beta_(false)
  , todoAutomatic_pres_(false)
  , beta0_(1./kbt_)
{
  plumed_massert(getNumberOfArguments()==2,"ENERGY and VOLUME should be given as ARG");

//set temp0
  const double Kb=plumed.getAtoms().getKBoltzmann();
  const double temp0=kbt_/Kb;

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
  double min_pres=std::numeric_limits<double>::quiet_NaN(); //-1 might be a meaningful pressure
  double max_pres=std::numeric_limits<double>::quiet_NaN();
  parse("MIN_PRESSURE",min_pres);
  parse("MAX_PRESSURE",max_pres);
  unsigned steps_pres=0;
  parse("STEPS_PRESSURE",steps_pres);
  parseVector("SET_ALL_PRESSURES",pres_);
//other
  std::vector<double> cut_corner;
  parseVector("CUT_CORNER",cut_corner);

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
      todoAutomatic_beta_=true;
  }
  const double tol=1e-3; //if temp is taken from MD engine it might be numerically slightly different
  if(temp0<(1-tol)*min_temp || temp0>(1+tol)*max_temp)
    log.printf(" +++ WARNING +++ running at TEMP=%g which is outside the chosen temperature range\n",temp0);

//set the intermediate pressures
  if(pres_.size()>0)
  {
    plumed_massert(steps_pres==0,"cannot set both STEPS_PRESSURE and SET_ALL_PRESSURES");
    plumed_massert(std::isnan(min_pres) && std::isnan(max_pres),"cannot set both SET_ALL_PRESSURES and MIN/MAX_PRESSURE");
    plumed_massert(pres_.size()>=2,"set at least 2 pressures");
    for(unsigned kk=0; kk<pres_.size()-1; kk++)
      plumed_massert(pres_[kk]<=pres_[kk+1],"SET_ALL_PRESSURES must be properly ordered");
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
    pres_.resize(2);
    pres_[0]=min_pres;
    pres_[1]=max_pres;
    if(min_pres==max_pres && steps_pres==0)
      steps_pres=1;
    if(steps_pres>0)
      setSteps(pres_,steps_pres,"PRESSURE");
    else
      todoAutomatic_pres_=true;
  }
  if(pres0_<min_pres || pres0_>max_pres)
    log.printf(" +++ WARNING +++ running at PRESSURE=%g which is outside the chosen pressure range\n",pres0_);

//set CUT_CORNER
  std::string cc_usage("CUT_CORNER=low_temp,low_pres,high_temp,high_pres");
  if(cut_corner.size()>0)
    plumed_massert(cut_corner.size()==4,"expected 4 values for "+cc_usage);
  if(cut_corner.size()==4) //this keeps cppcheck happy
  {
    const double low_temp=cut_corner[0];
    const double low_pres=cut_corner[1];
    const double high_temp=cut_corner[2];
    const double high_pres=cut_corner[3];
    plumed_massert(low_temp<high_temp,"low_temp="+std::to_string(low_temp)+" should be smaller than high_temp="+std::to_string(high_temp)+", "+cc_usage);
    plumed_massert(low_temp>=min_temp && low_temp<=max_temp,"low_temp="+std::to_string(low_temp)+" is out of temperature range. "+cc_usage);
    plumed_massert(high_temp>=min_temp && high_temp<=max_temp,"high_temp="+std::to_string(high_temp)+" is out of temperature range. "+cc_usage);
    plumed_massert(low_pres<high_pres,"low_pres="+std::to_string(low_pres)+" should be smaller than high_pres="+std::to_string(high_pres)+", "+cc_usage);
    plumed_massert(low_pres>=min_pres && low_pres<=max_pres,"low_pres="+std::to_string(low_pres)+" is out of pressure range. "+cc_usage);
    plumed_massert(high_pres>=min_pres && high_pres<=max_pres,"high_pres="+std::to_string(high_pres)+" is out of pressure range. "+cc_usage);
    low_temp_Kb_=low_temp*Kb;
    coeff_=(high_pres-low_pres)/(high_temp-low_temp)/Kb;
    plumed_massert(coeff_!=0,"this should not be possible");
    const double small_value=(high_temp-low_pres)/1e4;
    low_pres_=low_pres-small_value; //make sure max_pres is included
    plumed_massert(max_pres>=coeff_*(1./beta_[beta_.size()-1]-low_temp_Kb_)+low_pres_,"please chose a high_pres slightly smaller than MAX_PRESSURE in "+cc_usage);
  }
  else
  {
    coeff_=0;
    low_pres_=0;
    low_temp_Kb_=0;
  }

//print some info
  log.printf("  running at TEMP=%g and PRESSURE=%g\n",temp0,pres0_);
  log.printf("  targeting a temperature range from MIN_TEMP=%g to MAX_TEMP=%g\n",min_temp,max_temp);
  if(min_temp==max_temp)
    log.printf(" +++ WARNING +++ if you only need a multibaric simulation it is more efficient to set it up with ECV_LINEAR\n");
  log.printf("   and a pressure range from MIN_PRESSURE=%g to MAX_PRESSURE=%g\n",min_pres,max_pres);
  if(min_pres==max_pres)
    log.printf(" +++ WARNING +++ if you only need a multithermal simulation it is more efficient to set it up with ECV_MULTITHERMAL\n");
  if(coeff_!=0)
    log.printf(" -- CUT_CORNER: ignoring some high temperature and low pressure values\n");
}

void ECVmultiThermalBaric::calculateECVs(const double * ene_vol)
{
  unsigned i=0;
  for(unsigned k=0; k<beta_.size(); k++)
  {
    const double diff_k=(beta_[k]-beta0_);
    ECVs_beta_[k]=diff_k*ene_vol[0];
    derECVs_beta_[k]=diff_k;
    const double line_k=coeff_*(1./beta_[k]-low_temp_Kb_)+low_pres_;
    for(unsigned kk=0; kk<pres_.size(); kk++)
    {
      if(coeff_==0 || pres_[kk]>=line_k)
      {
        const double diff_i=(beta_[k]*pres_[kk]-beta0_*pres0_); //this is not great, each beta-pres combination must be stored separately
        ECVs_pres_[i]=diff_i*ene_vol[1];
        derECVs_pres_[i]=diff_i;
        i++;
      }
    }
  }
}

const double * ECVmultiThermalBaric::getPntrToECVs(unsigned j)
{
  plumed_massert(isReady_,"cannot access ECVs before initialization");
  plumed_massert(j==0 || j==1,getName()+" has only two CVs, the ENERGY and the VOLUME");
  if(j==0)
    return &ECVs_beta_[0];
  else //if (j==1)
    return &ECVs_pres_[0];
}

const double * ECVmultiThermalBaric::getPntrToDerECVs(unsigned j)
{
  plumed_massert(isReady_,"cannot access ECVs before initialization");
  plumed_massert(j==0 || j==1,getName()+" has only two CVs, the ENERGY and the VOLUME");
  if(j==0)
    return &derECVs_beta_[0];
  else //if (j==1)
    return &derECVs_pres_[0];
}

std::vector< std::vector<unsigned> > ECVmultiThermalBaric::getIndex_k() const
{
  plumed_massert(isReady_ && totNumECVs_>0,"cannot access getIndex_k() of ECV before initialization");
  std::vector< std::vector<unsigned> > index_k;
  unsigned i=0;
  for(unsigned k=0; k<beta_.size(); k++)
  {
    const double line_k=coeff_*(1./beta_[k]-low_temp_Kb_)+low_pres_;
    for(unsigned kk=0; kk<pres_.size(); kk++)
    {
      if(coeff_==0 || pres_[kk]>=line_k) //important to be inclusive, thus >=, not just >
      {
        index_k.emplace_back(std::vector<unsigned> {k,i});
        i++;
      }
    }
  }
  plumed_massert(totNumECVs_==index_k.size(),"this should not happen, is something wrong with CUT_CORNER ?");
  return index_k;
}

std::vector<std::string> ECVmultiThermalBaric::getLambdas() const
{
  plumed_massert(!todoAutomatic_beta_ && !todoAutomatic_pres_,"cannot access lambdas before initializing them");
  std::vector<std::string> lambdas;
  const double Kb=plumed.getAtoms().getKBoltzmann();
  for(unsigned k=0; k<beta_.size(); k++)
  {
    const double line_k=coeff_*(1./beta_[k]-low_temp_Kb_)+low_pres_;
    for(unsigned kk=0; kk<pres_.size(); kk++)
    {
      if(coeff_==0 || pres_[kk]>=line_k)
      {
        std::ostringstream subs;
        subs<<1./(Kb*beta_[k])<<"_"<<pres_[kk];
        lambdas.emplace_back(subs.str());
      }
    }
  }
  return lambdas;
}

void ECVmultiThermalBaric::initECVs()
{
  plumed_massert(!isReady_,"initialization should not be called twice");
  plumed_massert(!todoAutomatic_beta_ && !todoAutomatic_pres_,"this should not happen");
  totNumECVs_=getLambdas().size(); //slow, but runs only once
  plumed_massert(beta_.size()*pres_.size()>=totNumECVs_,"this should not happen, is something wrong with CUT_CORNER ?");
  ECVs_beta_.resize(beta_.size());
  derECVs_beta_.resize(beta_.size());
  ECVs_pres_.resize(totNumECVs_); //pres is mixed with temp (beta*p*V), thus we need to store all possible
  derECVs_pres_.resize(totNumECVs_);
  isReady_=true;
  log.printf("  *%4lu temperatures for %s\n",beta_.size(),getName().c_str());
  log.printf("  *%4lu pressures for %s\n",pres_.size(),getName().c_str());
  if(coeff_!=0)
    log.printf("    -- CUT_CORNER: %lu temp-pres points were excluded, thus total is %u\n",beta_.size()*pres_.size()-totNumECVs_,totNumECVs_);
}

void ECVmultiThermalBaric::initECVs_observ(const std::vector<double>& all_obs_cvs,const unsigned ncv,const unsigned index_j)
{
  if(todoAutomatic_beta_) //estimate the steps in beta from observations
  {
    plumed_massert(all_obs_cvs.size()%ncv==0 && index_j<ncv,"initECVs_observ parameters are inconsistent");
    std::vector<double> obs_ene(all_obs_cvs.size()/ncv); //copy only useful observations
    for(unsigned t=0; t<obs_ene.size(); t++)
      obs_ene[t]=all_obs_cvs[t*ncv+index_j]+pres0_*all_obs_cvs[t*ncv+index_j+1]; //U=E+pV
    const unsigned steps_temp=estimateSteps(beta_[0]-beta0_,beta_[1]-beta0_,obs_ene,"TEMP");
    log.printf("    (spacing is on beta, not on temperature)\n");
    setSteps(beta_,steps_temp,"TEMP");
    todoAutomatic_beta_=false;
  }
  if(todoAutomatic_pres_) //estimate the steps in pres from observations
  {
    plumed_massert(all_obs_cvs.size()%ncv==0 && index_j+1<ncv,"initECVs_observ parameters are inconsistent");
    std::vector<double> obs_vol(all_obs_cvs.size()/ncv); //copy only useful observations
    for(unsigned t=0; t<obs_vol.size(); t++)
      obs_vol[t]=all_obs_cvs[t*ncv+index_j+1];
    const unsigned steps_pres=estimateSteps(beta0_*(pres_[0]-pres0_),beta0_*(pres_[1]-pres0_),obs_vol,"PRESSURE");
    log.printf("    (spacing is in beta0 units)\n");
    setSteps(pres_,steps_pres,"PRESSURE");
    todoAutomatic_pres_=false;
  }
  initECVs();
  calculateECVs(&all_obs_cvs[index_j]);
}

void ECVmultiThermalBaric::initECVs_restart(const std::vector<std::string>& lambdas)
{
  std::size_t pos=lambdas[0].find("_");
  plumed_massert(pos!=std::string::npos,"this should not happen, two CVs are used in "+getName()+", not less");
  pos=lambdas[0].find("_",pos+1);
  plumed_massert(pos==std::string::npos,"this should not happen, two CVs are used in "+getName()+", not more");

  auto getPres=[&lambdas](const unsigned i) {return lambdas[i].substr(lambdas[i].find("_")+1);};
  if(todoAutomatic_pres_)
  {
    unsigned steps_pres=1;
    std::string pres_min=getPres(0);
    for(unsigned i=1; i<lambdas.size(); i++) //pres is second, thus increas by 1
    {
      if(getPres(i)==pres_min)
        break;
      steps_pres++;
    }
    setSteps(pres_,steps_pres,"PRESSURE");
    todoAutomatic_pres_=false;
  }
  if(todoAutomatic_beta_)
  {
    unsigned steps_temp=1;
    std::string pres_max=getPres(pres_.size()-1);
    for(unsigned i=pres_.size(); i<lambdas.size(); i++)
    { //even if CUT_CORNER, the max pressures are all present, for each temp
      if(getPres(i)==pres_max)
        steps_temp++;
    }
    setSteps(beta_,steps_temp,"TEMP");
    todoAutomatic_beta_=false;
  }
  std::vector<std::string> myLambdas=getLambdas();
  plumed_assert(myLambdas.size()==lambdas.size())<<"RESTART - mismatch in number of "<<getName()<<".\nFrom "<<lambdas.size()<<" labels "<<beta_.size()<<" temperatures and "<<pres_.size()<<" pressures were found, for a total of "<<myLambdas.size()<<" estimated steps.\nCheck if the CUT_CORNER option is consistent\n";
  plumed_massert(std::equal(myLambdas.begin(),myLambdas.end(),lambdas.begin()),"RESTART - mismatch in lambda values of "+getName());

  initECVs();
}

}
}
