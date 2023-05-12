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
#include "core/PlumedMain.h"
#include "core/Atoms.h"

namespace PLMD {
namespace opes {

//+PLUMEDOC OPES_EXPANSION_CV ECV_MULTITHERMAL_MULTIBARIC
/*
Expand a simulation to sample multiple temperatures and pressures.

The potential \ref ENERGY, \f$E\f$, and the \ref VOLUME, \f$V\f$, of the system should be used as ARG.
\f[
  \Delta u_{\beta',p'}=(\beta'-\beta) E + (\beta' p' -\beta p) V\, ,
\f]
where \f$\beta', p'\f$ are the temperatures and pressures to be sampled, while \f$\beta, p\f$ is the temperature and pressure at which the simulation is conducted.

If instead you wish to sample multiple temperatures and a single pressure, you should use \ref ECV_MULTITHERMAL with as ARG the internal energy \f$U=E+pV\f$.

The TEMP_STEPS and PRESSURE_STEPS are automatically guessed from the initial unbiased steps (see OBSERVATION_STEPS in \ref OPES_EXPANDED), unless explicitly set.
The algorithm for this guess is described in \cite Invernizzi2020unified should provide a rough estimate useful for most applications.
The pressures are uniformely spaced, while the temperatures steps are geometrically spaced.
Use instead the keyword NO_GEOM_SPACING for a linear spacing in inverse temperature (beta).
For more detailed control you can use instead TEMP_SET_ALL and/or PRESSURE_SET_ALL to explicitly set all of them.
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
  TEMP_MIN=270
  TEMP_MAX=800
  PRESSURE=0.06022140857*2000 #2 kbar
  PRESSURE_MIN=0.06022140857  #1 bar
  PRESSURE_MAX=0.06022140857*4000 #4 kbar
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
  bool geom_spacing_;
  double pres0_;
  std::vector<double> pres_;
  std::vector<double> ECVs_beta_;
  std::vector<double> ECVs_pres_;
  std::vector<double> derECVs_beta_; //(beta_k-beta0) or (temp0/temp_k-1)/kbt
  std::vector<double> derECVs_pres_; //(beta_k*pres_kk-beta0*pres0) or (temp0/temp_k*pres_kk-pres0)/kbt
  void initECVs();

//CUT_CORNER stuff
  double coeff_;
  double pres_low_;
  double kB_temp_low_;
//SET_ALL_TEMP_PRESSURE stuff
  std::vector<std::string> custom_lambdas_;

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
  keys.add("optional","TEMP_MIN","the minimum of the temperature range");
  keys.add("optional","TEMP_MAX","the maximum of the temperature range");
  keys.add("optional","TEMP_STEPS","the number of steps in temperature");
  keys.add("optional","TEMP_SET_ALL","manually set all the temperatures");
  keys.addFlag("NO_GEOM_SPACING",false,"do not use geometrical spacing in temperature, but instead linear spacing in inverse temperature");
//pressure
  keys.add("compulsory","PRESSURE","pressure. Use the proper units");
  keys.add("optional","PRESSURE_MIN","the minimum of the pressure range");
  keys.add("optional","PRESSURE_MAX","the maximum of the pressure range");
  keys.add("optional","PRESSURE_STEPS","the number of steps in pressure");
  keys.add("optional","PRESSURE_SET_ALL","manually set all the pressures");
//other
  keys.add("optional","SET_ALL_TEMP_PRESSURE","manually set all the target temperature_pressure pairs. An underscore separates temperature and pressure, while different points are comma-separated, e.g.: temp1_pres1,temp1_pres2,...");
  keys.add("optional","CUT_CORNER","avoid region of high temperature and low pressure. Exclude all points below a line in the temperature-pressure plane, defined by two points: \\f$T_{\\text{low}},P_{\\text{low}},T_{\\text{high}},P_{\\text{high}}\\f$");
}

ECVmultiThermalBaric::ECVmultiThermalBaric(const ActionOptions&ao)
  : Action(ao)
  , ExpansionCVs(ao)
  , todoAutomatic_beta_(false)
  , todoAutomatic_pres_(false)
  , coeff_(0)
  , pres_low_(0)
  , kB_temp_low_(0)
{
  plumed_massert(getNumberOfArguments()==2,"ENERGY and VOLUME should be given as ARG");

//set temp0
  const double kB=plumed.getAtoms().getKBoltzmann();
  const double temp0=kbt_/kB;

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
//parse pressures
  parse("PRESSURE",pres0_);
  double pres_min=std::numeric_limits<double>::quiet_NaN(); //-1 might be a meaningful pressure
  double pres_max=std::numeric_limits<double>::quiet_NaN();
  parse("PRESSURE_MIN",pres_min);
  parse("PRESSURE_MAX",pres_max);
  unsigned pres_steps=0;
  parse("PRESSURE_STEPS",pres_steps);
  parseVector("PRESSURE_SET_ALL",pres_);
//other
  std::vector<double> cut_corner;
  parseVector("CUT_CORNER",cut_corner);
  parseVector("SET_ALL_TEMP_PRESSURE",custom_lambdas_);

  checkRead();

  if(custom_lambdas_.size()>0)
  {
    //make sure no incompatible options are used
    plumed_massert(temps.size()==0,"cannot set both SET_ALL_TEMP_PRESSURE and TEMP_SET_ALL");
    plumed_massert(pres_.size()==0,"cannot set both SET_ALL_TEMP_PRESSURE and PRESSURE_SET_ALL");
    plumed_massert(temp_steps==0,"cannot set both SET_ALL_TEMP_PRESSURE and TEMP_STEPS");
    plumed_massert(pres_steps==0,"cannot set both SET_ALL_TEMP_PRESSURE and PRESSURE_STEPS");
    plumed_massert(temp_min==-1 && temp_max==-1,"cannot set both SET_ALL_TEMP_PRESSURE and TEMP_MIN/MAX");
    plumed_massert(std::isnan(pres_min) && std::isnan(pres_max),"cannot set both SET_ALL_TEMP_PRESSURE and PRESSURE_MIN/MAX");
    plumed_massert(cut_corner.size()==0,"cannot set both SET_ALL_TEMP_PRESSURE and CUT_CORNER");
//setup the target temperature-pressure grid
    derECVs_beta_.resize(custom_lambdas_.size());
    derECVs_pres_.resize(custom_lambdas_.size());
    const std::string error_msg="SET_ALL_TEMP_PRESSURE: two underscore-separated values are expected for each comma-separated point, cannot understand: ";
    for(unsigned i=0; i<custom_lambdas_.size(); i++)
    {
      try
      {
        std::size_t pos1;
        const double temp_i=std::stod(custom_lambdas_[i],&pos1);
        plumed_massert(pos1+1<custom_lambdas_[i].size(),error_msg+custom_lambdas_[i]);
        plumed_massert(custom_lambdas_[i][pos1]=='_',error_msg+custom_lambdas_[i]);
        std::size_t pos2;
        const double pres_i=std::stod(custom_lambdas_[i].substr(pos1+1),&pos2);
        plumed_massert(pos1+1+pos2==custom_lambdas_[i].size(),error_msg+custom_lambdas_[i]);

        derECVs_beta_[i]=(temp0/temp_i-1.)/kbt_;
        derECVs_pres_[i]=(temp0/temp_i*pres_i-pres0_)/kbt_;
      }
      catch (std::exception &ex)
      {
        plumed_merror(error_msg+custom_lambdas_[i]);
      }
    }
  }
  else
  {
    //set the intermediate temperatures
    if(temps.size()>0)
    {
      plumed_massert(temp_steps==0,"cannot set both TEMP_STEPS and TEMP_SET_ALL");
      plumed_massert(temp_min==-1 && temp_max==-1,"cannot set both TEMP_SET_ALL and TEMP_MIN/MAX");
      plumed_massert(temps.size()>=2,"set at least 2 temperatures");
      temp_min=temps[0];
      temp_max=temps[temps.size()-1];
      derECVs_beta_.resize(temps.size());
      for(unsigned k=0; k<derECVs_beta_.size(); k++)
      {
        derECVs_beta_[k]=(temp0/temps[k]-1.)/kbt_;
        if(k<derECVs_beta_.size()-1)
          plumed_massert(temps[k]<=temps[k+1],"TEMP_SET_ALL must be properly ordered");
      }
    }
    else
    { //get TEMP_MIN and TEMP_MAX
      if(temp_min==-1)
      {
        temp_min=temp0;
        log.printf("  no TEMP_MIN provided, using TEMP_MIN=TEMP\n");
      }
      if(temp_max==-1)
      {
        temp_max=temp0;
        log.printf("  no TEMP_MAX provided, using TEMP_MAX=TEMP\n");
      }
      plumed_massert(temp_max>=temp_min,"TEMP_MAX should be bigger than TEMP_MIN");
      derECVs_beta_.resize(2);
      derECVs_beta_[0]=(temp0/temp_min-1.)/kbt_;
      derECVs_beta_[1]=(temp0/temp_max-1.)/kbt_;
      if(temp_min==temp_max && temp_steps==0)
        temp_steps=1;
      if(temp_steps>0)
        derECVs_beta_=getSteps(derECVs_beta_[0],derECVs_beta_[1],temp_steps,"TEMP",geom_spacing_,1./kbt_);
      else
        todoAutomatic_beta_=true;
    }
    const double tol=1e-3; //if temp is taken from MD engine it might be numerically slightly different
    if(temp0<(1-tol)*temp_min || temp0>(1+tol)*temp_max)
      log.printf(" +++ WARNING +++ running at TEMP=%g which is outside the chosen temperature range\n",temp0);

    //set the intermediate pressures
    if(pres_.size()>0)
    {
      plumed_massert(pres_steps==0,"cannot set both PRESSURE_STEPS and PRESSURE_SET_ALL");
      plumed_massert(std::isnan(pres_min) && std::isnan(pres_max),"cannot set both PRESSURE_SET_ALL and PRESSURE_MIN/MAX");
      plumed_massert(pres_.size()>=2,"set at least 2 pressures");
      for(unsigned kk=0; kk<pres_.size()-1; kk++)
        plumed_massert(pres_[kk]<=pres_[kk+1],"PRESSURE_SET_ALL must be properly ordered");
      pres_min=pres_[0];
      pres_max=pres_[pres_.size()-1];
    }
    else
    { //get PRESSURE_MIN and PRESSURE_MAX
      if(std::isnan(pres_min))
      {
        pres_min=pres0_;
        log.printf("  no PRESSURE_MIN provided, using PRESSURE_MIN=PRESSURE\n");
      }
      if(std::isnan(pres_max))
      {
        pres_max=pres0_;
        log.printf("  no PRESSURE_MAX provided, using PRESSURE_MAX=PRESSURE\n");
      }
      plumed_massert(pres_max>=pres_min,"PRESSURE_MAX should be bigger than PRESSURE_MIN");
      if(pres_min==pres_max && pres_steps==0)
        pres_steps=1;
      if(pres_steps>0)
        pres_=getSteps(pres_min,pres_max,pres_steps,"PRESSURE",false,0);
      else
      {
        pres_.resize(2);
        pres_[0]=pres_min;
        pres_[1]=pres_max;
        todoAutomatic_pres_=true;
      }
    }
    if(pres0_<pres_min || pres0_>pres_max)
      log.printf(" +++ WARNING +++ running at PRESSURE=%g which is outside the chosen pressure range\n",pres0_);

    //set CUT_CORNER
    std::string cc_usage("CUT_CORNER=temp_low,pres_low,temp_high,pres_high");
    if(cut_corner.size()==4)
    {
      const double temp_low=cut_corner[0];
      const double pres_low=cut_corner[1];
      const double temp_high=cut_corner[2];
      const double pres_high=cut_corner[3];
      plumed_massert(temp_low<temp_high,"temp_low="+std::to_string(temp_low)+" should be smaller than temp_high="+std::to_string(temp_high)+", "+cc_usage);
      plumed_massert(temp_low>=temp_min && temp_low<=temp_max,"temp_low="+std::to_string(temp_low)+" is out of temperature range. "+cc_usage);
      plumed_massert(temp_high>=temp_min && temp_high<=temp_max,"temp_high="+std::to_string(temp_high)+" is out of temperature range. "+cc_usage);
      plumed_massert(pres_low<pres_high,"pres_low="+std::to_string(pres_low)+" should be smaller than pres_high="+std::to_string(pres_high)+", "+cc_usage);
      plumed_massert(pres_low>=pres_min && pres_low<=pres_max,"pres_low="+std::to_string(pres_low)+" is out of pressure range. "+cc_usage);
      plumed_massert(pres_high>=pres_min && pres_high<=pres_max,"pres_high="+std::to_string(pres_high)+" is out of pressure range. "+cc_usage);
      kB_temp_low_=kB*temp_low;
      coeff_=(pres_high-pres_low)/(temp_high-temp_low)/kB;
      plumed_massert(coeff_!=0,"this should not be possible");
      const double small_value=(temp_high-pres_low)/1e4;
      pres_low_=pres_low-small_value; //make sure pres_max is included
      plumed_massert(pres_max>=coeff_*(kB*temp_max-kB_temp_low_)+pres_low_,"please chose a pres_high slightly smaller than PRESSURE_MAX in "+cc_usage);
    }
    else
    {
      plumed_massert(cut_corner.size()==0,"expected 4 values: "+cc_usage);
    }
  }

//print some info
  log.printf("  running at TEMP=%g and PRESSURE=%g\n",temp0,pres0_);
  log.printf("  targeting a temperature range from TEMP_MIN=%g to TEMP_MAX=%g\n",temp_min,temp_max);
  if(temp_min==temp_max)
    log.printf(" +++ WARNING +++ if you only need a multibaric simulation it is more efficient to set it up with ECV_LINEAR\n");
  log.printf("   and a pressure range from PRESSURE_MIN=%g to PRESSURE_MAX=%g\n",pres_min,pres_max);
  if(pres_min==pres_max)
    log.printf(" +++ WARNING +++ if you only need a multithermal simulation it is more efficient to set it up with ECV_MULTITHERMAL\n");
  if(geom_spacing_)
    log.printf(" -- NO_GEOM_SPACING: inverse temperatures will be linearly spaced\n");
  if(coeff_!=0)
    log.printf(" -- CUT_CORNER: ignoring some high temperature and low pressure values\n");
}

void ECVmultiThermalBaric::calculateECVs(const double * ene_vol)
{
  for(unsigned k=0; k<derECVs_beta_.size(); k++)
    ECVs_beta_[k]=derECVs_beta_[k]*ene_vol[0];
  for(unsigned i=0; i<derECVs_pres_.size(); i++)
    ECVs_pres_[i]=derECVs_pres_[i]*ene_vol[1];
// derivatives are constant, as usual in linear expansions
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
  if(custom_lambdas_.size()>0)
  { //same as default getIndex_k() function
    plumed_massert(totNumECVs_==custom_lambdas_.size(),"this should not happen");
    for(unsigned i=0; i<totNumECVs_; i++)
      index_k.emplace_back(std::vector<unsigned> {i,i});
  }
  else
  {
    unsigned i=0;
    for(unsigned k=0; k<derECVs_beta_.size(); k++)
    {
      const double kB_temp_k=kbt_/(derECVs_beta_[k]*kbt_+1);
      const double line_k=coeff_*(kB_temp_k-kB_temp_low_)+pres_low_;
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
  }
  return index_k;
}

std::vector<std::string> ECVmultiThermalBaric::getLambdas() const
{
  if(custom_lambdas_.size()>0)
    return custom_lambdas_;

  plumed_massert(!todoAutomatic_beta_ && !todoAutomatic_pres_,"cannot access lambdas before initializing them");
  std::vector<std::string> lambdas;
  const double kB=plumed.getAtoms().getKBoltzmann();
  for(unsigned k=0; k<derECVs_beta_.size(); k++)
  {
    const double kB_temp_k=kbt_/(derECVs_beta_[k]*kbt_+1);
    const double line_k=coeff_*(kB_temp_k-kB_temp_low_)+pres_low_;
    for(unsigned kk=0; kk<pres_.size(); kk++)
    {
      if(coeff_==0 || pres_[kk]>=line_k)
      {
        std::ostringstream subs;
        subs<<kB_temp_k/kB<<"_"<<pres_[kk];
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
  if(custom_lambdas_.size()>0)
  {
    log.printf("  *%4lu temperatures for %s\n",derECVs_beta_.size(),getName().c_str());
    log.printf("  *%4lu beta-pressures for %s\n",derECVs_pres_.size(),getName().c_str());
    log.printf("    -- SET_ALL_TEMP_PRESSURE: total number of temp-pres points is %u\n",totNumECVs_);
  }
  else
  {
    plumed_massert(derECVs_beta_.size()*pres_.size()>=totNumECVs_,"this should not happen, is something wrong with CUT_CORNER ?");
    derECVs_pres_.resize(totNumECVs_); //pres is mixed with temp (beta*p*V), thus we need to store all possible
    //initialize the derECVs.
    //this could be done before and one could avoid storing also beta0, beta_k, etc. but this way the code should be more readable
    unsigned i=0;
    for(unsigned k=0; k<derECVs_beta_.size(); k++)
    {
      const double kB_temp_k=kbt_/(derECVs_beta_[k]*kbt_+1);
      const double line_k=coeff_*(kB_temp_k-kB_temp_low_)+pres_low_;
      for(unsigned kk=0; kk<pres_.size(); kk++)
      {
        if(coeff_==0 || pres_[kk]>=line_k)
        {
          derECVs_pres_[i]=(pres_[kk]/kB_temp_k-pres0_/kbt_);
          i++;
        }
      }
    }
    log.printf("  *%4lu temperatures for %s\n",derECVs_beta_.size(),getName().c_str());
    log.printf("  *%4lu pressures for %s\n",pres_.size(),getName().c_str());
    if(coeff_!=0)
      log.printf("    -- CUT_CORNER: %lu temp-pres points were excluded, thus total is %u\n",derECVs_beta_.size()*pres_.size()-totNumECVs_,totNumECVs_);
  }
  ECVs_beta_.resize(derECVs_beta_.size());
  ECVs_pres_.resize(derECVs_pres_.size());
  isReady_=true;
}

void ECVmultiThermalBaric::initECVs_observ(const std::vector<double>& all_obs_cvs,const unsigned ncv,const unsigned index_j)
{
  if(todoAutomatic_beta_) //estimate the steps in beta from observations
  {
    plumed_massert(all_obs_cvs.size()%ncv==0 && index_j<ncv,"initECVs_observ parameters are inconsistent");
    std::vector<double> obs_ene(all_obs_cvs.size()/ncv); //copy only useful observations
    for(unsigned t=0; t<obs_ene.size(); t++)
      obs_ene[t]=all_obs_cvs[t*ncv+index_j]+pres0_*all_obs_cvs[t*ncv+index_j+1]; //U=E+pV
    const unsigned temp_steps=estimateNumSteps(derECVs_beta_[0],derECVs_beta_[1],obs_ene,"TEMP");
    log.printf("    (spacing is on beta, not on temperature)\n");
    derECVs_beta_=getSteps(derECVs_beta_[0],derECVs_beta_[1],temp_steps,"TEMP",geom_spacing_,1./kbt_);
    todoAutomatic_beta_=false;
  }
  if(todoAutomatic_pres_) //estimate the steps in pres from observations
  {
    plumed_massert(all_obs_cvs.size()%ncv==0 && index_j+1<ncv,"initECVs_observ parameters are inconsistent");
    std::vector<double> obs_vol(all_obs_cvs.size()/ncv); //copy only useful observations
    for(unsigned t=0; t<obs_vol.size(); t++)
      obs_vol[t]=all_obs_cvs[t*ncv+index_j+1];
    const unsigned pres_steps=estimateNumSteps((pres_[0]-pres0_)/kbt_,(pres_[1]-pres0_)/kbt_,obs_vol,"PRESSURE");
    log.printf("    (spacing is in beta0 units)\n");
    pres_=getSteps(pres_[0],pres_[1],pres_steps,"PRESSURE",false,0);
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
    unsigned pres_steps=1;
    std::string pres_min=getPres(0);
    for(unsigned i=1; i<lambdas.size(); i++) //pres is second, thus increas by 1
    {
      if(getPres(i)==pres_min)
        break;
      pres_steps++;
    }
    pres_=getSteps(pres_[0],pres_[1],pres_steps,"PRESSURE",false,0);
    todoAutomatic_pres_=false;
  }
  if(todoAutomatic_beta_)
  {
    unsigned temp_steps=1;
    std::string pres_max=getPres(pres_.size()-1);
    for(unsigned i=pres_.size(); i<lambdas.size(); i++)
    { //even if CUT_CORNER, the max pressures are all present, for each temp
      if(getPres(i)==pres_max)
        temp_steps++;
    }
    derECVs_beta_=getSteps(derECVs_beta_[0],derECVs_beta_[1],temp_steps,"TEMP",geom_spacing_,1./kbt_);
    todoAutomatic_beta_=false;
  }
  std::vector<std::string> myLambdas=getLambdas();
  plumed_assert(myLambdas.size()==lambdas.size())<<"RESTART - mismatch in number of "<<getName()<<".\nFrom "<<lambdas.size()<<" labels "<<derECVs_beta_.size()<<" temperatures and "<<pres_.size()<<" pressures were found, for a total of "<<myLambdas.size()<<" estimated steps.\nCheck if the CUT_CORNER or the SET_ALL_TEMP_PRESSURE options are consistent\n";
  plumed_massert(std::equal(myLambdas.begin(),myLambdas.end(),lambdas.begin()),"RESTART - mismatch in lambda values of "+getName());

  initECVs();
}

}
}
