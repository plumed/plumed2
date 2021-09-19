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

#include "tools/OpenMP.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

namespace PLMD {
namespace opes {

void ExpansionCVs::registerKeywords(Keywords& keys)
{
  Action::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  ActionWithValue::useCustomisableComponents(keys);
  keys.add("compulsory","TEMP","-1","temperature. If not specified tries to get it from MD engine");
}

ExpansionCVs::ExpansionCVs(const ActionOptions&ao)
  : Action(ao)
  , ActionWithValue(ao)
  , ActionWithArguments(ao)
  , isReady_(false)
  , totNumECVs_(0)
{
//set kbt_
  const double kB=plumed.getAtoms().getKBoltzmann();
  kbt_=plumed.getAtoms().getKbT();
  double temp=-1;
  parse("TEMP",temp);
  if(temp>0)
  {
    if(kbt_>0 && std::abs(kbt_-kB*temp)>1e-4)
      log.printf(" +++ WARNING +++ using TEMP=%g while MD engine uses %g\n",temp,kbt_/kB);
    kbt_=kB*temp;
  }
  plumed_massert(kbt_>0,"your MD engine does not pass the temperature to plumed, you must specify it using TEMP");
  log.printf("  temperature = %g, beta = %g\n",kbt_/kB,1./kbt_);

//set components
  plumed_massert( getNumberOfArguments()!=0, "you must specify the underlying CV");
  for(unsigned j=0; j<getNumberOfArguments(); j++)
  {
    std::string name_j=getPntrToArgument(j)->getName();
    ActionWithValue::addComponentWithDerivatives(name_j);
    getPntrToComponent(j)->resizeDerivatives(1);
    if(getPntrToArgument(j)->isPeriodic()) //it should not be necessary, but why not
    {
      std::string min,max;
      getPntrToArgument(j)->getDomain(min,max);
      getPntrToComponent(j)->setDomain(min,max);
    }
    else
      getPntrToComponent(j)->setNotPeriodic();
  }
  plumed_massert((int)getNumberOfArguments()==getNumberOfComponents(),"Expansion CVs have same number of arguments and components");
}

void ExpansionCVs::calculate()
{
  std::vector<double> args(getNumberOfArguments());
  for(unsigned j=0; j<getNumberOfArguments(); j++)
  {
    args[j]=getArgument(j);
    getPntrToComponent(j)->set(args[j]); //components are equal to arguments
    getPntrToComponent(j)->addDerivative(0,1.); //the derivative of the identity is 1
  }
  if(isReady_)
    calculateECVs(&args[0]);
}

void ExpansionCVs::apply()
{
  for(unsigned j=0; j<getNumberOfArguments(); j++)
  {
    std::vector<double> force_j(1);
    if(getPntrToComponent(j)->applyForce(force_j)) //a bias is applied?
      getPntrToArgument(j)->addForce(force_j[0]); //just tell it to the CV!
  }
}

std::vector< std::vector<unsigned> > ExpansionCVs::getIndex_k() const
{
  plumed_massert(isReady_ && totNumECVs_>0,"cannot access getIndex_k() of ECV before initialization");
  std::vector< std::vector<unsigned> > index_k(totNumECVs_,std::vector<unsigned>(getNumberOfArguments()));
  for(unsigned k=0; k<totNumECVs_; k++)
    for(unsigned j=0; j<getNumberOfArguments(); j++)
      index_k[k][j]=k; //each CV gives rise to the same number of ECVs
  return index_k;
}

//following methods are meant to be used only in case of linear expansions
std::vector<double> ExpansionCVs::getSteps(double lambda_min,double lambda_max,const unsigned lambda_steps,const std::string& msg,const bool geom_spacing, const double shift) const
{
  plumed_massert(!(lambda_min==lambda_max && lambda_steps>1),"cannot have multiple "+msg+"_STEPS if "+msg+"_MIN=="+msg+"_MAX");
  std::vector<double> lambda(lambda_steps);
  if(lambda_steps==1)
  {
    lambda[0]=(lambda_min+lambda_max)/2.;
    log.printf(" +++ WARNING +++ using one single %s as target = %g\n",msg.c_str(),lambda[0]);
  }
  else
  {
    if(geom_spacing) //geometric spacing
    { //this way lambda[k]/lambda[k+1] is constant
      lambda_min+=shift;
      lambda_max+=shift;
      plumed_massert(lambda_min>0,"cannot use GEOM_SPACING when %s_MIN is not greater than zero");
      plumed_massert(lambda_max>0,"cannot use GEOM_SPACING when %s_MAX is not greater than zero");
      const double log_lambda_min=std::log(lambda_min);
      const double log_lambda_max=std::log(lambda_max);
      for(unsigned k=0; k<lambda.size(); k++)
        lambda[k]=std::exp(log_lambda_min+k*(log_lambda_max-log_lambda_min)/(lambda_steps-1))-shift;
    }
    else //linear spacing
      for(unsigned k=0; k<lambda.size(); k++)
        lambda[k]=lambda_min+k*(lambda_max-lambda_min)/(lambda_steps-1);
  }
  return lambda;
}

unsigned ExpansionCVs::estimateNumSteps(const double left_side,const double right_side,const std::vector<double>& obs,const std::string& msg) const
{ //for linear expansions only, it uses effective sample size (Neff) to estimate the grid spacing
  if(left_side==0 && right_side==0)
  {
    log.printf(" +++ WARNING +++ %s_MIN and %s_MAX are equal to %s, using single step\n",msg.c_str(),msg.c_str(),msg.c_str());
    return 1;
  }
  auto get_neff_HWHM=[](const double side,const std::vector<double>& obs,const double av_obs) //HWHM = half width at half maximum. neff is in general not symmetric
  {
    //func: Neff/N-0.5 is a function between -0.5 and 0.5
    auto func=[](const long double delta,const std::vector<double>& obs, const double av_obs)
    {
      long double sum_w=0;
      long double sum_w2=0;
      for(unsigned t=0; t<obs.size(); t++)
      {
        const long double w=std::exp(-delta*(obs[t]-av_obs));
        sum_w+=w;
        sum_w2+=w*w;
      }
      return sum_w*sum_w/sum_w2/obs.size()-0.5;
    };
    //here we find the root of func using the regula falsi (false position) method
    //but any method would be OK, not much precision is needed. src/tools/RootFindingBase.h looked complicated
    const double tolerance=1e-4; //seems to be a good default
    double a=0; //default is right side case
    double func_a=0.5;
    double b=side;
    double func_b=func(side,obs,av_obs);
    if(func_b>=0)
      return 0.0; //no zero is present!
    if(b<0) //left side case
    {
      std::swap(a,b);
      std::swap(func_a,func_b);
    }
    double c=a;
    double func_c=func_a;
    while(std::abs(func_c)>tolerance)
    {
      if(func_a*func_c>0)
      {
        a=c;
        func_a=func_c;
      }
      else
      {
        b=c;
        func_b=func_c;
      }
      c=(a*func_b-b*func_a)/(func_b-func_a);
      func_c=func(c,obs,av_obs); //func is evaluated only here, it might be expensive
    }
    return std::abs(c);
  };

//set average to zero, for numerical stability
  double av_obs=0;
  for(unsigned t=0; t<obs.size(); t++)
    av_obs+=obs[t];
  av_obs/=obs.size();

//estimation
  double left_HWHM=0;
  if(left_side!=0)
    left_HWHM=get_neff_HWHM(left_side,obs,av_obs);
  double right_HWHM=0;
  if(right_side!=0)
    right_HWHM=get_neff_HWHM(right_side,obs,av_obs);
  if(left_HWHM==0)
  {
    right_HWHM*=2;
    if(left_side==0)
      log.printf(" --- %s_MIN is equal to %s\n",msg.c_str(),msg.c_str());
    else
      log.printf(" +++ WARNING +++ %s_MIN is very close to %s\n",msg.c_str(),msg.c_str());
  }
  if(right_HWHM==0)
  {
    left_HWHM*=2;
    if(right_side==0)
      log.printf(" --- %s_MAX is equal to %s\n",msg.c_str(),msg.c_str());
    else
      log.printf(" +++ WARNING +++ %s_MAX is very close to %s\n",msg.c_str(),msg.c_str());
  }
  const double grid_spacing=left_HWHM+right_HWHM;
  log.printf("   estimated %s spacing = %g\n",msg.c_str(),grid_spacing);
  unsigned steps=std::ceil(std::abs(right_side-left_side)/grid_spacing);
  if(steps<2 || grid_spacing==0)
  {
    log.printf(" +++ WARNING +++ %s range is very narrow, using %s_MIN and %s_MAX as only steps\n",msg.c_str(),msg.c_str(),msg.c_str());
    steps=2;
  }
  return steps;
}

}
}
