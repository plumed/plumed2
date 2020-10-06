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
#include "tools/OpenMP.h"
#include "tools/Communicator.h"

namespace PLMD {
namespace opes {

void ExpansionCVs::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  ActionWithValue::useCustomisableComponents(keys);
  keys.add("optional","BARRIER","a guess of the free energy barrier to be overcome (better to stay higher than lower)");
//  keys.reserve("compulsory","PERIODIC","NO","if the output of your ECVs is periodic then you should specify the periodicity.");
}

ExpansionCVs::ExpansionCVs(const ActionOptions&ao)
  : Action(ao)
  , ActionWithValue(ao)
  , ActionWithArguments(ao)
  , isReady_(false)
  , totNumECVs_(0)
{
  plumed_massert( getNumberOfArguments()!=0, "you must specify the underlying CV");
  for(unsigned i=0; i<getNumberOfArguments(); i++)
  {
    std::string name_i=getPntrToArgument(i)->getName();
    ActionWithValue::addComponentWithDerivatives(name_i);
    getPntrToComponent(i)->setNotPeriodic();
    getPntrToComponent(i)->resizeDerivatives(getNumberOfArguments());
  }
  barrier_=std::numeric_limits<double>::infinity();
  parse("BARRIER",barrier_);
  if(barrier_!=std::numeric_limits<double>::infinity())
    log.printf("  guess for free energy BARRIER = %g\n",barrier_);
}

void ExpansionCVs::calculate()
{
  std::vector<double> args(getNumberOfArguments());
  for(unsigned i=0; i<getNumberOfArguments(); i++)
  {
    Value *v_i=getPntrToComponent(i);
    args[i]=getArgument(i);
    v_i->set(args[i]);
    for(unsigned j=0; j<getNumberOfArguments(); j++)
    {
      const double der_ij=(i==j?1:0);
      v_i->addDerivative(j,der_ij);
    }
  }
  if(isReady_)
    calculateECVs(&args[0]);
}

void ExpansionCVs::apply()
{
  const unsigned noa=getNumberOfArguments();
  const unsigned ncp=getNumberOfComponents();
  const unsigned cgs=comm.Get_size();

  std::vector<double> f(noa,0.0);

  unsigned stride=1;
  unsigned rank=0;
  if(ncp>4*cgs) {
    stride=comm.Get_size();
    rank=comm.Get_rank();
  }

  unsigned at_least_one_forced=0;
  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(f)
  {
    std::vector<double> omp_f(noa,0.0);
    std::vector<double> forces(noa);
    #pragma omp for reduction( + : at_least_one_forced)
    for(unsigned i=rank; i<ncp; i+=stride) {
      if(getPntrToComponent(i)->applyForce(forces)) {
        at_least_one_forced+=1;
        for(unsigned j=0; j<noa; j++) omp_f[j]+=forces[j];
      }
    }
    #pragma omp critical
    for(unsigned j=0; j<noa; j++) f[j]+=omp_f[j];
  }

  if(noa>0&&ncp>4*cgs) { comm.Sum(&f[0],noa); comm.Sum(at_least_one_forced); }

  if(at_least_one_forced>0) for(unsigned i=0; i<noa; ++i) getPntrToArgument(i)->addForce(f[i]);
}

unsigned ExpansionCVs::estimate_steps(const double left_side,const double right_side,const std::vector<double>& obs,const std::string msg) const
{ //for linear expansions only, it uses Neff to estimate the grid spacing
  if(left_side==0 && right_side==0)
  {
    log.printf(" +++ WARNING +++ MIN_%s and MAX_%s are equal to %s, using single step\n",msg.c_str(),msg.c_str(),msg.c_str());
    return 1;
  }
  auto get_neff_HWHM=[](const double side,const std::vector<double>& obs,const double av_obs) //HWHM = half width at half maximum. neff is in general not symmetric
  {
    //func: Neff/N-0.5 is a function between -0.5 and 0.5
    auto func=[](const long double delta,const std::vector<double> obs, const double av_obs)
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
    if(b<0)//left side case
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
      func_c=func(c,obs,av_obs);//func is evaluated only here, it might be expensive
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
      log.printf(" --- MIN_%s is equal to %s\n",msg.c_str(),msg.c_str());
    else
      log.printf(" +++ WARNING +++ MIN_%s is very close to %s\n",msg.c_str(),msg.c_str());
  }
  if(right_HWHM==0)
  {
    left_HWHM*=2;
    if(right_side==0)
      log.printf(" --- MAX_%s is equal to %s\n",msg.c_str(),msg.c_str());
    else
      log.printf(" +++ WARNING +++ MAX_%s is very close to %s\n",msg.c_str(),msg.c_str());
  }
  if(left_HWHM==0 && right_HWHM==0)
  {
    log.printf(" +++ WARNING +++ %s range is very narrow, using MIN_%s and MAX_%s as only steps\n",msg.c_str(),msg.c_str(),msg.c_str());
    return 2;
  }
  const double grid_spacing=left_HWHM+right_HWHM;
  log.printf("   estimated %s spacing = %g\n",msg.c_str(),grid_spacing);
  unsigned steps=std::ceil(std::abs(right_side-left_side)/grid_spacing);
  plumed_massert(steps>1,"something went wrong and estimated grid spacing for "+msg+" gives a step="+std::to_string(steps));
  return steps;
}

std::vector< std::vector<unsigned> > ExpansionCVs::getIndex_k() const
{
  plumed_massert(isReady_ && totNumECVs_>0,"cannot access getIndex_k() of ECV before initialization");
  plumed_massert(getNumberOfArguments()==1,"buggy ECV: you should override getIndex_k() if you have more than one ARG");
  std::vector< std::vector<unsigned> > index_k(totNumECVs_,std::vector<unsigned>(1));
  for(unsigned k=0; k<totNumECVs_; k++)
    index_k[k][0]=k; //it is trivial when only one ARG is used
  return index_k;
}

}
}
