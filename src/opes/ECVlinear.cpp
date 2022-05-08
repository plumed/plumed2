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

//+PLUMEDOC OPES_EXPANSION_CV ECV_LINEAR
/*
Linear expansion, according to a parameter lambda.

This can be used e.g. for thermodynamic integration, or for multibaric simulations, in which case lambda=pressure.
It can also be used for multithermal simulations, but for simplicity it is more convenient to use \ref ECV_MULTITHERMAL.

The difference in Hamiltonian \f$\Delta U\f$ is expected as ARG.
\f[
  \Delta u_\lambda=\beta \lambda \Delta U\, .
\f]
Use the DIMENSIONLESS flag to avoid multiplying for the inverse temperature \f$\beta\f$.

By defauly the needed steps in lambda are automatically guessed from few initial unbiased MD steps, as descibed in \cite Invernizzi2020unified.
Otherwise one can set this number with LAMBDA_STEPS.
In both cases the steps will be uniformly distriuted.
Finally, one can use instead the keyword LAMBDA_SET_ALL and explicitly provide each lambda value.

\par Examples

Typical multibaric simulation:

\plumedfile
vol: VOLUME
ecv: ECV_LINEAR ...
  ARG=vol
  TEMP=300
  LAMBDA=0.06022140857*2000 #2 kbar
  LAMBDA_MIN=0.06022140857  #1 bar
  LAMBDA_MAX=0.06022140857*4000 #4 kbar
...
opes: OPES_EXPANDED ARG=ecv.vol PACE=500
\endplumedfile

Typical thermodynamic integration:

\plumedfile
DeltaU: EXTRACV NAME=energy_difference
ecv: ECV_LINEAR ARG=DeltaU TEMP=300
opes: OPES_EXPANDED ARG=ecv.* PACE=100
\endplumedfile

Notice that by defauly LAMBDA=0, LAMBDA_MIN=0 and LAMBDA_MAX=1, which is the typical case for thermodynamic integration.

*/
//+ENDPLUMEDOC

class ECVlinear :
  public ExpansionCVs
{
private:
  bool todoAutomatic_;
  bool geom_spacing_;
  double beta0_;
  double lambda0_;
  std::vector<double> ECVs_;
  std::vector<double> derECVs_; //beta0*(lambda_k-lambda0)
  void initECVs();

public:
  explicit ECVlinear(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  void calculateECVs(const double *) override;
  const double * getPntrToECVs(unsigned) override;
  const double * getPntrToDerECVs(unsigned) override;
  std::vector<std::string> getLambdas() const override;
  void initECVs_observ(const std::vector<double>&,const unsigned,const unsigned) override;
  void initECVs_restart(const std::vector<std::string>&) override;
};

PLUMED_REGISTER_ACTION(ECVlinear,"ECV_LINEAR")

void ECVlinear::registerKeywords(Keywords& keys)
{
  ExpansionCVs::registerKeywords(keys);
  keys.remove("ARG");
  keys.add("compulsory","ARG","the label of the Hamiltonian difference. \\f$\\Delta U\\f$");
  keys.add("compulsory","LAMBDA","0","the lambda at which the underlying simulation runs");
  keys.add("optional","LAMBDA_MIN","( default=0 ) the minimum of the lambda range");
  keys.add("optional","LAMBDA_MAX","( default=1 ) the maximum of the lambda range");
  keys.add("optional","LAMBDA_STEPS","uniformly place the lambda values, for a total of LAMBDA_STEPS");
  keys.add("optional","LAMBDA_SET_ALL","manually set all the lamdbas");
  keys.addFlag("DIMENSIONLESS",false,"ARG is considered dimensionless rather than an energy, thus is not multiplied by \\f$\\beta\\f$");
  keys.addFlag("GEOM_SPACING",false,"use geometrical spacing in lambda instead of linear spacing");
}

ECVlinear::ECVlinear(const ActionOptions&ao)
  : Action(ao)
  , ExpansionCVs(ao)
  , todoAutomatic_(false)
  , beta0_(1./kbt_)
{
  plumed_massert(getNumberOfArguments()==1,"only DeltaU should be given as ARG");

//set beta0_
  bool dimensionless;
  parseFlag("DIMENSIONLESS",dimensionless);
  if(dimensionless)
    beta0_=1;

//parse lambda info
  parse("LAMBDA",lambda0_);
  double lambda_min=std::numeric_limits<double>::quiet_NaN();
  double lambda_max=std::numeric_limits<double>::quiet_NaN();
  parse("LAMBDA_MIN",lambda_min);
  parse("LAMBDA_MAX",lambda_max);
  unsigned lambda_steps=0;
  parse("LAMBDA_STEPS",lambda_steps);
  std::vector<double> lambdas;
  parseVector("LAMBDA_SET_ALL",lambdas);
  parseFlag("GEOM_SPACING",geom_spacing_);

  checkRead();

//set the diff vector using lambdas
  if(lambdas.size()>0)
  {
    plumed_massert(lambda_steps==0,"cannot set both LAMBDA_STEPS and LAMBDA_SET_ALL");
    plumed_massert(std::isnan(lambda_min) && std::isnan(lambda_max),"cannot set both LAMBDA_SET_ALL and LAMBDA_MIN/MAX");
    plumed_massert(lambdas.size()>=2,"set at least 2 lambdas with LAMBDA_SET_ALL");
    for(unsigned k=0; k<lambdas.size()-1; k++)
      plumed_massert(lambdas[k]<=lambdas[k+1],"LAMBDA_SET_ALL must be properly ordered");
    lambda_min=lambdas[0];
    lambda_max=lambdas[lambdas.size()-1];
    derECVs_.resize(lambdas.size());
    for(unsigned k=0; k<derECVs_.size(); k++)
      derECVs_[k]=beta0_*(lambdas[k]-lambda0_);
  }
  else
  { //get LAMBDA_MIN and LAMBDA_MAX
    if(std::isnan(lambda_min))
    {
      lambda_min=0;
      log.printf("  no LAMBDA_MIN provided, using LAMBDA_MIN = %g\n",lambda_min);
    }
    if(std::isnan(lambda_max))
    {
      lambda_max=1;
      log.printf("  no LAMBDA_MAX provided, using LAMBDA_MAX = %g\n",lambda_max);
    }
    plumed_massert(lambda_max>=lambda_min,"LAMBDA_MAX should be bigger than LAMBDA_MIN");
    derECVs_.resize(2);
    derECVs_[0]=beta0_*(lambda_min-lambda0_);
    derECVs_[1]=beta0_*(lambda_max-lambda0_);
    if(lambda_min==lambda_max && lambda_steps==0)
      lambda_steps=1;
    if(lambda_steps>0)
      derECVs_=getSteps(derECVs_[0],derECVs_[1],lambda_steps,"LAMBDA",geom_spacing_,beta0_*lambda0_);
    else
      todoAutomatic_=true;
  }
  if(lambda0_<lambda_min || lambda0_>lambda_max)
    log.printf(" +++ WARNING +++ running at LAMBDA=%g which is outside the chosen lambda range\n",lambda0_);

//print some info
  log.printf("  running at LAMBDA=%g\n",lambda0_);
  log.printf("  targeting a lambda range from LAMBDA_MIN=%g to LAMBDA_MAX=%g\n",lambda_min,lambda_max);
  if(dimensionless)
    log.printf(" -- DIMENSIONLESS: the ARG is not multiplied by beta\n");
  if(geom_spacing_)
    log.printf(" -- GEOM_SPACING: lambdas will be geometrically spaced\n");
}

void ECVlinear::calculateECVs(const double * DeltaU)
{
  for(unsigned k=0; k<derECVs_.size(); k++)
    ECVs_[k]=derECVs_[k]*DeltaU[0];
// derivatives never change: derECVs_k=beta0*(lambda_k-lambda0)
}

const double * ECVlinear::getPntrToECVs(unsigned j)
{
  plumed_massert(isReady_,"cannot access ECVs before initialization");
  plumed_massert(j==0,getName()+" has only one CV, the DeltaU");
  return &ECVs_[0];
}

const double * ECVlinear::getPntrToDerECVs(unsigned j)
{
  plumed_massert(isReady_,"cannot access ECVs before initialization");
  plumed_massert(j==0,getName()+" has only one CV, the DeltaU");
  return &derECVs_[0];
}

std::vector<std::string> ECVlinear::getLambdas() const
{
  plumed_massert(!todoAutomatic_,"cannot access lambdas before initializing them");
  std::vector<std::string> lambdas(derECVs_.size());
  for(unsigned k=0; k<derECVs_.size(); k++)
  {
    std::ostringstream subs;
    subs<<derECVs_[k]/beta0_+lambda0_; //lambda_k
    lambdas[k]=subs.str();
  }
  return lambdas;
}

void ECVlinear::initECVs()
{
  plumed_massert(!isReady_,"initialization should not be called twice");
  plumed_massert(!todoAutomatic_,"this should not happen");
  totNumECVs_=derECVs_.size();
  ECVs_.resize(derECVs_.size());
  isReady_=true;
  log.printf("  *%4lu lambdas for %s\n",derECVs_.size(),getName().c_str());
}

void ECVlinear::initECVs_observ(const std::vector<double>& all_obs_cvs,const unsigned ncv,const unsigned index_j)
{
  if(todoAutomatic_) //estimate the steps in lambda from observations
  {
    plumed_massert(all_obs_cvs.size()%ncv==0 && index_j<ncv,"initECVs_observ parameters are inconsistent");
    std::vector<double> obs_cv(all_obs_cvs.size()/ncv); //copy only useful observation (would be better not to copy...)
    for(unsigned t=0; t<obs_cv.size(); t++)
      obs_cv[t]=all_obs_cvs[t*ncv+index_j];
    const unsigned lambda_steps=estimateNumSteps(derECVs_[0],derECVs_[1],obs_cv,"LAMBDA");
    if(beta0_!=1)
      log.printf("    (spacing is in beta0 units)\n");
    derECVs_=getSteps(derECVs_[0],derECVs_[1],lambda_steps,"LAMBDA",geom_spacing_,beta0_*lambda0_);
    todoAutomatic_=false;
  }
  initECVs();
  calculateECVs(&all_obs_cvs[index_j]);
}

void ECVlinear::initECVs_restart(const std::vector<std::string>& lambdas)
{
  std::size_t pos=lambdas[0].find("_");
  plumed_massert(pos==std::string::npos,"this should not happen, only one CV is used in "+getName());
  if(todoAutomatic_)
  {
    derECVs_=getSteps(derECVs_[0],derECVs_[1],lambdas.size(),"LAMBDA",geom_spacing_,beta0_*lambda0_);
    todoAutomatic_=false;
  }
  std::vector<std::string> myLambdas=getLambdas();
  plumed_massert(myLambdas.size()==lambdas.size(),"RESTART - mismatch in number of "+getName());
  plumed_massert(std::equal(myLambdas.begin(),myLambdas.end(),lambdas.begin()),"RESTART - mismatch in lambda values of "+getName());

  initECVs();
}

}
}
