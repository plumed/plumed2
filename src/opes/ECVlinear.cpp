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

namespace PLMD {
namespace opes {

//+PLUMEDOC EXPANSION_CV ECV_LINEAR
/*
Linear expansion, according to a parameter lambda.

This can be used e.g. for thermodynamic integration, or for multibaric simulations, in which case lambda=pressure.
It can also be used for multithermal simulations, but for simplicity it is more convenient to use \ref ECV_MULTITHERMAL.

The difference in Hamiltonian \f$\Delta U\f$ is expected as ARG.
\f[
  \Delta u_\lambda=\beta \lambda \Delta U\, .
\f]
Use the DIMENSIONLESS flag to avoid multiplying for the inverse temperature \f$\beta\f$.

\par Examples

Typical multibaric simulation:

\plumedfile
vol: VOLUME
ecv: ECV_LINEAR ...
  ARG=vol
  TEMP=300
  LAMBDA=0.06022140857*2000 #2 kbar
  MIN_LAMBDA=0.06022140857  #1 bar
  MAX_LAMBDA=0.06022140857*4000 #4 kbar
...
opes: OPES_EXPANDED ARG=ecv.vol PACE=500
\endplumedfile

Typical thermodynamic integration:

\plumedfile
DeltaU: EXTRACV NAME=energy_difference
ecv: ECV_LINEAR ARG=DeltaU TEMP=300
opes: OPES_EXPANDED ARG=ecv.* PACE=100
\endplumedfile

Notice that by defauly LAMBDA=0, MIN_LAMBDA=0 and MAX_LAMBDA=1, which is the typical case for thermodynamic integration.

*/
//+ENDPLUMEDOC

class ECVlinear :
  public ExpansionCVs
{
private:
  bool todoAutomatic_;
  double beta0_;
  double lambda0_;
  std::vector<double> lambda_;
  std::vector<double> ECVs_;
  std::vector<double> derECVs_;
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
  keys.add("compulsory","ARG","the label of the Hamiltonian difference \\f$\\Delta U\\f$");
  keys.add("compulsory","LAMBDA","0","the lambda at which the underlying simulation runs");
  keys.add("optional","MIN_LAMBDA","( default=0 ) the minimum of the lambda range");
  keys.add("optional","MAX_LAMBDA","( default=1 ) the maximum of the lambda range");
  keys.add("optional","STEPS_LAMBDA","uniformly place the lambda values, for a total of STEPS_LAMBDA");
  keys.add("optional","SET_ALL_LAMBDAS","manually set all the lamdbas");
  keys.addFlag("DIMENSIONLESS",false,"ARG is considered dimensionless rather than an energy, thus is not multiplied by \\f$\\beta\\f$");
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
  double min_lambda=std::numeric_limits<double>::quiet_NaN();
  double max_lambda=std::numeric_limits<double>::quiet_NaN();
  parse("MIN_LAMBDA",min_lambda);
  parse("MAX_LAMBDA",max_lambda);
  unsigned steps_lambda=0;
  parse("STEPS_LAMBDA",steps_lambda);
  parseVector("SET_ALL_LAMBDAS",lambda_);

  checkRead();

//set the lambdas
  if(lambda_.size()>0)
  {
    plumed_massert(steps_lambda==0,"cannot set both STEPS_LAMBDA and SET_ALL_LAMBDAS");
    plumed_massert(std::isnan(min_lambda) && std::isnan(max_lambda),"cannot set both SET_ALL_LAMBDAS and MIN/MAX_LAMBDA");
    plumed_massert(lambda_.size()>=2,"set at least 2 lambdas with SET_ALL_LAMBDAS");
    for(unsigned k=0; k<lambda_.size()-1; k++)
      plumed_massert(lambda_[k]<=lambda_[k+1],"SET_ALL_LAMBDAS must be properly ordered");
    min_lambda=lambda_[0];
    max_lambda=lambda_[lambda_.size()-1];
  }
  else
  { //get MIN_LAMBDA and MAX_LAMBDA
    if(std::isnan(min_lambda))
    {
      min_lambda=0;
      log.printf("  no MIN_LAMBDA provided, using MIN_LAMBDA = %g\n",min_lambda);
    }
    if(std::isnan(max_lambda))
    {
      max_lambda=1;
      log.printf("  no MAX_LAMBDA provided, using MAX_LAMBDA = %g\n",max_lambda);
    }
    plumed_massert(max_lambda>=min_lambda,"MAX_LAMBDA should be bigger than MIN_LAMBDA");
    lambda_.resize(2);
    lambda_[0]=min_lambda;
    lambda_[1]=max_lambda;
    if(min_lambda==max_lambda && steps_lambda==0)
      steps_lambda=1;
    if(steps_lambda>0)
      setSteps(lambda_,steps_lambda,"LAMBDA");
    else
      todoAutomatic_=true;
  }
  if(lambda0_<min_lambda || lambda0_>max_lambda)
    log.printf(" +++ WARNING +++ running at LAMBDA=%g which is outside the chosen lambda range\n",lambda0_);

//print some info
  log.printf("  running at LAMBDA=%g\n",lambda0_);
  log.printf("  targeting a lambda range from MIN_LAMBDA=%g to MAX_LAMBDA=%g\n",min_lambda,max_lambda);
  if(dimensionless)
    log.printf(" -- DIMENSIONLESS: the ARG is not multiplied by beta\n");
}

void ECVlinear::calculateECVs(const double * DeltaU)
{
  for(unsigned k=0; k<lambda_.size(); k++)
  {
    const double diff_k=beta0_*(lambda_[k]-lambda0_);
    ECVs_[k]=diff_k*DeltaU[0];
    derECVs_[k]=diff_k;
  }
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
  std::vector<std::string> lambdas(lambda_.size());
  for(unsigned k=0; k<lambda_.size(); k++)
  {
    std::ostringstream subs;
    subs<<lambda_[k];
    lambdas[k]=subs.str();
  }
  return lambdas;
}

void ECVlinear::initECVs()
{
  plumed_massert(!isReady_,"initialization should not be called twice");
  plumed_massert(!todoAutomatic_,"this should not happen");
  totNumECVs_=lambda_.size();
  ECVs_.resize(lambda_.size());
  derECVs_.resize(lambda_.size());
  isReady_=true;
  log.printf("  *%4lu lambdas for %s\n",lambda_.size(),getName().c_str());
}

void ECVlinear::initECVs_observ(const std::vector<double>& all_obs_cvs,const unsigned ncv,const unsigned index_j)
{
  if(todoAutomatic_) //estimate the steps in lambda from observations
  {
    plumed_massert(all_obs_cvs.size()%ncv==0 && index_j<ncv,"initECVs_observ parameters are inconsistent");
    std::vector<double> obs_cv(all_obs_cvs.size()/ncv); //copy only useful observation (would be better not to copy...)
    for(unsigned t=0; t<obs_cv.size(); t++)
      obs_cv[t]=all_obs_cvs[t*ncv+index_j];
    const unsigned steps_lambda=estimateSteps(beta0_*(lambda_[0]-lambda0_),beta0_*(lambda_[1]-lambda0_),obs_cv,"LAMBDA");
    if(beta0_!=1)
      log.printf("    (spacing is in beta0 units)\n");
    setSteps(lambda_,steps_lambda,"LAMBDA");
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
    setSteps(lambda_,lambdas.size(),"LAMBDA");
    todoAutomatic_=false;
  }
  std::vector<std::string> myLambdas=getLambdas();
  plumed_massert(myLambdas.size()==lambdas.size(),"RESTART - mismatch in number of "+getName());
  plumed_massert(std::equal(myLambdas.begin(),myLambdas.end(),lambdas.begin()),"RESTART - mismatch in lambda values of "+getName());

  initECVs();
}

}
}
