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

//+PLUMEDOC EXPANSION_CV ECV_CUSTOM
/*
Use some given CVs as a set of expansion collective variables (ECVs).

This can be useful e.g. for quickly testing new ECVs, but from a performance point of view it is probably better to implement a new ECV class.

By default the ARGs are expeted to be energies, \f$\Delta U_i\f$, and are then multiplied by the inverse temperature \f$\beta\f$
\f[
  \Delta u_i=\beta \Delta U_i\, .
\f]
Use the DIMENSIONLESS flag to avoid this multiplication.

The flag ADD_P0 adds also the unbiased distribution to the target.
It is possible to specify a BARRIER as in \ref ECV_UMBRELLAS_LINE, to avoid a too high initial bias.

\par Examples

\plumedfile
ene: ENERGY
t1: CUSTOM PERIODIC=NO ARG=ene FUNC=(300/500-1)*x
t2: CUSTOM PERIODIC=NO ARG=ene FUNC=(300/1000-1)*x
ecv: ECV_CUSTOM ARG=t1,t2 TEMP=300 ADD_P0
opes: OPES_EXPANDED ARG=ecv.* PACE=500
\endplumedfile

It is equivalent to the following:

\plumedfile
ene: ENERGY
ecv: ECV_MULTITHERMAL ARG=ene TEMP=300 SET_ALL_TEMPS=300,500,1000
opes: OPES_EXPANDED ARG=ecv.* PACE=500
\endplumedfile

*/
//+ENDPLUMEDOC

class ECVcustom :
  public ExpansionCVs
{
private:
  unsigned P0_contribution_;
  double barrier_;
  double beta0_;

  std::vector< std::vector<double> > ECVs_;
  std::vector< std::vector<double> > derECVs_;
  void initECVs();

public:
  explicit ECVcustom(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  void calculateECVs(const double *) override;
  const double * getPntrToECVs(unsigned) override;
  const double * getPntrToDerECVs(unsigned) override;
  std::vector< std::vector<unsigned> > getIndex_k() const override;
  std::vector<std::string> getLambdas() const override;
  void initECVs_observ(const std::vector<double>&,const unsigned,const unsigned) override;
  void initECVs_restart(const std::vector<std::string>&) override;
};

PLUMED_REGISTER_ACTION(ECVcustom,"ECV_CUSTOM")

void ECVcustom::registerKeywords(Keywords& keys)
{
  ExpansionCVs::registerKeywords(keys);
  keys.remove("ARG");
  keys.add("compulsory","ARG","the labels of the single ECVs, \\f$\\Delta U_i\\f$, in energy units");
  keys.addFlag("ADD_P0",false,"add the unbiased Boltzmann distribution to the target distribution, to make sure to sample it");
  keys.addFlag("DIMENSIONLESS",false,"consider ARG as dimensionless rather than an energy, thus do not multiply it by \\f$\\beta\\f$");
  keys.add("optional","BARRIER","a guess of the free energy barrier to be overcome (better to stay higher than lower)");
}

ECVcustom::ECVcustom(const ActionOptions&ao)
  : Action(ao)
  , ExpansionCVs(ao)
  , beta0_(1./kbt_)
{
//set beta0_
  bool dimensionless;
  parseFlag("DIMENSIONLESS",dimensionless);
  if(dimensionless)
    beta0_=1;

//set P0_contribution_
  bool add_P0=false;
  parseFlag("ADD_P0",add_P0);
  if(add_P0)
    P0_contribution_=1;
  else
    P0_contribution_=0;

//set barrier_
  barrier_=std::numeric_limits<double>::infinity();
  parse("BARRIER",barrier_);

  checkRead();

//set ECVs stuff
  totNumECVs_=getNumberOfArguments()+P0_contribution_;
  ECVs_.resize(getNumberOfArguments(),std::vector<double>(totNumECVs_));
  derECVs_.resize(getNumberOfArguments(),std::vector<double>(totNumECVs_));
  for(unsigned j=0; j<getNumberOfArguments(); j++)
    derECVs_[j][j+P0_contribution_]=beta0_; //always constant

//print some info
  if(dimensionless)
    log.printf(" -- DIMENSIONLESS: the ARG is not multiplied by beta\n");
  if(barrier_!=std::numeric_limits<double>::infinity())
  {
    log.printf("  guess for free energy BARRIER = %g\n",barrier_);
    if(dimensionless)
      log.printf("    also the BARRIER is considered to be DIMENSIONLESS\n");
  }
  if(P0_contribution_==1)
    log.printf(" -- ADD_P0: the target includes also the unbiased probability itself\n");
}

void ECVcustom::calculateECVs(const double * cv)
{
  for(unsigned j=0; j<getNumberOfArguments(); j++)
    ECVs_[j][j+P0_contribution_]=beta0_*cv[j];
  //derivative is constant
}

const double * ECVcustom::getPntrToECVs(unsigned j)
{
  plumed_massert(isReady_,"cannot access ECVs before initialization");
  plumed_massert(j<getNumberOfArguments(),getName()+" has fewer CVs");
  return &ECVs_[j][0];
}

const double * ECVcustom::getPntrToDerECVs(unsigned j)
{
  plumed_massert(isReady_,"cannot access ECVs before initialization");
  plumed_massert(j<getNumberOfArguments(),getName()+" has fewer CVs");
  return &derECVs_[j][0];
}

std::vector< std::vector<unsigned> > ECVcustom::getIndex_k() const
{
  plumed_massert(isReady_ && totNumECVs_>0,"cannot access getIndex_k() of ECV before initialization");
  std::vector< std::vector<unsigned> > index_k(totNumECVs_,std::vector<unsigned>(getNumberOfArguments()));
  for(unsigned k=0; k<totNumECVs_; k++)
    for(unsigned j=0; j<getNumberOfArguments(); j++)
      if(k==j+P0_contribution_)
        index_k[k][j]=k;
  return index_k;
}

std::vector<std::string> ECVcustom::getLambdas() const
{
  std::vector<std::string> lambdas(totNumECVs_);
  if(P0_contribution_==1)
  {
    std::ostringstream subs;
    subs<<"P0";
    for(unsigned j=1; j<getNumberOfArguments(); j++)
      subs<<"_P0";
    lambdas[0]=subs.str();
  }
  for(unsigned k=P0_contribution_; k<totNumECVs_; k++)
  {
    const unsigned kk=k-P0_contribution_;
    std::ostringstream subs;
//the getLambdas method is const, so it complains if one tries to access a non-const pointer, hence the const_cast
    if(kk==0)
      subs<<const_cast<ECVcustom *>(this)->getPntrToArgument(kk)->getName();
    else
      subs<<"NaN";
    for(unsigned j=1; j<getNumberOfArguments(); j++)
    {
      if(kk==j)
        subs<<"_"<<const_cast<ECVcustom *>(this)->getPntrToArgument(kk)->getName();
      else
        subs<<"_NaN";
    }
    lambdas[k]=subs.str();
  }
  return lambdas;
}

void ECVcustom::initECVs()
{
  plumed_massert(!isReady_,"initialization should not be called twice");
  isReady_=true;
  log.printf("  *%4u ECVs for %s\n",totNumECVs_,getName().c_str());
}

void ECVcustom::initECVs_observ(const std::vector<double>& all_obs_cvs,const unsigned ncv,const unsigned index_j)
{
  initECVs();
  calculateECVs(&all_obs_cvs[index_j]);
  for(unsigned j=0; j<getNumberOfArguments(); j++)
    ECVs_[j][j+P0_contribution_]=std::min(barrier_*beta0_,ECVs_[j][j+P0_contribution_]);
}

void ECVcustom::initECVs_restart(const std::vector<std::string>& lambdas)
{
  std::size_t pos=0;
  for(unsigned j=0; j<getNumberOfArguments()-1; j++)
    pos=lambdas[0].find("_",pos+1); //checking only lambdas[0] is hopefully enough
  plumed_massert(pos<lambdas[0].length(),"this should not happen, fewer '_' than expected in "+getName());
  pos=lambdas[0].find("_",pos+1);
  plumed_massert(pos>lambdas[0].length(),"this should not happen, more '_' than expected in "+getName());

  std::vector<std::string> myLambdas=getLambdas();
  plumed_massert(myLambdas.size()==lambdas.size(),"RESTART - mismatch in number of "+getName());
  plumed_massert(std::equal(myLambdas.begin(),myLambdas.end(),lambdas.begin()),"RESTART - mismatch in lambda values of "+getName());

  initECVs();
}

}
}
