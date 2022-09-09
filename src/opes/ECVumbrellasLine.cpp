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

//+PLUMEDOC OPES_EXPANSION_CV ECV_UMBRELLAS_LINE
/*
Target a multiumbrella ensemble, by combining systems each with a parabolic bias potential at a different location.

Any set of collective variables \f$\mathbf{s}\f$ can be used as ARG.
\f[
  \Delta u_{\mathbf{s}_i}(\mathbf{s})=\sum_j^{\text{dim}}\frac{([s]_j-[s_i]_j)^2}{2\sigma^2}\, .
\f]
The Gaussian umbrellas are placed along a line, from CV_MIN to CV_MAX.
The umbrellas are placed at a distance SIGMA*SPACING from each other, if you need more flexibility use \ref ECV_UMBRELLAS_FILE.
The unbiased fluctuations in the basin usually are a reasonable guess for the value of SIGMA.
Typically, a SPACING greater than 1 can lead to faster convergence, by reducing the total number of umbrellas.
The umbrellas can be multidimensional, but the CVs dimensions should be rescaled since a single SIGMA must be used.

The keyword BARRIER can be helpful to avoid breaking your system due to a too strong initial bias.
If you think the placed umbrellas will not cover the whole unbiased probability distribution you should add it explicitly to the target, with the flag ADD_P0, for more robust convergence.
See also Appendix B of Ref.\cite Invernizzi2020unified for more details on these last two options.

The flag LOWER_HALF_ONLY modifies the ECVs so that they are set to zero when \f$\mathbf{s}>\mathbf{s}_i\f$, as in \ref LOWER_WALLS.
This can be useful e.g. when the CV used is the \ref ENERGY and one wants to sample a broad range of high energy values, similar to \ref ECV_MULTITHERMAL but with a flat energy distribution as target.
By pushing only from below one can avoid too extreme forces that could crash the simulation.

\par Examples

\plumedfile
cv: DISTANCE ATOMS=1,2
ecv: ECV_UMBRELLAS_LINE ARG=cv CV_MIN=1.2 CV_MAX=4.3 SIGMA=0.5 SPACING=1.5
opes: OPES_EXPANDED ARG=ecv.* PACE=500
\endplumedfile

It is also possible to combine different ECV_UMBRELLAS_LINE to build a grid of CV values that will be sampled.
For example the following code will sample a whole 2D region of cv1 and cv2.

\plumedfile
cv1: DISTANCE ATOMS=1,2
ecv2: ECV_UMBRELLAS_LINE ARG=cv1 CV_MIN=1.2 CV_MAX=4.3 SIGMA=0.5

cv2: DISTANCE ATOMS=3,4
ecv1: ECV_UMBRELLAS_LINE ARG=cv2 CV_MIN=13.8 CV_MAX=21.4 SIGMA=4.3

opes: OPES_EXPANDED ARG=ecv1.*,ecv2.* PACE=500
\endplumedfile

*/
//+ENDPLUMEDOC

class ECVumbrellasLine :
  public ExpansionCVs
{
private:
  double barrier_;
  unsigned P0_contribution_;
  bool lower_only_;

  std::vector< std::vector<double> > centers_;
  double sigma_;

  std::vector< std::vector<double> > ECVs_;
  std::vector< std::vector<double> > derECVs_;
  void initECVs();

public:
  explicit ECVumbrellasLine(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  void calculateECVs(const double *) override;
  const double * getPntrToECVs(unsigned) override;
  const double * getPntrToDerECVs(unsigned) override;
  std::vector<std::string> getLambdas() const override;
  void initECVs_observ(const std::vector<double>&,const unsigned,const unsigned) override;
  void initECVs_restart(const std::vector<std::string>&) override;
};

PLUMED_REGISTER_ACTION(ECVumbrellasLine,"ECV_UMBRELLAS_LINE")

void ECVumbrellasLine::registerKeywords(Keywords& keys)
{
  ExpansionCVs::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","CV_MIN","the minimum of the CV range to be explored");
  keys.add("compulsory","CV_MAX","the maximum of the CV range to be explored");
  keys.add("compulsory","SIGMA","sigma of the umbrella Gaussians");
  keys.add("compulsory","SPACING","1","the distance between umbrellas, in units of SIGMA");
  keys.add("optional","BARRIER","a guess of the free energy barrier to be overcome (better to stay higher than lower)");
  keys.addFlag("ADD_P0",false,"add the unbiased Boltzmann distribution to the target distribution, to make sure to sample it");
  keys.addFlag("LOWER_HALF_ONLY",false,"use only the lower half of each umbrella potentials");
}

ECVumbrellasLine::ECVumbrellasLine(const ActionOptions&ao):
  Action(ao),
  ExpansionCVs(ao)
{
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
  parseFlag("LOWER_HALF_ONLY",lower_only_);

//set umbrellas
  parse("SIGMA",sigma_);
  std::vector<double> cv_min;
  std::vector<double> cv_max;
  parseVector("CV_MIN",cv_min);
  parseVector("CV_MAX",cv_max);
  plumed_massert(cv_min.size()==getNumberOfArguments(),"wrong number of CV_MINs");
  plumed_massert(cv_max.size()==getNumberOfArguments(),"wrong number of CV_MAXs");
  double spacing;
  parse("SPACING",spacing);
  double length=0;
  for(unsigned j=0; j<getNumberOfArguments(); j++)
    length+=std::pow(cv_max[j]-cv_min[j],2);
  length=std::sqrt(length);
  unsigned sizeUmbrellas=1+std::round(length/(sigma_*spacing));
  centers_.resize(getNumberOfArguments()); //centers_[cv][umbrellas]
  unsigned full_period=0;
  for(unsigned j=0; j<getNumberOfArguments(); j++)
  {
    centers_[j].resize(sizeUmbrellas);
    if(sizeUmbrellas>1)
      for(unsigned k=0; k<sizeUmbrellas; k++)
        centers_[j][k]=cv_min[j]+k*(cv_max[j]-cv_min[j])/(sizeUmbrellas-1);
    else
      centers_[j][0]=(cv_min[j]+cv_max[j])/2.;
    if(getPntrToArgument(j)->isPeriodic())
    {
      double min,max;
      std::string min_str,max_str;
      getPntrToArgument(j)->getDomain(min,max);
      getPntrToArgument(j)->getDomain(min_str,max_str);
      plumed_massert(cv_min[j]>=min,"ARG "+std::to_string(j)+": CV_MIN cannot be smaller than the periodic bound "+min_str);
      plumed_massert(cv_max[j]<=max,"ARG "+std::to_string(j)+": CV_MAX cannot be greater than the periodic bound "+max_str);
      if(cv_min[j]==min && cv_max[j]==max)
        full_period++;
    }
  }
  if(full_period==getNumberOfArguments() && sizeUmbrellas>1) //first and last are the same point
  {
    sizeUmbrellas--;
    for(unsigned j=0; j<getNumberOfArguments(); j++)
      centers_[j].pop_back();
  }

  checkRead();

//set ECVs stuff
  totNumECVs_=sizeUmbrellas+P0_contribution_;
  ECVs_.resize(getNumberOfArguments(),std::vector<double>(totNumECVs_));
  derECVs_.resize(getNumberOfArguments(),std::vector<double>(totNumECVs_));

//printing some info
  log.printf("  total number of umbrellas = %u\n",sizeUmbrellas);
  log.printf("    with SIGMA = %g\n",sigma_);
  log.printf("    and SPACING = %g\n",spacing);
  if(barrier_!=std::numeric_limits<double>::infinity())
    log.printf("  guess for free energy BARRIER = %g\n",barrier_);
  if(P0_contribution_==1)
    log.printf(" -- ADD_P0: the target includes also the unbiased probability itself\n");
  if(lower_only_)
    log.printf(" -- LOWER_HALF_ONLY: the ECVs are set to zero for values of the CV above the respective center\n");
}

void ECVumbrellasLine::calculateECVs(const double * cv)
{
  if(lower_only_)
  {
    for(unsigned j=0; j<getNumberOfArguments(); j++)
    {
      for(unsigned k=P0_contribution_; k<totNumECVs_; k++) //if ADD_P0, the first ECVs=0
      {
        const unsigned kk=k-P0_contribution_;
        const double dist_jk=difference(j,centers_[j][kk],cv[j])/sigma_; //PBC might be present
        if(dist_jk>=0)
        {
          ECVs_[j][k]=0;
          derECVs_[j][k]=0;
        }
        else
        {
          ECVs_[j][k]=0.5*std::pow(dist_jk,2);
          derECVs_[j][k]=dist_jk/sigma_;
        }
      }
    }
  }
  else
  {
    for(unsigned j=0; j<getNumberOfArguments(); j++)
    {
      for(unsigned k=P0_contribution_; k<totNumECVs_; k++) //if ADD_P0, the first ECVs=0
      {
        const unsigned kk=k-P0_contribution_;
        const double dist_jk=difference(j,centers_[j][kk],cv[j])/sigma_; //PBC might be present
        ECVs_[j][k]=0.5*std::pow(dist_jk,2);
        derECVs_[j][k]=dist_jk/sigma_;
      }
    }
  }
}

const double * ECVumbrellasLine::getPntrToECVs(unsigned j)
{
  plumed_massert(isReady_,"cannot access ECVs before initialization");
  plumed_massert(j<getNumberOfArguments(),getName()+" has fewer CVs");
  return &ECVs_[j][0];
}

const double * ECVumbrellasLine::getPntrToDerECVs(unsigned j)
{
  plumed_massert(isReady_,"cannot access ECVs before initialization");
  plumed_massert(j<getNumberOfArguments(),getName()+" has fewer CVs");
  return &derECVs_[j][0];
}

std::vector<std::string> ECVumbrellasLine::getLambdas() const
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
    subs<<centers_[0][kk];
    for(unsigned j=1; j<getNumberOfArguments(); j++)
      subs<<"_"<<centers_[j][kk];
    lambdas[k]=subs.str();
  }
  return lambdas;
}

void ECVumbrellasLine::initECVs()
{
  plumed_massert(!isReady_,"initialization should not be called twice");
  isReady_=true;
  log.printf("  *%4u windows for %s\n",totNumECVs_,getName().c_str());
}

void ECVumbrellasLine::initECVs_observ(const std::vector<double>& all_obs_cvs,const unsigned ncv,const unsigned index_j)
{
  //this non-linear exansion never uses automatic initialization
  initECVs();
  calculateECVs(&all_obs_cvs[index_j]); //use only first obs point
  for(unsigned j=0; j<getNumberOfArguments(); j++)
    for(unsigned k=P0_contribution_; k<totNumECVs_; k++)
      ECVs_[j][k]=std::min(barrier_/kbt_,ECVs_[j][k]);
}

void ECVumbrellasLine::initECVs_restart(const std::vector<std::string>& lambdas)
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
