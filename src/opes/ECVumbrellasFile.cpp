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
#include "tools/File.h"

namespace PLMD {
namespace opes {

//+PLUMEDOC OPES_EXPANSION_CV ECV_UMBRELLAS_FILE
/*
Target a multiumbrella ensemble, by combining systems each with a parabolic bias potential at a different location.

Any set of collective variables $\mathbf{s}$ can be used as ARG.
The positions $\mathbf{s}_i$ and dimension $\mathbf{\sigma}_i$ of the umbrellas are read from file.

$$
  \Delta u_{\mathbf{s}_i}(\mathbf{s})=\sum_j^{\text{dim}}\frac{([s]_j-[s_i]_j)^2}{2[\sigma_i]_j^2}\, .
$$

Notice that $\mathbf{\sigma}_i$ is diagonal, thus only one SIGMA per CV has to be specified for each umbrella.
You can choose the umbrellas manually, or place them on a grid, or along a path, similar to [PATH](PATH.md).
They must cover all the CV space that one wishes to sample.

The first column of the umbrellas file is always ignored and must be called "time".
You can also use as input file a STATE file from an earlier [OPES_METAD](OPES_METAD.md) run (or an [OPES_METAD_EXPLORE](OPES_METAD_EXPLORE.md) run, if you combine it with other ECVs).

Similarly to [ECV_UMBRELLAS_LINE](ECV_UMBRELLAS_LINE.md), you should set the flag ADD_P0 if you think your umbrellas might not properly cover all the CV region relevant for the unbiased distribution.
You can also use BARRIER to set the maximum barrier height to be explored, and avoid huge biases at the beginning of your simulation.
See also Appendix B of the paper cited below for more details on these last two options.

## Examples

```plumed
#SETTINGS INPUTFILES=extras/Umbrellas.data

cv1: DISTANCE ATOMS=1,2
cv2: DISTANCE ATOMS=3,4
cv3: DISTANCE ATOMS=4,1
ecv: ECV_UMBRELLAS_FILE ...
   ARG=cv1,cv2,cv3
   FILE=extras/Umbrellas.data
   ADD_P0 BARRIER=70
...
opes: OPES_EXPANDED ARG=ecv.* PACE=500
PRINT FILE=COLVAR STRIDE=500 ARG=cv1,cv2,cv3,opes.bias
```

*/
//+ENDPLUMEDOC

class ECVumbrellasFile :
  public ExpansionCVs {
private:
  double barrier_;
  unsigned P0_contribution_;
  bool lower_only_;

  std::vector< std::vector<double> > centers_;
  std::vector< std::vector<double> > sigmas_;

  std::vector< std::vector<double> > ECVs_;
  std::vector< std::vector<double> > derECVs_;
  void initECVs();

public:
  explicit ECVumbrellasFile(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  void calculateECVs(const double *) override;
  const double * getPntrToECVs(unsigned) override;
  const double * getPntrToDerECVs(unsigned) override;
  std::vector<std::string> getLambdas() const override;
  void initECVs_observ(const std::vector<double>&,const unsigned,const unsigned) override;
  void initECVs_restart(const std::vector<std::string>&) override;
};

PLUMED_REGISTER_ACTION(ECVumbrellasFile,"ECV_UMBRELLAS_FILE")

void ECVumbrellasFile::registerKeywords(Keywords& keys) {
  ExpansionCVs::registerKeywords(keys);
  keys.add("compulsory","FILE","the name of the file containing the umbrellas");
  keys.add("optional","BARRIER","a guess of the free energy barrier to be overcome (better to stay higher than lower)");
  keys.addFlag("ADD_P0",false,"add the unbiased Boltzmann distribution to the target distribution, to make sure to sample it");
  keys.addFlag("LOWER_HALF_ONLY",false,"use only the lower half of each umbrella potentials");
  keys.addDOI("10.1103/PhysRevX.10.041034");
}

ECVumbrellasFile::ECVumbrellasFile(const ActionOptions&ao):
  Action(ao),
  ExpansionCVs(ao) {
//get number of CVs
  const unsigned ncv=getNumberOfArguments();
  centers_.resize(ncv);
  sigmas_.resize(ncv);

//set P0_contribution_
  bool add_P0=false;
  parseFlag("ADD_P0",add_P0);
  if(add_P0) {
    P0_contribution_=1;
  } else {
    P0_contribution_=0;
  }

//set barrier_
  barrier_=std::numeric_limits<double>::infinity();
  parse("BARRIER",barrier_);
  parseFlag("LOWER_HALF_ONLY",lower_only_);

//set umbrellas
  std::string umbrellasFileName;
  parse("FILE",umbrellasFileName);
  IFile ifile;
  ifile.link(*this);
  if(ifile.FileExist(umbrellasFileName)) {
    log.printf("  reading from FILE '%s'\n",umbrellasFileName.c_str());
    ifile.open(umbrellasFileName);
    ifile.allowIgnoredFields();
    double time; //first field is ignored
    while(ifile.scanField("time",time)) {
      for(unsigned j=0; j<ncv; j++) {
        double centers_j;
        ifile.scanField(getPntrToArgument(j)->getName(),centers_j);
        centers_[j].push_back(centers_j); //this might be slow
      }
      for(unsigned j=0; j<ncv; j++) {
        double sigmas_j;
        ifile.scanField("sigma_"+getPntrToArgument(j)->getName(),sigmas_j);
        sigmas_[j].push_back(sigmas_j);
      }
      ifile.scanField();
    }
  } else {
    plumed_merror("Umbrellas FILE '"+umbrellasFileName+"' not found");
  }

  checkRead();

//extra consistency checks
  const unsigned sizeUmbrellas=centers_[0].size();
  for(unsigned j=0; j<ncv; j++) {
    plumed_massert(centers_[j].size()==sizeUmbrellas,"mismatch in the number of centers read from file");
    plumed_massert(sigmas_[j].size()==sizeUmbrellas,"mismatch in the number of sigmas read from file");
  }

//set ECVs stuff
  totNumECVs_=sizeUmbrellas+P0_contribution_;
  ECVs_.resize(ncv,std::vector<double>(totNumECVs_));
  derECVs_.resize(ncv,std::vector<double>(totNumECVs_));

//printing some info
  log.printf("  total number of umbrellas = %u\n",sizeUmbrellas);
  if(barrier_!=std::numeric_limits<double>::infinity()) {
    log.printf("  guess for free energy BARRIER = %g\n",barrier_);
  }
  if(P0_contribution_==1) {
    log.printf(" -- ADD_P0: the target includes also the unbiased probability itself\n");
  }
  if(lower_only_) {
    log.printf(" -- LOWER_HALF_ONLY: the ECVs are set to zero for values of the CV above the respective center\n");
  }
}

void ECVumbrellasFile::calculateECVs(const double * cv) {
  if(lower_only_) {
    for(unsigned j=0; j<getNumberOfArguments(); j++) {
      for(unsigned k=P0_contribution_; k<totNumECVs_; k++) { //if ADD_P0, the first ECVs=0
        const unsigned kk=k-P0_contribution_;
        const double dist_jk=difference(j,centers_[j][kk],cv[j])/sigmas_[j][kk]; //PBC might be present
        if(dist_jk>=0) {
          ECVs_[j][k]=0;
          derECVs_[j][k]=0;
        } else {
          ECVs_[j][k]=0.5*std::pow(dist_jk,2);
          derECVs_[j][k]=dist_jk/sigmas_[j][kk];
        }
      }
    }
  } else {
    for(unsigned j=0; j<getNumberOfArguments(); j++) {
      for(unsigned k=P0_contribution_; k<totNumECVs_; k++) { //if ADD_P0, the first ECVs=0
        const unsigned kk=k-P0_contribution_;
        const double dist_jk=difference(j,centers_[j][kk],cv[j])/sigmas_[j][kk]; //PBC might be present
        ECVs_[j][k]=0.5*std::pow(dist_jk,2);
        derECVs_[j][k]=dist_jk/sigmas_[j][kk];
      }
    }
  }
}

const double * ECVumbrellasFile::getPntrToECVs(unsigned j) {
  plumed_massert(isReady_,"cannot access ECVs before initialization");
  plumed_massert(j<getNumberOfArguments(),getName()+" has fewer CVs");
  return &ECVs_[j][0];
}

const double * ECVumbrellasFile::getPntrToDerECVs(unsigned j) {
  plumed_massert(isReady_,"cannot access ECVs before initialization");
  plumed_massert(j<getNumberOfArguments(),getName()+" has fewer CVs");
  return &derECVs_[j][0];
}

std::vector<std::string> ECVumbrellasFile::getLambdas() const {
  //notice that sigmas are not considered!
  std::vector<std::string> lambdas(totNumECVs_);
  if(P0_contribution_==1) {
    std::ostringstream subs;
    subs<<"P0";
    for(unsigned j=1; j<getNumberOfArguments(); j++) {
      subs<<"_P0";
    }
    lambdas[0]=subs.str();
  }
  for(unsigned k=P0_contribution_; k<totNumECVs_; k++) {
    const unsigned kk=k-P0_contribution_;
    std::ostringstream subs;
    subs<<centers_[0][kk];
    for(unsigned j=1; j<getNumberOfArguments(); j++) {
      subs<<"_"<<centers_[j][kk];
    }
    lambdas[k]=subs.str();
  }
  return lambdas;
}

void ECVumbrellasFile::initECVs() {
  plumed_massert(!isReady_,"initialization should not be called twice");
  isReady_=true;
  log.printf("  *%4u windows for %s\n",totNumECVs_,getName().c_str());
}

void ECVumbrellasFile::initECVs_observ(const std::vector<double>& all_obs_cvs,const unsigned ncv,const unsigned index_j) {
  //this non-linear exansion never uses automatic initialization
  initECVs();
  calculateECVs(&all_obs_cvs[index_j]); //use only first obs point
  for(unsigned j=0; j<getNumberOfArguments(); j++)
    for(unsigned k=P0_contribution_; k<totNumECVs_; k++) {
      ECVs_[j][k]=std::min(barrier_/kbt_,ECVs_[j][k]);
    }
}

void ECVumbrellasFile::initECVs_restart(const std::vector<std::string>& lambdas) {
  std::size_t pos=0;
  for(unsigned j=0; j<getNumberOfArguments()-1; j++) {
    pos=lambdas[0].find("_",pos+1);  //checking only lambdas[0] is hopefully enough
  }
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
