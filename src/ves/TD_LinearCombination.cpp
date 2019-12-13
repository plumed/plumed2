/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The VES code team
   (see the PEOPLE-VES file at the root of this folder for a list of names)

   See http://www.ves-code.org for more information.

   This file is part of VES code module.

   The VES code module is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The VES code module is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the VES code module.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include "TargetDistribution.h"
#include "VesTools.h"

#include "core/ActionRegister.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"
#include "tools/Grid.h"


namespace PLMD {

// class Grid;
class Action;

namespace ves {

//+PLUMEDOC VES_TARGETDIST TD_LINEAR_COMBINATION
/*
Target distribution given by linear combination of distributions (static or dynamic).

Employ a target distribution that is a linear combination of the other
distributions, defined as
\f[
p(\mathbf{s}) = \sum_{i} w_{i} \, p_{i}(\mathbf{s})
\f]
where the weights \f$w_{i}\f$ are normalized to 1, \f$\sum_{i}w_{i}=1\f$.

The labels of the distributions \f$p_{i}(\mathbf{s})\f$ to be used in the
linear combination are given in the DISTRIBUTIONS keyword.

The weights \f$w_{i}\f$ can be given using
the WEIGHTS keyword. The distributions are weighted equally if no weights are given.

It is assumed that all the distributions \f$p_{i}(\mathbf{s})\f$ are normalized.
If that is not the case for some reason should you
normalize each distribution separately by using the NORMALIZE
keyword when defining them in the input file (i.e. before the
TD_LINEAR_COMBINATION action).
Note that normalizing the overall
linear combination will generally lead to different results than normalizing
each distribution separately.

The linear combination will be a dynamic target distribution if one or more
of the distributions used is a dynamic distribution, otherwise it will be a
static distribution.

\par Examples

Here we employ a linear combination of a uniform and a Gaussian distribution.
No weights are given so the two distributions will be weighted equally.
\plumedfile
td_uni: TD_UNIFORM

td_gauss: TD_GAUSSIAN CENTER1=-2.0 SIGMA1=0.5

td_comb: TD_LINEAR_COMBINATION DISTRIBUTIONS=td_uni,td_gauss
\endplumedfile

Here we employ a linear combination of a uniform and two Gaussian distribution.
The weights are automatically normalized to 1 such that giving
WEIGHTS=1.0,1.0,2.0 as we do here is equal to giving WEIGHTS=0.25,0.25,0.50.
\plumedfile
td_uni: TD_UNIFORM

td_gauss1: TD_GAUSSIAN CENTER1=-2.0,-2.0 SIGMA1=0.5,0.3

td_gauss2: TD_GAUSSIAN CENTER1=+2.0,+2.0 SIGMA1=0.3,0.5

TD_LINEAR_COMBINATION ...
 DISTRIBUTIONS=td_uni,td_gauss1,td_gauss2
 WEIGHTS=1.0,1.0,2.0
 LABEL=td_comb
... TD_LINEAR_COMBINATION
\endplumedfile

In the above example the two Gaussian kernels are given using two separate
DISTRIBUTION keywords. As the \ref TD_GAUSSIAN target distribution allows multiple
centers is it also possible to use just one DISTRIBUTION keyword for the two
Gaussian kernels. This is shown in the following example which will give the
exact same result as the one above as the weights have been appropriately
adjusted
\plumedfile
td_uni: TD_UNIFORM

TD_GAUSSIAN ...
 CENTER1=-2.0,-2.0  SIGMA1=0.5,0.3
 CENTER2=+2.0,+2.0  SIGMA2=0.3,0.5
 WEIGHTS=1.0,2.0
 LABEL=td_gauss
... TD_GAUSSIAN

TD_LINEAR_COMBINATION ...
 DISTRIBUTIONS=td_uni,td_gauss
 WEIGHTS=0.25,0.75
 LABEL=td_comb
... TD_LINEAR_COMBINATION
\endplumedfile

*/
//+ENDPLUMEDOC

class VesBias;

class TD_LinearCombination: public TargetDistribution {
private:
  std::vector<TargetDistribution*> distribution_pntrs_;
  std::vector<Grid*> grid_pntrs_;
  std::vector<double> weights_;
  unsigned int ndist_;
  void setupAdditionalGrids(const std::vector<Value*>&, const std::vector<std::string>&, const std::vector<std::string>&, const std::vector<unsigned int>&) override;
public:
  static void registerKeywords(Keywords&);
  explicit TD_LinearCombination(const ActionOptions& ao);
  void updateGrid() override;
  double getValue(const std::vector<double>&) const override;
  //
  void linkVesBias(VesBias*) override;
  void linkAction(Action*) override;
  //
  void linkBiasGrid(Grid*) override;
  void linkBiasWithoutCutoffGrid(Grid*) override;
  void linkFesGrid(Grid*) override;
  //
};


PLUMED_REGISTER_ACTION(TD_LinearCombination,"TD_LINEAR_COMBINATION")


void TD_LinearCombination::registerKeywords(Keywords& keys) {
  TargetDistribution::registerKeywords(keys);
  keys.add("compulsory","DISTRIBUTIONS","The labels of the target distribution actions to be used in the linear combination.");
  keys.add("optional","WEIGHTS","The weights of target distributions. Have to be as many as the number of target distribution labels given in DISTRIBUTIONS. If no weights are given the distributions are weighted equally. The weights are automatically normalized to 1.");
  keys.use("WELLTEMPERED_FACTOR");
  //keys.use("SHIFT_TO_ZERO");
  keys.use("NORMALIZE");
}


TD_LinearCombination::TD_LinearCombination(const ActionOptions& ao):
  PLUMED_VES_TARGETDISTRIBUTION_INIT(ao),
  distribution_pntrs_(0),
  grid_pntrs_(0),
  weights_(0),
  ndist_(0)
{
  std::vector<std::string> targetdist_labels;
  parseVector("DISTRIBUTIONS",targetdist_labels);

  std::string error_msg = "";
  distribution_pntrs_ = VesTools::getPointersFromLabels<TargetDistribution*>(targetdist_labels,plumed.getActionSet(),error_msg);
  if(error_msg.size()>0) {plumed_merror("Error in keyword DISTRIBUTIONS of "+getName()+": "+error_msg);}

  for(unsigned int i=0; i<distribution_pntrs_.size(); i++) {
    if(distribution_pntrs_[i]->isDynamic()) {setDynamic();}
    if(distribution_pntrs_[i]->fesGridNeeded()) {setFesGridNeeded();}
    if(distribution_pntrs_[i]->biasGridNeeded()) {setBiasGridNeeded();}
  }

  ndist_ = distribution_pntrs_.size();
  grid_pntrs_.assign(ndist_,NULL);
  if(ndist_==0) {plumed_merror(getName()+ ": no distributions are given.");}
  if(ndist_==1) {plumed_merror(getName()+ ": giving only one distribution does not make sense.");}
  //
  parseVector("WEIGHTS",weights_);
  if(weights_.size()==0) {weights_.assign(distribution_pntrs_.size(),1.0);}
  if(distribution_pntrs_.size()!=weights_.size()) {
    plumed_merror(getName()+ ": there has to be as many weights given in WEIGHTS as the number of target distribution labels given in DISTRIBUTIONS");
  }
  //
  double sum_weights=0.0;
  for(unsigned int i=0; i<weights_.size(); i++) {sum_weights+=weights_[i];}
  for(unsigned int i=0; i<weights_.size(); i++) {weights_[i]/=sum_weights;}
  checkRead();
}


double TD_LinearCombination::getValue(const std::vector<double>& argument) const {
  plumed_merror("getValue not implemented for TD_LinearCombination");
  return 0.0;
}


void TD_LinearCombination::setupAdditionalGrids(const std::vector<Value*>& arguments, const std::vector<std::string>& min, const std::vector<std::string>& max, const std::vector<unsigned int>& nbins) {
  for(unsigned int i=0; i<ndist_; i++) {
    distribution_pntrs_[i]->setupGrids(arguments,min,max,nbins);
    if(distribution_pntrs_[i]->getDimension()!=this->getDimension()) {
      plumed_merror(getName() + ": all target distribution must have the same dimension");
    }
    grid_pntrs_[i]=distribution_pntrs_[i]->getTargetDistGridPntr();
  }
}


void TD_LinearCombination::updateGrid() {
  for(unsigned int i=0; i<ndist_; i++) {
    distribution_pntrs_[i]->updateTargetDist();
  }
  for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++) {
    double value = 0.0;
    for(unsigned int i=0; i<ndist_; i++) {
      value += weights_[i]*grid_pntrs_[i]->getValue(l);
    }
    targetDistGrid().setValue(l,value);
    logTargetDistGrid().setValue(l,-std::log(value));
  }
  logTargetDistGrid().setMinToZero();
}


void TD_LinearCombination::linkVesBias(VesBias* vesbias_pntr_in) {
  TargetDistribution::linkVesBias(vesbias_pntr_in);
  for(unsigned int i=0; i<ndist_; i++) {
    distribution_pntrs_[i]->linkVesBias(vesbias_pntr_in);
  }
}


void TD_LinearCombination::linkAction(Action* action_pntr_in) {
  TargetDistribution::linkAction(action_pntr_in);
  for(unsigned int i=0; i<ndist_; i++) {
    distribution_pntrs_[i]->linkAction(action_pntr_in);
  }
}


void TD_LinearCombination::linkBiasGrid(Grid* bias_grid_pntr_in) {
  TargetDistribution::linkBiasGrid(bias_grid_pntr_in);
  for(unsigned int i=0; i<ndist_; i++) {
    distribution_pntrs_[i]->linkBiasGrid(bias_grid_pntr_in);
  }
}


void TD_LinearCombination::linkBiasWithoutCutoffGrid(Grid* bias_withoutcutoff_grid_pntr_in) {
  TargetDistribution::linkBiasWithoutCutoffGrid(bias_withoutcutoff_grid_pntr_in);
  for(unsigned int i=0; i<ndist_; i++) {
    distribution_pntrs_[i]->linkBiasWithoutCutoffGrid(bias_withoutcutoff_grid_pntr_in);
  }
}


void TD_LinearCombination::linkFesGrid(Grid* fes_grid_pntr_in) {
  TargetDistribution::linkFesGrid(fes_grid_pntr_in);
  for(unsigned int i=0; i<ndist_; i++) {
    distribution_pntrs_[i]->linkFesGrid(fes_grid_pntr_in);
  }
}


}
}
