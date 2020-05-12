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
namespace ves {

//+PLUMEDOC VES_TARGETDIST TD_PRODUCT_DISTRIBUTION
/*
Target distribution given by a separable product of one-dimensional distributions (static or dynamic).

Employ a target distribution that is a separable product
of one-dimensional distributions, defined as
\f[
p(\mathbf{s}) =
\prod_{k}^{d} p_{k}(s_{k})
\f]
where \f$d\f$ is the number of arguments used and \f$p_{k}(s_{k})\f$ is the
one-dimensional distribution corresponding to the \f$k\f$-th argument.

Note the difference between this target distribution and the one defined in
\ref TD_PRODUCT_COMBINATION. Here we have a separable distribution given as a
product of one-dimensional distribution \f$p_{k}(s_{k})\f$.

The labels of the one-dimensional distributions \f$p_{k}(s_{k})\f$ to be
used in the product distribution are given in the DISTRIBUTIONS keyword.
Note that the order of the labels is very important.

It is assumed that all the distributions to be used in the product distribution
are normalized. If that is not the case you need to
normalize the distributions by using the NORMALIZE keyword.
Here it does not matter if you normalize each distribution separately
or the overall product, it will give the same results.

The product distribution will be a dynamic target distribution if one or more
of the distributions used is a dynamic distribution. Otherwise it will be a
static distribution.

\par Examples

In the following example we employ a uniform distribution for
argument 1 and a Gaussian distribution for argument 2.
\plumedfile
target_uniform: TD_UNIFORM

target_Gaussian: TD_GAUSSIAN CENTER1=-2.0 SIGMA1=0.5

td_pd: TD_PRODUCT_DISTRIBUTION DISTRIBUTIONS=target_uniform,target_Gaussian
\endplumedfile
Note that order of the labels is important, using DISTRIBUTIONS=target_Gaussian,target_uniform
would mean that we would employ a Gaussian distribution for argument 1 and a uniform
distribution for argument 2, which would lead to completely different results.

*/
//+ENDPLUMEDOC

class TD_ProductDistribution: public TargetDistribution {
private:
  std::vector<TargetDistribution*> distribution_pntrs_;
  std::vector<Grid*> grid_pntrs_;
  unsigned int ndist_;
  void setupAdditionalGrids(const std::vector<Value*>&, const std::vector<std::string>&, const std::vector<std::string>&, const std::vector<unsigned int>&) override;
public:
  static void registerKeywords(Keywords&);
  explicit TD_ProductDistribution(const ActionOptions& ao);
  void updateGrid() override;
  double getValue(const std::vector<double>&) const override;
  //
  void linkVesBias(VesBias*) override;
  void linkAction(Action*) override;
  void linkBiasGrid(Grid*) override;
  void linkBiasWithoutCutoffGrid(Grid*) override;
  void linkFesGrid(Grid*) override;
};


PLUMED_REGISTER_ACTION(TD_ProductDistribution,"TD_PRODUCT_DISTRIBUTION")


void TD_ProductDistribution::registerKeywords(Keywords& keys) {
  TargetDistribution::registerKeywords(keys);
  keys.add("compulsory","DISTRIBUTIONS","Labels of the one-dimensional target distribution actions for each argument to be used in the product distribution. Note that order of the labels is important.");
  keys.use("WELLTEMPERED_FACTOR");
  keys.use("SHIFT_TO_ZERO");
  keys.use("NORMALIZE");
}


TD_ProductDistribution::TD_ProductDistribution(const ActionOptions& ao):
  PLUMED_VES_TARGETDISTRIBUTION_INIT(ao),
  distribution_pntrs_(0),
  grid_pntrs_(0),
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
  setDimension(ndist_);

  checkRead();
}


double TD_ProductDistribution::getValue(const std::vector<double>& argument) const {
  plumed_merror("getValue not implemented for TD_ProductDistribution");
  return 0.0;
}


void TD_ProductDistribution::setupAdditionalGrids(const std::vector<Value*>& arguments, const std::vector<std::string>& min, const std::vector<std::string>& max, const std::vector<unsigned int>& nbins) {
  for(unsigned int i=0; i<ndist_; i++) {
    std::vector<Value*> arg1d(1);
    std::vector<std::string> min1d(1);
    std::vector<std::string> max1d(1);
    std::vector<unsigned int> nbins1d(1);
    arg1d[0]=arguments[i];
    min1d[0]=min[i];
    max1d[0]=max[i];
    nbins1d[0]=nbins[i];
    distribution_pntrs_[i]->setupGrids(arg1d,min1d,max1d,nbins1d);
    grid_pntrs_[i]=distribution_pntrs_[i]->getTargetDistGridPntr();
    if(distribution_pntrs_[i]->getDimension()!=1 || grid_pntrs_[i]->getDimension()!=1) {
      plumed_merror(getName() + ": all target distributions must be one dimensional");
    }
  }
}


void TD_ProductDistribution::updateGrid() {
  for(unsigned int i=0; i<ndist_; i++) {
    distribution_pntrs_[i]->updateTargetDist();
  }
  for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++) {
    std::vector<unsigned int> indices = targetDistGrid().getIndices(l);
    double value = 1.0;
    for(unsigned int i=0; i<ndist_; i++) {
      value *= grid_pntrs_[i]->getValue(indices[i]);
    }
    targetDistGrid().setValue(l,value);
    logTargetDistGrid().setValue(l,-std::log(value));
  }
  logTargetDistGrid().setMinToZero();
}


void TD_ProductDistribution::linkVesBias(VesBias* vesbias_pntr_in) {
  TargetDistribution::linkVesBias(vesbias_pntr_in);
  for(unsigned int i=0; i<ndist_; i++) {
    distribution_pntrs_[i]->linkVesBias(vesbias_pntr_in);
  }
}


void TD_ProductDistribution::linkAction(Action* action_pntr_in) {
  TargetDistribution::linkAction(action_pntr_in);
  for(unsigned int i=0; i<ndist_; i++) {
    distribution_pntrs_[i]->linkAction(action_pntr_in);
  }
}


void TD_ProductDistribution::linkBiasGrid(Grid* bias_grid_pntr_in) {
  TargetDistribution::linkBiasGrid(bias_grid_pntr_in);
  for(unsigned int i=0; i<ndist_; i++) {
    distribution_pntrs_[i]->linkBiasGrid(bias_grid_pntr_in);
  }
}


void TD_ProductDistribution::linkBiasWithoutCutoffGrid(Grid* bias_withoutcutoff_grid_pntr_in) {
  TargetDistribution::linkBiasWithoutCutoffGrid(bias_withoutcutoff_grid_pntr_in);
  for(unsigned int i=0; i<ndist_; i++) {
    distribution_pntrs_[i]->linkBiasWithoutCutoffGrid(bias_withoutcutoff_grid_pntr_in);
  }
}


void TD_ProductDistribution::linkFesGrid(Grid* fes_grid_pntr_in) {
  TargetDistribution::linkFesGrid(fes_grid_pntr_in);
  for(unsigned int i=0; i<ndist_; i++) {
    distribution_pntrs_[i]->linkFesGrid(fes_grid_pntr_in);
  }
}


}
}
