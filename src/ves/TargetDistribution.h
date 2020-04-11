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
#ifndef __PLUMED_ves_TargetDistribution_h
#define __PLUMED_ves_TargetDistribution_h

#include "core/Action.h"

#include <vector>
#include <string>
#include <cmath>

#define PLUMED_VES_TARGETDISTRIBUTION_INIT(ao) TargetDistribution(ao)

namespace PLMD {

/**
\ingroup INHERIT
Abstract base class for implenting new target distributions.
*/

class Action;
class Grid;
class Value;
class Keywords;

namespace ves {

class TargetDistModifer;
class VesBias;

class TargetDistribution :
  public Action
{
private:
  enum TargetDistType {
    static_targetdist,
    dynamic_targetdist
  } type_;
  //
  bool force_normalization_;
  bool check_normalization_;
  bool check_nonnegative_;
  bool check_nan_inf_;
  bool shift_targetdist_to_zero_;
  // dimension of the distribution
  unsigned int dimension_;
  // grid parameters
  std::vector<Value*> grid_args_;
  //
  Grid* targetdist_grid_pntr_;
  Grid* log_targetdist_grid_pntr_;
  //
  std::vector<TargetDistModifer*> targetdist_modifer_pntrs_;
  //
  Action* action_pntr_;
  VesBias* vesbias_pntr_;
  //
  bool needs_bias_grid_;
  bool needs_bias_withoutcutoff_grid_;
  bool needs_fes_grid_;
  //
  Grid* bias_grid_pntr_;
  Grid* bias_withoutcutoff_grid_pntr_;
  Grid* fes_grid_pntr_;
  //
  bool static_grid_calculated;
  //
  bool allow_bias_cutoff_;
  bool bias_cutoff_active_;
  //
  void calculateStaticDistributionGrid();
  void updateBiasCutoffForTargetDistGrid();
  void checkNanAndInf();
protected:
  void setStatic() {type_=static_targetdist;}
  void setDynamic() {type_=dynamic_targetdist;}
  // set the that target distribution is normalized
  void setForcedNormalization() {force_normalization_=true; check_normalization_=false;}
  void unsetForcedNormalization() {force_normalization_=false; check_normalization_=true;};
  //
  void setBiasGridNeeded() {needs_bias_grid_=true;}
  void setBiasWithoutCutoffGridNeeded() {needs_bias_withoutcutoff_grid_=true;}
  void setFesGridNeeded() {needs_fes_grid_=true;}
  //
  VesBias* getPntrToVesBias() const;
  Action* getPntrToAction() const;
  //
  virtual void setupAdditionalGrids(const std::vector<Value*>&, const std::vector<std::string>&, const std::vector<std::string>&, const std::vector<unsigned int>&) {}
  //
  void normalizeTargetDistGrid();
  //
  Grid& targetDistGrid() const {return *targetdist_grid_pntr_;}
  Grid& logTargetDistGrid() const {return *log_targetdist_grid_pntr_;}
  //
  Grid* getBiasGridPntr() const {return bias_grid_pntr_;}
  Grid* getBiasWithoutCutoffGridPntr() const {return bias_withoutcutoff_grid_pntr_;}
  Grid* getFesGridPntr() const {return fes_grid_pntr_;}
  //
  double getBeta() const;
  //
  void applyTargetDistModiferToGrid(TargetDistModifer* modifer_pntr);
  //
  void setMinimumOfTargetDistGridToZero();
  void updateLogTargetDistGrid();
  //
  virtual void updateGrid() {calculateStaticDistributionGrid();}
public:
  static void registerKeywords(Keywords&);
  explicit TargetDistribution(const ActionOptions&);
  ~TargetDistribution();
  //
  bool isStatic() const {return type_==static_targetdist;}
  bool isDynamic() const {return type_==dynamic_targetdist;}
  // is the target distribution normalize or not
  bool forcedNormalization() const {return force_normalization_;};
  bool isTargetDistGridShiftedToZero() const {return shift_targetdist_to_zero_;}
  //
  bool biasGridNeeded() const {return needs_bias_grid_;}
  bool biasWithoutCutoffGridNeeded() const {return needs_bias_withoutcutoff_grid_;}
  bool fesGridNeeded()  const {return needs_fes_grid_;}
  //
  void allowBiasCutoff() {allow_bias_cutoff_=true;}
  void doNotAllowBiasCutoff() {allow_bias_cutoff_=false;}
  bool isBiasCutoffAllowed() const {return allow_bias_cutoff_;}
  bool biasCutoffActive() const {return bias_cutoff_active_;}
  //
  void setDimension(const unsigned int dimension);
  unsigned getDimension() const {return dimension_;}
  //
  virtual void linkVesBias(VesBias*);
  virtual void linkAction(Action*);
  //
  virtual void linkBiasGrid(Grid*);
  virtual void linkBiasWithoutCutoffGrid(Grid*);
  virtual void linkFesGrid(Grid*);
  //
  void setupBiasCutoff();
  //
  Grid* getTargetDistGridPntr() const {return targetdist_grid_pntr_;}
  Grid* getLogTargetDistGridPntr() const {return log_targetdist_grid_pntr_;}
  //
  void clearLogTargetDistGrid();
  // calculate the target distribution itself
  virtual double getValue(const std::vector<double>&) const = 0;
  //
  void setupGrids(const std::vector<Value*>&, const std::vector<std::string>&, const std::vector<std::string>&, const std::vector<unsigned int>&);
  //
  Grid getMarginal(const std::vector<std::string>&);
  //
  void updateTargetDist();
  //
  void readInRestartTargetDistGrid(const std::string&);
  //
  static double integrateGrid(const Grid*);
  static double normalizeGrid(Grid*);
  static Grid getMarginalDistributionGrid(Grid*, const std::vector<std::string>&);
  // empty standard action stuff
  void update() override {};
  void apply() override {};
  void calculate() override {};
};


inline
VesBias* TargetDistribution::getPntrToVesBias() const {
  plumed_massert(vesbias_pntr_!=NULL,"the VES bias has not been linked");
  return vesbias_pntr_;
}


inline
Action* TargetDistribution::getPntrToAction() const {
  plumed_massert(action_pntr_!=NULL,"the action has not been linked");
  return action_pntr_;
}


inline
void TargetDistribution::normalizeTargetDistGrid() {
  double normalization = normalizeGrid(targetdist_grid_pntr_);
  if(normalization<0.0) {plumed_merror(getName()+": something went wrong trying to normalize the target distribution, integrating over it gives a negative value.");}
}




}
}
#endif
