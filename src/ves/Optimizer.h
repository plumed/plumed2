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
#ifndef __PLUMED_ves_Optimizer_h
#define __PLUMED_ves_Optimizer_h

#include "VesBias.h"

#include "core/ActionPilot.h"
#include "core/ActionWithValue.h"

#include <vector>
#include <string>
#include <cmath>


#define PLUMED_VES_OPTIMIZER_INIT(ao) Action(ao),Optimizer(ao)

namespace PLMD {

/**
\ingroup INHERIT
Abstract base class for implenting new optimization methods
*/

class OFile;

namespace ves {

class CoeffsVector;
class VesBias;


class Optimizer :
  public ActionPilot,
  public ActionWithValue
{
private:
  const std::string description_;
  const std::string type_;
  //
  std::vector<double> stepsizes_;
  std::vector<double> current_stepsizes;
  bool fixed_stepsize_;
  //
  unsigned int iter_counter;
  //
  bool use_hessian_;
  bool diagonal_hessian_;
  bool hessian_covariance_from_averages_;
  //
  bool monitor_instantaneous_gradient_;
  //
  bool use_mwalkers_mpi_;
  bool mwalkers_mpi_single_files_;
  //
  std::vector<bool> dynamic_targetdists_;
  unsigned int ustride_targetdist_;
  //
  unsigned int ustride_reweightfactor_;
  //
  std::string coeffssetid_prefix_;
  //
  unsigned int coeffs_wstride_;
  std::vector<OFile*> coeffsOFiles_;
  std::string coeffs_output_fmt_;
  //
  unsigned int gradient_wstride_;
  std::vector<OFile*> gradientOFiles_;
  std::string gradient_output_fmt_;
  //
  unsigned int hessian_wstride_;
  std::vector<OFile*> hessianOFiles_;
  std::string hessian_output_fmt_;
  //
  unsigned int targetdist_averages_wstride_;
  std::vector<OFile*> targetdist_averagesOFiles_;
  std::string targetdist_averages_output_fmt_;
  //
  unsigned int nbiases_;
  std::vector<VesBias*> bias_pntrs_;
  //
  unsigned int ncoeffssets_;
  std::vector<CoeffsVector*> coeffs_pntrs_;
  std::vector<CoeffsVector*> aux_coeffs_pntrs_;
  std::vector<CoeffsVector*> gradient_pntrs_;
  std::vector<CoeffsVector*> aver_gradient_pntrs_;
  std::vector<CoeffsMatrix*> hessian_pntrs_;
  std::vector<CoeffsVector*> coeffs_mask_pntrs_;
  std::vector<CoeffsVector*> targetdist_averages_pntrs_;
  //
  bool identical_coeffs_shape_;
  //
  bool bias_output_active_;
  unsigned int bias_output_stride_;
  bool fes_output_active_;
  unsigned int fes_output_stride_;
  bool fesproj_output_active_;
  unsigned int fesproj_output_stride_;
  bool targetdist_output_active_;
  unsigned int targetdist_output_stride_;
  bool targetdist_proj_output_active_;
  unsigned int targetdist_proj_output_stride_;
  //
  bool isFirstStep;
  //
private:
  void updateOutputComponents();
  void writeOutputFiles(const unsigned int coeffs_id = 0);
  void readCoeffsFromFiles(const std::vector<std::string>&, const bool);
  void setAllCoeffsSetIterationCounters();
protected:
  void turnOnHessian();
  void turnOffHessian();
  std::vector<CoeffsMatrix*> enableHessian(VesBias*, const bool diagonal_hessian=false);
  // CoeffsMatrix* switchToDiagonalHessian(VesBias*);
  // CoeffsMatrix* switchToFullHessian(VesBias*);
  //
  CoeffsVector& Coeffs(const unsigned int coeffs_id = 0) const;
  CoeffsVector& AuxCoeffs(const unsigned int coeffs_id = 0) const;
  CoeffsVector& Gradient(const unsigned int coeffs_id = 0) const;
  CoeffsMatrix& Hessian(const unsigned int coeffs_id = 0) const;
  CoeffsVector& CoeffsMask(const unsigned int coeffs_id = 0) const;
  CoeffsVector& TargetDistAverages(const unsigned int coeffs_id = 0) const;
  double StepSize(const unsigned int coeffs_id = 0) const;
  virtual void coeffsUpdate(const unsigned int coeffs_id = 0) = 0;
  void setCurrentStepSize(const double,const unsigned int i = 0);
  void setCurrentStepSizes(const std::vector<double>&);
  //
  void turnOffCoeffsOutputFiles();
  //
  template<class T>
  bool parseMultipleValues(const std::string&, std::vector<T>&);
  template<class T>
  bool parseMultipleValues(const std::string&, std::vector<T>&, const T&);
  void parseFilenames(const std::string&, std::vector<std::string>&, const std::string&);
  void parseFilenames(const std::string&, std::vector<std::string>&);
  void addCoeffsSetIDsToFilenames(std::vector<std::string>&, std::string&);
  void setupOFiles(std::vector<std::string>&, std::vector<OFile*>&, const bool multi_sim_single_files=false);
public:
  static void registerKeywords(Keywords&);
  static void useMultipleWalkersKeywords(Keywords&);
  static void useHessianKeywords(Keywords&);
  static void useFixedStepSizeKeywords(Keywords&);
  static void useDynamicStepSizeKeywords(Keywords&);
  static void useMaskKeywords(Keywords&);
  static void useRestartKeywords(Keywords&);
  static void useMonitorAverageGradientKeywords(Keywords&);
  static void useDynamicTargetDistributionKeywords(Keywords&);
  static void useReweightFactorKeywords(Keywords&);
  //
  explicit Optimizer(const ActionOptions&ao);
  ~Optimizer();
  std::string getType() const {return type_;}
  std::string getDescription() const {return description_;}
  //
  unsigned int numberOfBiases() const {return nbiases_;}
  unsigned int numberOfCoeffsSets() const {return ncoeffssets_;}
  //
  std::vector<double> getStepSizes() const;
  std::vector<double> getCurrentStepSizes() const;
  double getStepSize(const unsigned int coeffs_id = 0) const;
  double getCurrentStepSize(const unsigned int coeffs_id = 0) const;
  void setStepSizes(const std::vector<double>&);
  void setStepSize(const double, const unsigned int coeffs_id = 0);
  //
  unsigned int getIterationCounter() const;
  double getIterationCounterDbl() const;
  std::string getIterationCounterStr(const int offset=0) const;
  void setIterationCounter(const unsigned int);
  void increaseIterationCounter();
  //
  void apply() override {};
  void calculate() override {};
  void update() override;
  unsigned int getNumberOfDerivatives() override {return 0;}
  //
  bool fixedStepSize() const {return fixed_stepsize_;}
  bool dynamicStepSize() const {return !fixed_stepsize_;}
  //
  bool useHessian() const {return use_hessian_;}
  bool diagonalHessian() const {return diagonal_hessian_;}
  //
  bool useMultipleWalkers() const {return use_mwalkers_mpi_;}
  //
  std::vector<VesBias*> getBiasPntrs() const {return bias_pntrs_;}
  std::vector<CoeffsVector*> getCoeffsPntrs() const {return coeffs_pntrs_;}
  std::vector<CoeffsVector*> getAuxCoeffsPntrs() const {return aux_coeffs_pntrs_;}
  std::vector<CoeffsVector*> getGradientPntrs()const {return gradient_pntrs_;}
  std::vector<CoeffsMatrix*> getHessianPntrs() const {return hessian_pntrs_;}
  std::vector<CoeffsVector*> getCoeffsMaskPntrs() const {return coeffs_mask_pntrs_;}
  std::vector<CoeffsVector*> getTargetDistAveragesPntrs() const {return targetdist_averages_pntrs_;}
  //
  bool isBiasOutputActive() const {return bias_output_active_;}
  unsigned int getBiasOutputStride() const {return bias_output_stride_;}
  void setBiasOutputStride(unsigned int stride) {bias_output_stride_=stride;}
  void writeBiasOutputFiles() const;
  //
  bool isFesOutputActive() const {return fes_output_active_;}
  unsigned int getFesOutputStride() const {return fes_output_stride_;}
  void setFesOutputStride(unsigned int stride) {fes_output_stride_=stride;}
  void writeFesOutputFiles() const;
  //
  bool isFesProjOutputActive() const {return fesproj_output_active_;}
  unsigned int getFesProjOutputStride() const {return fesproj_output_stride_;}
  void setFesProjOutputStride(unsigned int stride) {fesproj_output_stride_=stride;}
  void writeFesProjOutputFiles() const;
  //
  bool isTargetDistOutputActive() const {return targetdist_output_active_;}
  unsigned int getTargetDistOutputStride() const {return targetdist_output_stride_;}
  void setTargetDistOutputStride(unsigned int stride) {targetdist_output_stride_=stride;}
  void writeTargetDistOutputFiles() const;
  //
  bool isTargetDistProjOutputActive() const {return targetdist_proj_output_active_;}
  unsigned int getTargetDistProjOutputStride() const {return targetdist_proj_output_stride_;}
  void setTargetDistProjOutputStride(unsigned int stride) {targetdist_proj_output_stride_=stride;}
  void writeTargetDistProjOutputFiles() const;
  //
};

inline
double Optimizer::StepSize(const unsigned int coeffs_id) const {return stepsizes_[coeffs_id];}

inline
CoeffsVector& Optimizer::Coeffs(const unsigned int coeffs_id) const {return *coeffs_pntrs_[coeffs_id];}

inline
CoeffsVector& Optimizer::AuxCoeffs(const unsigned int coeffs_id) const {return *aux_coeffs_pntrs_[coeffs_id];}

inline
CoeffsVector& Optimizer::Gradient(const unsigned int coeffs_id) const {return *gradient_pntrs_[coeffs_id];}

inline
CoeffsMatrix& Optimizer::Hessian(const unsigned int coeffs_id) const {
  plumed_massert(use_hessian_,"You cannot use the Hessian without asking for before");
  return *hessian_pntrs_[coeffs_id];
}

inline
CoeffsVector& Optimizer::CoeffsMask(const unsigned int coeffs_id) const {return *coeffs_mask_pntrs_[coeffs_id];}

inline
std::vector<double> Optimizer::getStepSizes() const {return stepsizes_;}

inline
std::vector<double> Optimizer::getCurrentStepSizes() const {return current_stepsizes;}

inline
double Optimizer::getStepSize(const unsigned int coeffs_id) const {return stepsizes_[coeffs_id];}

inline
double Optimizer::getCurrentStepSize(const unsigned int coeffs_id) const {return current_stepsizes[coeffs_id];}

inline
void Optimizer::setStepSizes(const std::vector<double>& stepsizes_in) {
  plumed_assert(stepsizes_in.size()==ncoeffssets_);
  stepsizes_ = stepsizes_in;
}

inline
void Optimizer::setStepSize(const double stepsize_in, const unsigned int coeffs_id) {
  stepsizes_[coeffs_id] = stepsize_in;
}

inline
void Optimizer::setCurrentStepSize(const double current_stepsize_in, const unsigned int coeffs_id) {
  current_stepsizes[coeffs_id] = current_stepsize_in;
}

inline
void Optimizer::setCurrentStepSizes(const std::vector<double>& current_stepsizes_in) {
  plumed_assert(current_stepsizes_in.size()==ncoeffssets_);
  current_stepsizes = current_stepsizes_in;
}

inline
unsigned int Optimizer::getIterationCounter() const {return iter_counter;}

inline
double Optimizer::getIterationCounterDbl() const {return static_cast<double>(iter_counter);}

inline
void Optimizer::increaseIterationCounter() {iter_counter++;}

inline
void Optimizer::setIterationCounter(const unsigned int iter_counter_in) {iter_counter = iter_counter_in;}


template<class T>
bool Optimizer::parseMultipleValues(const std::string& keyword, std::vector<T>& values) {
  plumed_assert(ncoeffssets_>0);
  plumed_assert(values.size()==0);
  bool identical_values=false;
  //
  parseVector(keyword,values);
  if(values.size()==1 && ncoeffssets_>1) {
    values.resize(ncoeffssets_,values[0]);
    identical_values=true;
  }
  if(values.size()>0 && values.size()!=ncoeffssets_) {
    std::string s1; Tools::convert(ncoeffssets_,s1);
    plumed_merror("Error in " + keyword + " keyword: either give 1 common value for all coefficient sets or " + s1 + " separate value for each set");
  }
  return identical_values;
}

template<class T>
bool Optimizer::parseMultipleValues(const std::string& keyword, std::vector<T>& values, const T& default_value) {
  bool identical_values = parseMultipleValues(keyword,values);
  if(values.size()==0) {
    values.resize(ncoeffssets_,default_value);
    identical_values=true;
  }
  return identical_values;
}

inline
void Optimizer::parseFilenames(const std::string& keyword, std::vector<std::string>& fnames, const std::string& default_fname) {
  if(parseMultipleValues<std::string>(keyword,fnames,default_fname)) {
    addCoeffsSetIDsToFilenames(fnames,coeffssetid_prefix_);
  }
}

inline
void Optimizer::parseFilenames(const std::string& keyword, std::vector<std::string>& fnames) {
  if(parseMultipleValues<std::string>(keyword,fnames)) {
    addCoeffsSetIDsToFilenames(fnames,coeffssetid_prefix_);
  }
}


}
}

#endif
