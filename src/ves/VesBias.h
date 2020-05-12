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
#ifndef __PLUMED_ves_VesBias_h
#define __PLUMED_ves_VesBias_h

#include "CoeffsVector.h"
#include "CoeffsMatrix.h"

#include "core/ActionPilot.h"
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "bias/Bias.h"

#include <vector>
#include <string>
#include <cmath>


#define PLUMED_VES_VESBIAS_INIT(ao) Action(ao),VesBias(ao)

namespace PLMD {

class Value;

namespace ves {

class CoeffsVector;
class CoeffsMatrix;
class BasisFunctions;
class Optimizer;
class TargetDistribution;
class FermiSwitchingFunction;

/**
\ingroup INHERIT
Abstract base class for implementing biases the extents the normal Bias.h class
to include functions related to the variational approach.
*/

class VesBias:
  public bias::Bias
{
private:
  unsigned int ncoeffssets_;
  std::vector<CoeffsVector*> coeffs_pntrs_;
  std::vector<CoeffsVector*> targetdist_averages_pntrs_;
  std::vector<CoeffsVector*> gradient_pntrs_;
  std::vector<CoeffsMatrix*> hessian_pntrs_;
  std::vector<std::vector<double> > sampled_averages;
  std::vector<std::vector<double> > sampled_cross_averages;
  bool use_multiple_coeffssets_;
  //
  std::vector<std::string> coeffs_fnames;
  //
  size_t ncoeffs_total_;
  //
  Optimizer* optimizer_pntr_;
  bool optimize_coeffs_;
  //
  bool compute_hessian_;
  bool diagonal_hessian_;
  //
  std::vector<unsigned int> aver_counters;
  //
  double kbt_;
  //
  std::vector<TargetDistribution*> targetdist_pntrs_;
  bool dynamic_targetdist_;
  //
  std::vector<unsigned int> grid_bins_;
  std::vector<double> grid_min_;
  std::vector<double> grid_max_;
  //
  std::string bias_filename_;
  std::string fes_filename_;
  std::string targetdist_filename_;
  std::string targetdist_averages_filename_;
  std::string coeffs_id_prefix_;
  //
  std::string bias_file_fmt_;
  std::string fes_file_fmt_;
  std::string targetdist_file_fmt_;
  std::string targetdist_restart_file_fmt_;
  //
  bool filenames_have_iteration_number_;
  //
  bool bias_fileoutput_active_;
  bool fes_fileoutput_active_;
  bool fesproj_fileoutput_active_;
  bool dynamic_targetdist_fileoutput_active_;
  bool static_targetdist_fileoutput_active_;
  //

  bool bias_cutoff_active_;
  double bias_cutoff_value_;
  double bias_current_max_value;
  FermiSwitchingFunction* bias_cutoff_swfunc_pntr_;
  //
  std::vector< std::vector<std::string> > projection_args_;
  //
  bool calc_reweightfactor_;
private:
  void initializeCoeffs(CoeffsVector*);
  std::vector<double> computeCovarianceFromAverages(const unsigned int) const;
  void multiSimSumAverages(const unsigned int, const double walker_weight=1.0);
protected:
  //
  void checkThatTemperatureIsGiven();
  //
  void addCoeffsSet(const std::vector<std::string>&,const std::vector<unsigned int>&);
  void addCoeffsSet(std::vector<Value*>&,std::vector<BasisFunctions*>&);
  void addCoeffsSet(CoeffsVector*);
  //
  std::string getCoeffsSetLabelString(const std::string&, const unsigned int coeffs_id = 0) const;
  void clearCoeffsPntrsVector() {coeffs_pntrs_.clear();}
  void addToSampledAverages(const std::vector<double>&, const unsigned int c_id = 0);
  void setTargetDistAverages(const std::vector<double>&, const unsigned int coeffs_id = 0);
  void setTargetDistAverages(const CoeffsVector&, const unsigned int coeffs_id= 0);
  void setTargetDistAveragesToZero(const unsigned int coeffs_id= 0);
  //
  bool readCoeffsFromFiles();
  //
  template<class T>
  bool parseMultipleValues(const std::string&, std::vector<T>&, unsigned int);
  template<class T>
  bool parseMultipleValues(const std::string&, std::vector<T>&, unsigned int, const T&);
  //
public:
  static void registerKeywords(Keywords&);
  explicit VesBias(const ActionOptions&ao);
  ~VesBias();
  //
  static void useInitialCoeffsKeywords(Keywords&);
  static void useTargetDistributionKeywords(Keywords&);
  static void useMultipleTargetDistributionKeywords(Keywords&);
  static void useGridBinKeywords(Keywords&);
  static void useGridLimitsKeywords(Keywords&);
  static void useBiasCutoffKeywords(Keywords&);
  static void useProjectionArgKeywords(Keywords&);
  static void useReweightFactorKeywords(Keywords&);
  //
  std::vector<CoeffsVector*> getCoeffsPntrs() const {return coeffs_pntrs_;}
  std::vector<CoeffsVector*> getTargetDistAveragesPntrs() const {return targetdist_averages_pntrs_;}
  std::vector<CoeffsVector*> getGradientPntrs()const {return gradient_pntrs_;}
  std::vector<CoeffsMatrix*> getHessianPntrs() const {return hessian_pntrs_;}
  std::vector<TargetDistribution*> getTargetDistributionPntrs() const {return targetdist_pntrs_;}
  //
  CoeffsVector* getCoeffsPntr(const unsigned int coeffs_id = 0) const {return coeffs_pntrs_[coeffs_id];}
  CoeffsVector* getTargetDistAveragesPntr(const unsigned int coeffs_id = 0) const {return targetdist_averages_pntrs_[coeffs_id];}
  CoeffsVector* getGradientPntr(const unsigned int coeffs_id = 0)const {return gradient_pntrs_[coeffs_id];}
  CoeffsMatrix* getHessianPntr(const unsigned int coeffs_id = 0) const {return hessian_pntrs_[coeffs_id];}
  //
  unsigned int getNumberOfTargetDistributionPntrs() const {return targetdist_pntrs_.size();}
  //
  size_t numberOfCoeffs(const unsigned int coeffs_id = 0) const {return coeffs_pntrs_[coeffs_id]->numberOfCoeffs();}
  size_t totalNumberOfCoeffs() const {return ncoeffs_total_;}
  unsigned int numberOfCoeffsSets() const {return ncoeffssets_;}
  double getKbT() const {return kbt_;}
  double getBeta() const;
  //
  CoeffsVector& Coeffs(const unsigned int coeffs_id = 0) const {return *coeffs_pntrs_[coeffs_id];}
  CoeffsVector& TargetDistAverages(const unsigned int coeffs_id = 0) const {return *targetdist_averages_pntrs_[coeffs_id];}
  CoeffsVector& Gradient(const unsigned int coeffs_id = 0) const {return *gradient_pntrs_[coeffs_id];}
  CoeffsMatrix& Hessian(const unsigned int coeffs_id = 0) const {return *hessian_pntrs_[coeffs_id];}
  //
  size_t getCoeffsIndex(const std::vector<unsigned int>& indices, const unsigned int coeffs_id = 0) const;
  std::vector<unsigned int> getCoeffsIndices(const size_t index, const unsigned int coeffs_id = 0) const;
  size_t getHessianIndex(const size_t index1, const size_t index2, const unsigned int coeffs_id = 0) const;
  //
  bool computeHessian() const {return compute_hessian_;}
  bool diagonalHessian() const {return diagonal_hessian_;}
  //
  bool optimizeCoeffs() const {return optimize_coeffs_;}
  Optimizer* getOptimizerPntr() const {return optimizer_pntr_;}
  bool useMultipleWalkers() const;
  //
  unsigned int getIterationCounter() const;
  //
  void updateGradientAndHessian(const bool);
  void clearGradientAndHessian() {};
  //
  virtual void updateTargetDistributions() {};
  virtual void restartTargetDistributions() {};
  //
  void linkOptimizer(Optimizer*);
  void enableHessian(const bool diagonal_hessian=true);
  void disableHessian();
  //
  void enableMultipleCoeffsSets() {use_multiple_coeffssets_=true;}
  //
  void enableDynamicTargetDistribution() {dynamic_targetdist_=true;}
  void disableDynamicTargetDistribution() {dynamic_targetdist_=false;}
  bool dynamicTargetDistribution() const {return dynamic_targetdist_;}
  //
  std::vector<unsigned int> getGridBins() const {return grid_bins_;}
  void setGridBins(const std::vector<unsigned int>&);
  void setGridBins(const unsigned int);
  std::vector<double> getGridMax() const {return grid_max_;}
  void setGridMax(const std::vector<double>&);
  std::vector<double> getGridMin() const {return grid_min_;}
  void setGridMin(const std::vector<double>&);
  //
  bool filenamesIncludeIterationNumber() const {return filenames_have_iteration_number_;}
  void enableIterationNumberInFilenames() {filenames_have_iteration_number_=true;}
  //
  std::string getIterationFilenameSuffix() const;
  std::string getCoeffsSetFilenameSuffix(const unsigned int coeffs_id) const;
  std::string getCurrentOutputFilename(const std::string&, const std::string& suffix="") const;
  std::string getBiasOutputFilename() const {return bias_filename_;}
  std::string getCurrentBiasOutputFilename(const std::string& suffix="") const;
  std::string getFesOutputFilename() const {return fes_filename_;}
  std::string getCurrentFesOutputFilename(const std::string& suffix="") const;
  std::string getTargetDistOutputFilename() const {return targetdist_filename_;}
  std::string getCurrentTargetDistOutputFilename(const std::string& suffix="") const;
  //
  void enableBiasFileOutput() {bias_fileoutput_active_=true;}
  void disableBiasFileOutput() {bias_fileoutput_active_=false;}
  bool isBiasFileOutputActive() const {return bias_fileoutput_active_;}
  std::string getBiasFileFmt() const {return bias_file_fmt_;}
  //
  void enableFesFileOutput() {fes_fileoutput_active_=true;}
  void disableFesFileOutput() {fes_fileoutput_active_=false;}
  bool isFesFileOutputActive() const {return fes_fileoutput_active_;}
  std::string getFesFileFmt() const {return fes_file_fmt_;}
  //
  void enableFesProjFileOutput() {fesproj_fileoutput_active_=true;}
  void disableFesFileProjOutput() {fesproj_fileoutput_active_=false;}
  bool isFesProjFileOutputActive() const {return fesproj_fileoutput_active_;}
  //
  void enableDynamicTargetDistFileOutput() {dynamic_targetdist_fileoutput_active_=true;}
  void disableDynamicTargetDistFileOutput() {dynamic_targetdist_fileoutput_active_=false;}
  bool isDynamicTargetDistFileOutputActive() const {return dynamic_targetdist_fileoutput_active_;}
  std::string getTargetDistFileFmt() const {return targetdist_file_fmt_;}
  std::string getTargetDistRestartFileFmt() const {return targetdist_restart_file_fmt_;}
  //
  void enableStaticTargetDistFileOutput() {static_targetdist_fileoutput_active_=true;}
  void disableStaticTargetDistFileOutput() {static_targetdist_fileoutput_active_=false;}
  bool isStaticTargetDistFileOutputActive() const {return static_targetdist_fileoutput_active_;}
  //
  std::vector< std::vector<std::string> > getProjectionArguments() const {return projection_args_;}
  std::vector<std::string> getProjectionArgument(unsigned int i) const {return projection_args_[i];}
  unsigned int getNumberOfProjectionArguments() const {return projection_args_.size();}
  //
  void setupBiasCutoff(const double, const double);
  bool biasCutoffActive() const {return bias_cutoff_active_;}
  double getBiasCutoffValue() const {return bias_cutoff_value_;}
  void setCurrentBiasMaxValue(const double max_value) {bias_current_max_value=max_value;}
  double getCurrentBiasMaxValue() const {return bias_current_max_value;}
  double getBiasCutoffSwitchingFunction(const double, double&) const;
  double getBiasCutoffSwitchingFunction(const double) const;
  void applyBiasCutoff(double&, std::vector<double>&) const;
  void applyBiasCutoff(double&, std::vector<double>&, std::vector<double>&) const;
  //
  OFile* getOFile(const std::string& filename, const bool multi_sim_single_file=false, const bool enforce_backup=true);
  //
  virtual void setupBiasFileOutput() {};
  virtual void writeBiasToFile() {};
  virtual void resetBiasFileOutput() {};
  //
  virtual void setupFesFileOutput() {};
  virtual void writeFesToFile() {};
  virtual void resetFesFileOutput() {};
  //
  virtual void setupFesProjFileOutput() {};
  virtual void writeFesProjToFile() {};
  virtual void resetFesProjFileOutput() {};
  //
  virtual void setupTargetDistFileOutput() {};
  virtual void writeTargetDistToFile() {};
  virtual void resetTargetDistFileOutput() {};
  //
  virtual void setupTargetDistProjFileOutput() {};
  virtual void writeTargetDistProjToFile() {};
  virtual void resetTargetDistProjFileOutput() {};
  //
  void updateReweightFactor();
  virtual double calculateReweightFactor() const;
  bool isReweightFactorCalculated() const {return calc_reweightfactor_;}
};


inline
size_t VesBias::getCoeffsIndex(const std::vector<unsigned int>& indices, const unsigned int coeffs_id) const {return coeffs_pntrs_[coeffs_id]->getIndex(indices);}

inline
std::vector<unsigned int> VesBias::getCoeffsIndices(const size_t index, const unsigned int coeffs_id) const {return coeffs_pntrs_[coeffs_id]->getIndices(index);}

inline
size_t VesBias::getHessianIndex(const size_t index1, const size_t index2, const unsigned int coeffs_id) const {return hessian_pntrs_[coeffs_id]->getMatrixIndex(index1,index2);}


inline
double VesBias::getBeta() const {
  plumed_massert(kbt_!=0.0,"you are requesting beta=1/(kB*T) when kB*T has not been defined. You need to give the temperature using the TEMP keyword as the MD engine does not pass it to PLUMED.");
  return 1.0/kbt_;
}


inline
std::string VesBias::getCurrentBiasOutputFilename(const std::string& suffix) const {
  return getCurrentOutputFilename(bias_filename_,suffix);
}


inline
std::string VesBias::getCurrentFesOutputFilename(const std::string& suffix) const {
  return getCurrentOutputFilename(fes_filename_,suffix);
}


inline
double VesBias::getBiasCutoffSwitchingFunction(const double bias) const {
  double dummy=0.0;
  return getBiasCutoffSwitchingFunction(bias,dummy);
}


inline
void VesBias::applyBiasCutoff(double& bias, std::vector<double>& forces) const {
  std::vector<double> dummy(0);
  applyBiasCutoff(bias,forces,dummy);
}


inline
void VesBias::applyBiasCutoff(double& bias, std::vector<double>& forces, std::vector<double>& coeffsderivs_values) const {
  double deriv_factor_sf=0.0;
  double value_sf = getBiasCutoffSwitchingFunction(bias,deriv_factor_sf);
  bias *= value_sf;
  for(unsigned int i=0; i<forces.size(); i++) {
    forces[i] *= deriv_factor_sf;
  }
  //
  for(unsigned int i=0; i<coeffsderivs_values.size(); i++) {
    coeffsderivs_values[i] *= deriv_factor_sf;
  }
}


inline
std::vector<double> VesBias::computeCovarianceFromAverages(const unsigned int c_id) const {
  size_t ncoeffs = numberOfCoeffs(c_id);
  std::vector<double> covariance(sampled_cross_averages[c_id].size(),0.0);
  // diagonal part
  for(size_t i=0; i<ncoeffs; i++) {
    size_t midx = getHessianIndex(i,i,c_id);
    covariance[midx] = sampled_cross_averages[c_id][midx] - sampled_averages[c_id][i]*sampled_averages[c_id][i];
  }
  if(!diagonal_hessian_) {
    for(size_t i=0; i<ncoeffs; i++) {
      for(size_t j=(i+1); j<ncoeffs; j++) {
        size_t midx = getHessianIndex(i,j,c_id);
        covariance[midx] = sampled_cross_averages[c_id][midx] - sampled_averages[c_id][i]*sampled_averages[c_id][j];
      }
    }
  }
  return covariance;
}


template<class T>
bool VesBias::parseMultipleValues(const std::string& keyword, std::vector<T>& values, unsigned int nvalues) {
  plumed_assert(nvalues>0);
  plumed_assert(values.size()==0);
  bool identical_values=false;
  //
  parseVector(keyword,values);
  if(values.size()==1 && nvalues>1) {
    values.resize(nvalues,values[0]);
    identical_values=true;
  }
  if(values.size()>0 && values.size()!=nvalues) {
    std::string s1; Tools::convert(nvalues,s1);
    plumed_merror("Error in " + keyword + " keyword: either give 1 common parameter value or " + s1 + " separate parameter values");
  }
  return identical_values;
}

template<class T>
bool VesBias::parseMultipleValues(const std::string& keyword, std::vector<T>& values, unsigned int nvalues, const T& default_value) {
  bool identical_values = parseMultipleValues(keyword,values,nvalues);
  if(values.size()==0) {
    values.resize(nvalues,default_value);
    identical_values=true;
  }
  return identical_values;
}


}
}

#endif
