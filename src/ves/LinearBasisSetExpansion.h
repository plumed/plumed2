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
#ifndef __PLUMED_ves_LinearBasisSetExpansion_h
#define __PLUMED_ves_LinearBasisSetExpansion_h

#include <vector>
#include <string>


namespace PLMD {

class Action;
class Keywords;
class Value;
class Communicator;
class Grid;
class OFile;


namespace ves {

class CoeffsVector;
class BasisFunctions;
class TargetDistribution;
class VesBias;


class LinearBasisSetExpansion {
  LinearBasisSetExpansion& operator=(const LinearBasisSetExpansion&) = delete;
private:
  std::string label_;
  //
  Action* action_pntr_;
  VesBias* vesbias_pntr_;
  Communicator& mycomm_;
  bool serial_;
  //
  double beta_;
  double kbt_;
  //
  std::vector<Value*> args_pntrs_;
  unsigned int nargs_;
  //
  std::vector<BasisFunctions*> basisf_pntrs_;
  std::vector<unsigned int> nbasisf_;
  //
  CoeffsVector* bias_coeffs_pntr_;
  size_t ncoeffs_;
  CoeffsVector* targetdist_averages_pntr_;
  //
  std::vector<std::string> grid_min_;
  std::vector<std::string> grid_max_;
  std::vector<unsigned int> grid_bins_;
  //
  std::string targetdist_grid_label_;
  //
  long int step_of_last_biasgrid_update;
  long int step_of_last_biaswithoutcutoffgrid_update;
  long int step_of_last_fesgrid_update;
  //
  Grid* bias_grid_pntr_;
  Grid* bias_withoutcutoff_grid_pntr_;
  Grid* fes_grid_pntr_;
  Grid* log_targetdist_grid_pntr_;
  Grid* targetdist_grid_pntr_;
  //
  TargetDistribution* targetdist_pntr_;
public:
  static void registerKeywords( Keywords& keys );
  // Constructor
  explicit LinearBasisSetExpansion(
    const std::string&,
    const double,
    Communicator&,
    std::vector<Value*>&,
    std::vector<BasisFunctions*>&,
    CoeffsVector* bias_coeffs_pntr_in=NULL);
  //
private:
  // copy constructor is disabled (private and unimplemented)
  explicit LinearBasisSetExpansion(const LinearBasisSetExpansion&);
public:
  ~LinearBasisSetExpansion();
  //
  std::vector<Value*> getPntrsToArguments() const {return args_pntrs_;}
  std::vector<BasisFunctions*> getPntrsToBasisFunctions() const {return basisf_pntrs_;}
  CoeffsVector* getPntrToBiasCoeffs() const {return bias_coeffs_pntr_;}
  Grid* getPntrToBiasGrid() const {return bias_grid_pntr_;};
  //
  unsigned int getNumberOfArguments() const {return nargs_;};
  std::vector<unsigned int> getNumberOfBasisFunctions() const {return nbasisf_;};
  size_t getNumberOfCoeffs() const {return ncoeffs_;};
  //
  CoeffsVector& BiasCoeffs() const {return *bias_coeffs_pntr_;};
  CoeffsVector& TargetDistAverages() const {return *targetdist_averages_pntr_;};
  //
  void setSerial() {serial_=true;}
  void setParallel() {serial_=false;}
  //
  void linkVesBias(VesBias*);
  void linkAction(Action*);
  // calculate bias and derivatives
  static double getBiasAndForces(const std::vector<double>&, bool&, std::vector<double>&, std::vector<double>&, std::vector<BasisFunctions*>&, CoeffsVector*, Communicator* comm_in=NULL);
  double getBiasAndForces(const std::vector<double>&, bool&, std::vector<double>&, std::vector<double>&);
  double getBiasAndForces(const std::vector<double>&, bool&, std::vector<double>&);
  double getBias(const std::vector<double>&, bool&, const bool parallel=true);
  //
  static void getBasisSetValues(const std::vector<double>&, std::vector<double>&, std::vector<BasisFunctions*>&, CoeffsVector*, Communicator* comm_in=NULL);
  void getBasisSetValues(const std::vector<double>&, std::vector<double>&, const bool parallel=true);
  //
  static double getBasisSetValue(const std::vector<double>&, const size_t, std::vector<BasisFunctions*>&, CoeffsVector*);
  double getBasisSetValue(const std::vector<double>&, const size_t);
  double getBasisSetConstant();
  // Bias grid and output stuff
  void setupBiasGrid(const bool usederiv=false);
  void updateBiasGrid();
  void resetStepOfLastBiasGridUpdate() {step_of_last_biasgrid_update = -1000;}
  void setStepOfLastBiasGridUpdate(long int step) {step_of_last_biasgrid_update = step;}
  long int getStepOfLastBiasGridUpdate() const {return step_of_last_biasgrid_update;}
  void writeBiasGridToFile(OFile&, const bool append=false) const;
  //
  void updateBiasWithoutCutoffGrid();
  void resetStepOfLastBiasWithoutCutoffGridUpdate() {step_of_last_biaswithoutcutoffgrid_update = -1000;}
  void setStepOfLastBiasWithoutCutoffGridUpdate(long int step) {step_of_last_biaswithoutcutoffgrid_update = step;}
  long int getStepOfLastBiasWithoutCutoffGridUpdate() const {return step_of_last_biaswithoutcutoffgrid_update;}
  void writeBiasWithoutCutoffGridToFile(OFile&, const bool append=false) const;
  //
  void setBiasMinimumToZero();
  void setBiasMaximumToZero();
  //
  void setupFesGrid();
  void updateFesGrid();
  void resetStepOfLastFesGridUpdate() {step_of_last_fesgrid_update = -1000;}
  void setStepOfLastFesGridUpdate(long int step) {step_of_last_fesgrid_update = step;}
  long int getStepOfLastFesGridUpdate() const {return step_of_last_fesgrid_update;}
  void writeFesGridToFile(OFile&, const bool append=false) const;
  //
  void setupFesProjGrid();
  void writeFesProjGridToFile(const std::vector<std::string>&, OFile&, const bool append=false) const;
  //
  void writeTargetDistGridToFile(OFile&, const bool append=false) const;
  void writeLogTargetDistGridToFile(OFile&, const bool append=false) const;
  void writeTargetDistProjGridToFile(const std::vector<std::string>&, OFile&, const bool append=false) const;
  void writeTargetDistributionToFile(const std::string&) const;
  //
  std::vector<unsigned int> getGridBins() const {return grid_bins_;}
  void setGridBins(const std::vector<unsigned int>&);
  void setGridBins(const unsigned int);
  //
  double getBeta() const {return beta_;}
  double getKbT() const {return kbt_;}
  double beta() const {return beta_;}
  double kBT() const {return kbt_;}
  //
  void setupUniformTargetDistribution();
  void setupTargetDistribution(TargetDistribution*);
  void updateTargetDistribution();
  //
  void readInRestartTargetDistribution(const std::string&);
  void restartTargetDistribution();
  //
  bool biasCutoffActive() const;
  //
  double calculateReweightFactor() const;
  //
private:
  //
  Grid* setupGeneralGrid(const std::string&, const bool usederiv=false);
  //
  void calculateTargetDistAveragesFromGrid(const Grid*);
  //
  bool isStaticTargetDistFileOutputActive() const;
};


inline
double LinearBasisSetExpansion::getBiasAndForces(const std::vector<double>& args_values, bool& all_inside, std::vector<double>& forces, std::vector<double>& coeffsderivs_values) {
  return getBiasAndForces(args_values,all_inside,forces,coeffsderivs_values,basisf_pntrs_, bias_coeffs_pntr_, &mycomm_);
}


inline
double LinearBasisSetExpansion::getBiasAndForces(const std::vector<double>& args_values, bool& all_inside, std::vector<double>& forces) {
  std::vector<double> coeffsderivs_values_dummy(ncoeffs_);
  return getBiasAndForces(args_values,all_inside,forces,coeffsderivs_values_dummy,basisf_pntrs_, bias_coeffs_pntr_, &mycomm_);
}


inline
double LinearBasisSetExpansion::getBias(const std::vector<double>& args_values, bool& all_inside, const bool parallel) {
  std::vector<double> forces_dummy(nargs_);
  std::vector<double> coeffsderivs_values_dummy(ncoeffs_);
  if(parallel) {
    return getBiasAndForces(args_values,all_inside,forces_dummy,coeffsderivs_values_dummy,basisf_pntrs_, bias_coeffs_pntr_, &mycomm_);
  }
  else {
    return getBiasAndForces(args_values,all_inside,forces_dummy,coeffsderivs_values_dummy,basisf_pntrs_, bias_coeffs_pntr_, NULL);
  }
}


inline
void LinearBasisSetExpansion::getBasisSetValues(const std::vector<double>& args_values, std::vector<double>& basisset_values, const bool parallel) {
  if(parallel) {
    getBasisSetValues(args_values,basisset_values,basisf_pntrs_, bias_coeffs_pntr_, &mycomm_);
  }
  else {
    getBasisSetValues(args_values,basisset_values,basisf_pntrs_, bias_coeffs_pntr_, NULL);
  }
}


inline
double LinearBasisSetExpansion::getBasisSetValue(const std::vector<double>& args_values, const size_t basisset_index) {
  return getBasisSetValue(args_values,basisset_index,basisf_pntrs_, bias_coeffs_pntr_);
}


inline
double LinearBasisSetExpansion::getBasisSetConstant() {
  std::vector<double> args_dummy(nargs_,0.0);
  return getBasisSetValue(args_dummy,0,basisf_pntrs_, bias_coeffs_pntr_);
}


}

}

#endif
