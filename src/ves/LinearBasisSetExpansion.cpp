/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2017 The VES code team
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

#include "LinearBasisSetExpansion.h"
#include "VesBias.h"
#include "CoeffsVector.h"
#include "VesTools.h"
#include "GridIntegrationWeights.h"
#include "BasisFunctions.h"
#include "TargetDistribution.h"


#include "tools/Keywords.h"
#include "tools/Grid.h"
#include "tools/Communicator.h"

#include "GridProjWeights.h"

namespace PLMD {
namespace ves {

void LinearBasisSetExpansion::registerKeywords(Keywords& keys) {
}


LinearBasisSetExpansion::LinearBasisSetExpansion(
  const std::string& label,
  const double beta_in,
  Communicator& cc,
  std::vector<Value*>& args_pntrs_in,
  std::vector<BasisFunctions*>& basisf_pntrs_in,
  CoeffsVector* bias_coeffs_pntr_in):
  label_(label),
  action_pntr_(NULL),
  vesbias_pntr_(NULL),
  mycomm_(cc),
  serial_(false),
  beta_(beta_in),
  kbt_(1.0/beta_),
  args_pntrs_(args_pntrs_in),
  nargs_(args_pntrs_.size()),
  basisf_pntrs_(basisf_pntrs_in),
  nbasisf_(basisf_pntrs_.size()),
  bias_coeffs_pntr_(bias_coeffs_pntr_in),
  ncoeffs_(0),
  targetdist_averages_pntr_(NULL),
  grid_min_(nargs_),
  grid_max_(nargs_),
  grid_bins_(nargs_,100),
  targetdist_grid_label_("targetdist"),
  step_of_last_biasgrid_update(-1000),
  step_of_last_biaswithoutcutoffgrid_update(-1000),
  step_of_last_fesgrid_update(-1000),
  bias_grid_pntr_(NULL),
  bias_withoutcutoff_grid_pntr_(NULL),
  fes_grid_pntr_(NULL),
  log_targetdist_grid_pntr_(NULL),
  targetdist_grid_pntr_(NULL),
  targetdist_pntr_(NULL)
{
  plumed_massert(args_pntrs_.size()==basisf_pntrs_.size(),"number of arguments and basis functions do not match");
  for(unsigned int k=0; k<nargs_; k++) {nbasisf_[k]=basisf_pntrs_[k]->getNumberOfBasisFunctions();}
  //
  if(bias_coeffs_pntr_==NULL) {
    bias_coeffs_pntr_ = new CoeffsVector(label_+".coeffs",args_pntrs_,basisf_pntrs_,mycomm_,true);
  }
  plumed_massert(bias_coeffs_pntr_->numberOfDimensions()==basisf_pntrs_.size(),"dimension of coeffs does not match with number of basis functions ");
  //
  ncoeffs_ = bias_coeffs_pntr_->numberOfCoeffs();
  targetdist_averages_pntr_ = new CoeffsVector(*bias_coeffs_pntr_);

  std::string targetdist_averages_label = bias_coeffs_pntr_->getLabel();
  if(targetdist_averages_label.find("coeffs")!=std::string::npos) {
    targetdist_averages_label.replace(targetdist_averages_label.find("coeffs"), std::string("coeffs").length(), "targetdist_averages");
  }
  else {
    targetdist_averages_label += "_targetdist_averages";
  }
  targetdist_averages_pntr_->setLabels(targetdist_averages_label);
  //
  for(unsigned int k=0; k<nargs_; k++) {
    grid_min_[k] = basisf_pntrs_[k]->intervalMinStr();
    grid_max_[k] = basisf_pntrs_[k]->intervalMaxStr();
  }
}

LinearBasisSetExpansion::~LinearBasisSetExpansion() {
  if(targetdist_averages_pntr_!=NULL) {
    delete targetdist_averages_pntr_;
  }
  if(bias_grid_pntr_!=NULL) {
    delete bias_grid_pntr_;
  }
  if(bias_withoutcutoff_grid_pntr_!=NULL) {
    delete bias_withoutcutoff_grid_pntr_;
  }
  if(fes_grid_pntr_!=NULL) {
    delete fes_grid_pntr_;
  }
}


bool LinearBasisSetExpansion::isStaticTargetDistFileOutputActive() const {
  bool output_static_targetdist_files=true;
  if(vesbias_pntr_!=NULL) {
    output_static_targetdist_files = vesbias_pntr_->isStaticTargetDistFileOutputActive();
  }
  return output_static_targetdist_files;
}


void LinearBasisSetExpansion::linkVesBias(VesBias* vesbias_pntr_in) {
  vesbias_pntr_ = vesbias_pntr_in;
  action_pntr_ = static_cast<Action*>(vesbias_pntr_in);

}


void LinearBasisSetExpansion::linkAction(Action* action_pntr_in) {
  action_pntr_ = action_pntr_in;
}


void LinearBasisSetExpansion::setGridBins(const std::vector<unsigned int>& grid_bins_in) {
  plumed_massert(grid_bins_in.size()==nargs_,"the number of grid bins given doesn't match the number of arguments");
  plumed_massert(bias_grid_pntr_==NULL,"setGridBins should be used before setting up the grids, otherwise it doesn't work");
  plumed_massert(fes_grid_pntr_==NULL,"setGridBins should be used before setting up the grids, otherwise it doesn't work");
  grid_bins_=grid_bins_in;
}


void LinearBasisSetExpansion::setGridBins(const unsigned int nbins) {
  std::vector<unsigned int> grid_bins_in(nargs_,nbins);
  setGridBins(grid_bins_in);
}


Grid* LinearBasisSetExpansion::setupGeneralGrid(const std::string& label_suffix, const bool usederiv) {
  bool use_spline = false;
  Grid* grid_pntr = new Grid(label_+"."+label_suffix,args_pntrs_,grid_min_,grid_max_,grid_bins_,use_spline,usederiv);
  return grid_pntr;
}


void LinearBasisSetExpansion::setupBiasGrid(const bool usederiv) {
  if(bias_grid_pntr_!=NULL) {return;}
  bias_grid_pntr_ = setupGeneralGrid("bias",usederiv);
  if(biasCutoffActive()) {
    bias_withoutcutoff_grid_pntr_ = setupGeneralGrid("bias_withoutcutoff",usederiv);
  }
}


void LinearBasisSetExpansion::setupFesGrid() {
  if(fes_grid_pntr_!=NULL) {return;}
  if(bias_grid_pntr_==NULL) {
    setupBiasGrid(true);
  }
  fes_grid_pntr_ = setupGeneralGrid("fes",false);
}


void LinearBasisSetExpansion::setupFesProjGrid() {
  if(fes_grid_pntr_==NULL) {
    setupFesGrid();
  }
}


void LinearBasisSetExpansion::updateBiasGrid() {
  plumed_massert(bias_grid_pntr_!=NULL,"the bias grid is not defined");
  if(action_pntr_!=NULL &&  getStepOfLastBiasGridUpdate()==action_pntr_->getStep()) {
    return;
  }
  for(Grid::index_t l=0; l<bias_grid_pntr_->getSize(); l++) {
    std::vector<double> forces(nargs_);
    std::vector<double> args = bias_grid_pntr_->getPoint(l);
    bool all_inside=true;
    double bias=getBiasAndForces(args,all_inside,forces);
    //
    if(biasCutoffActive()) {
      vesbias_pntr_->applyBiasCutoff(bias,forces);
    }
    //
    if(bias_grid_pntr_->hasDerivatives()) {
      bias_grid_pntr_->setValueAndDerivatives(l,bias,forces);
    }
    else {
      bias_grid_pntr_->setValue(l,bias);
    }
    //
  }
  if(vesbias_pntr_!=NULL) {
    vesbias_pntr_->setCurrentBiasMaxValue(bias_grid_pntr_->getMaxValue());
  }
  if(action_pntr_!=NULL) {
    setStepOfLastBiasGridUpdate(action_pntr_->getStep());
  }
}


void LinearBasisSetExpansion::updateBiasWithoutCutoffGrid() {
  plumed_massert(bias_withoutcutoff_grid_pntr_!=NULL,"the bias without cutoff grid is not defined");
  plumed_massert(biasCutoffActive(),"the bias cutoff has to be active");
  plumed_massert(vesbias_pntr_!=NULL,"has to be linked to a VesBias to work");
  if(action_pntr_!=NULL &&  getStepOfLastBiasWithoutCutoffGridUpdate()==action_pntr_->getStep()) {
    return;
  }
  //
  for(Grid::index_t l=0; l<bias_withoutcutoff_grid_pntr_->getSize(); l++) {
    std::vector<double> forces(nargs_);
    std::vector<double> args = bias_withoutcutoff_grid_pntr_->getPoint(l);
    bool all_inside=true;
    double bias=getBiasAndForces(args,all_inside,forces);
    if(bias_withoutcutoff_grid_pntr_->hasDerivatives()) {
      bias_withoutcutoff_grid_pntr_->setValueAndDerivatives(l,bias,forces);
    }
    else {
      bias_withoutcutoff_grid_pntr_->setValue(l,bias);
    }
  }
  //
  double bias_max = bias_withoutcutoff_grid_pntr_->getMaxValue();
  double bias_min = bias_withoutcutoff_grid_pntr_->getMinValue();
  double shift = 0.0;
  bool bias_shifted=false;
  if(bias_min < 0.0) {
    shift += -bias_min;
    bias_shifted=true;
    BiasCoeffs()[0] -= bias_min;
    bias_max -= bias_min;
  }
  if(bias_max > vesbias_pntr_->getBiasCutoffValue()) {
    shift += -(bias_max-vesbias_pntr_->getBiasCutoffValue());
    bias_shifted=true;
    BiasCoeffs()[0] -= (bias_max-vesbias_pntr_->getBiasCutoffValue());
    bias_max -= (bias_max-vesbias_pntr_->getBiasCutoffValue());
  }
  if(bias_shifted) {
    // this should be done inside a grid function really,
    // need to define my grid class for that
    for(Grid::index_t l=0; l<bias_withoutcutoff_grid_pntr_->getSize(); l++) {
      if(bias_withoutcutoff_grid_pntr_->hasDerivatives()) {
        std::vector<double> zeros(nargs_,0.0);
        bias_withoutcutoff_grid_pntr_->addValueAndDerivatives(l,shift,zeros);
      }
      else {
        bias_withoutcutoff_grid_pntr_->addValue(l,shift);
      }
    }
  }
  if(vesbias_pntr_!=NULL) {
    vesbias_pntr_->setCurrentBiasMaxValue(bias_max);
  }
  if(action_pntr_!=NULL) {
    setStepOfLastBiasWithoutCutoffGridUpdate(action_pntr_->getStep());
  }
}


void LinearBasisSetExpansion::updateFesGrid() {
  plumed_massert(fes_grid_pntr_!=NULL,"the FES grid is not defined");
  updateBiasGrid();
  if(action_pntr_!=NULL && getStepOfLastFesGridUpdate() == action_pntr_->getStep()) {
    return;
  }
  //
  double bias2fes_scalingf = -1.0;
  for(Grid::index_t l=0; l<fes_grid_pntr_->getSize(); l++) {
    double fes_value = bias2fes_scalingf*bias_grid_pntr_->getValue(l);
    if(log_targetdist_grid_pntr_!=NULL) {
      fes_value += kBT()*log_targetdist_grid_pntr_->getValue(l);
    }
    fes_grid_pntr_->setValue(l,fes_value);
  }
  fes_grid_pntr_->setMinToZero();
  if(action_pntr_!=NULL) {
    setStepOfLastFesGridUpdate(action_pntr_->getStep());
  }
}


void LinearBasisSetExpansion::writeBiasGridToFile(OFile& ofile, const bool append_file) const {
  plumed_massert(bias_grid_pntr_!=NULL,"the bias grid is not defined");
  if(append_file) {ofile.enforceRestart();}
  bias_grid_pntr_->writeToFile(ofile);
}


void LinearBasisSetExpansion::writeBiasWithoutCutoffGridToFile(OFile& ofile, const bool append_file) const {
  plumed_massert(bias_withoutcutoff_grid_pntr_!=NULL,"the bias without cutoff grid is not defined");
  if(append_file) {ofile.enforceRestart();}
  bias_withoutcutoff_grid_pntr_->writeToFile(ofile);
}


void LinearBasisSetExpansion::writeFesGridToFile(OFile& ofile, const bool append_file) const {
  plumed_massert(fes_grid_pntr_!=NULL,"the FES grid is not defined");
  if(append_file) {ofile.enforceRestart();}
  fes_grid_pntr_->writeToFile(ofile);
}


void LinearBasisSetExpansion::writeFesProjGridToFile(const std::vector<std::string>& proj_arg, OFile& ofile, const bool append_file) const {
  plumed_massert(fes_grid_pntr_!=NULL,"the FES grid is not defined");
  FesWeight* Fw = new FesWeight(beta_);
  Grid proj_grid = fes_grid_pntr_->project(proj_arg,Fw);
  proj_grid.setMinToZero();
  if(append_file) {ofile.enforceRestart();}
  proj_grid.writeToFile(ofile);
  delete Fw;
}


void LinearBasisSetExpansion::writeTargetDistGridToFile(OFile& ofile, const bool append_file) const {
  if(targetdist_grid_pntr_==NULL) {return;}
  if(append_file) {ofile.enforceRestart();}
  targetdist_grid_pntr_->writeToFile(ofile);
}


void LinearBasisSetExpansion::writeLogTargetDistGridToFile(OFile& ofile, const bool append_file) const {
  if(log_targetdist_grid_pntr_==NULL) {return;}
  if(append_file) {ofile.enforceRestart();}
  log_targetdist_grid_pntr_->writeToFile(ofile);
}


void LinearBasisSetExpansion::writeTargetDistProjGridToFile(const std::vector<std::string>& proj_arg, OFile& ofile, const bool append_file) const {
  if(targetdist_grid_pntr_==NULL) {return;}
  if(append_file) {ofile.enforceRestart();}
  Grid proj_grid = TargetDistribution::getMarginalDistributionGrid(targetdist_grid_pntr_,proj_arg);
  proj_grid.writeToFile(ofile);
}


void LinearBasisSetExpansion::writeTargetDistributionToFile(const std::string& filename) const {
  OFile of1; OFile of2;
  if(action_pntr_!=NULL) {
    of1.link(*action_pntr_); of2.link(*action_pntr_);
  }
  of1.enforceBackup(); of2.enforceBackup();
  of1.open(filename);
  of2.open(FileBase::appendSuffix(filename,".log"));
  writeTargetDistGridToFile(of1);
  writeLogTargetDistGridToFile(of2);
  of1.close(); of2.close();
}


double LinearBasisSetExpansion::getBiasAndForces(const std::vector<double>& args_values, bool& all_inside, std::vector<double>& forces, std::vector<double>& coeffsderivs_values, std::vector<BasisFunctions*>& basisf_pntrs_in, CoeffsVector* coeffs_pntr_in, Communicator* comm_in) {
  unsigned int nargs = args_values.size();
  plumed_assert(coeffs_pntr_in->numberOfDimensions()==nargs);
  plumed_assert(basisf_pntrs_in.size()==nargs);
  plumed_assert(forces.size()==nargs);
  plumed_assert(coeffsderivs_values.size()==coeffs_pntr_in->numberOfCoeffs());

  std::vector<double> args_values_trsfrm(nargs);
  // std::vector<bool>   inside_interval(nargs,true);
  all_inside = true;
  //
  std::vector< std::vector <double> > bf_values(nargs);
  std::vector< std::vector <double> > bf_derivs(nargs);
  //
  for(unsigned int k=0; k<nargs; k++) {
    bf_values[k].assign(basisf_pntrs_in[k]->getNumberOfBasisFunctions(),0.0);
    bf_derivs[k].assign(basisf_pntrs_in[k]->getNumberOfBasisFunctions(),0.0);
    bool curr_inside=true;
    basisf_pntrs_in[k]->getAllValues(args_values[k],args_values_trsfrm[k],curr_inside,bf_values[k],bf_derivs[k]);
    // inside_interval[k]=curr_inside;
    if(!curr_inside) {all_inside=false;}
    forces[k]=0.0;
  }
  //
  size_t stride=1;
  size_t rank=0;
  if(comm_in!=NULL)
  {
    stride=comm_in->Get_size();
    rank=comm_in->Get_rank();
  }
  // loop over coeffs
  double bias=0.0;
  for(size_t i=rank; i<coeffs_pntr_in->numberOfCoeffs(); i+=stride) {
    std::vector<unsigned int> indices=coeffs_pntr_in->getIndices(i);
    double coeff = coeffs_pntr_in->getValue(i);
    double bf_curr=1.0;
    for(unsigned int k=0; k<nargs; k++) {
      bf_curr*=bf_values[k][indices[k]];
    }
    bias+=coeff*bf_curr;
    coeffsderivs_values[i] = bf_curr;
    for(unsigned int k=0; k<nargs; k++) {
      double der = 1.0;
      for(unsigned int l=0; l<nargs; l++) {
        if(l!=k) {der*=bf_values[l][indices[l]];}
        else {der*=bf_derivs[l][indices[l]];}
      }
      forces[k]-=coeff*der;
      // maybe faster but dangerous
      // forces[k]-=coeff*bf_curr*(bf_derivs[k][indices[k]]/bf_values[k][indices[k]]);
    }
  }
  //
  if(comm_in!=NULL) {
    comm_in->Sum(bias);
    comm_in->Sum(forces);
  }
  return bias;
}


void LinearBasisSetExpansion::getBasisSetValues(const std::vector<double>& args_values, std::vector<double>& basisset_values, std::vector<BasisFunctions*>& basisf_pntrs_in, CoeffsVector* coeffs_pntr_in, Communicator* comm_in) {
  unsigned int nargs = args_values.size();
  plumed_assert(coeffs_pntr_in->numberOfDimensions()==nargs);
  plumed_assert(basisf_pntrs_in.size()==nargs);

  std::vector<double> args_values_trsfrm(nargs);
  std::vector< std::vector <double> > bf_values;
  //
  for(unsigned int k=0; k<nargs; k++) {
    std::vector<double> tmp_val(basisf_pntrs_in[k]->getNumberOfBasisFunctions());
    std::vector<double> tmp_der(tmp_val.size());
    bool inside=true;
    basisf_pntrs_in[k]->getAllValues(args_values[k],args_values_trsfrm[k],inside,tmp_val,tmp_der);
    bf_values.push_back(tmp_val);
  }
  //
  size_t stride=1;
  size_t rank=0;
  if(comm_in!=NULL)
  {
    stride=comm_in->Get_size();
    rank=comm_in->Get_rank();
  }
  // loop over basis set
  for(size_t i=rank; i<coeffs_pntr_in->numberOfCoeffs(); i+=stride) {
    std::vector<unsigned int> indices=coeffs_pntr_in->getIndices(i);
    double bf_curr=1.0;
    for(unsigned int k=0; k<nargs; k++) {
      bf_curr*=bf_values[k][indices[k]];
    }
    basisset_values[i] = bf_curr;
  }
  //
  if(comm_in!=NULL) {
    comm_in->Sum(basisset_values);
  }
}


void LinearBasisSetExpansion::setupUniformTargetDistribution() {
  std::vector< std::vector <double> > bf_integrals(0);
  std::vector<double> targetdist_averages(ncoeffs_,0.0);
  //
  for(unsigned int k=0; k<nargs_; k++) {
    bf_integrals.push_back(basisf_pntrs_[k]->getUniformIntegrals());
  }
  //
  for(size_t i=0; i<ncoeffs_; i++) {
    std::vector<unsigned int> indices=bias_coeffs_pntr_->getIndices(i);
    double value = 1.0;
    for(unsigned int k=0; k<nargs_; k++) {
      value*=bf_integrals[k][indices[k]];
    }
    targetdist_averages[i]=value;
  }
  TargetDistAverages() = targetdist_averages;
}


void LinearBasisSetExpansion::setupTargetDistribution(TargetDistribution* targetdist_pntr_in) {
  targetdist_pntr_ = targetdist_pntr_in;
  //
  targetdist_pntr_->setupGrids(args_pntrs_,grid_min_,grid_max_,grid_bins_);
  targetdist_grid_pntr_      = targetdist_pntr_->getTargetDistGridPntr();
  log_targetdist_grid_pntr_  = targetdist_pntr_->getLogTargetDistGridPntr();
  //
  if(targetdist_pntr_->isDynamic()) {
    vesbias_pntr_->enableDynamicTargetDistribution();
  }
  //
  if(targetdist_pntr_->biasGridNeeded()) {
    setupBiasGrid(true);
    targetdist_pntr_->linkBiasGrid(bias_grid_pntr_);
  }
  if(targetdist_pntr_->biasWithoutCutoffGridNeeded()) {
    setupBiasGrid(true);
    targetdist_pntr_->linkBiasWithoutCutoffGrid(bias_withoutcutoff_grid_pntr_);
  }
  if(targetdist_pntr_->fesGridNeeded()) {
    setupFesGrid();
    targetdist_pntr_->linkFesGrid(fes_grid_pntr_);
  }
  //
  targetdist_pntr_->updateTargetDist();
  calculateTargetDistAveragesFromGrid(targetdist_grid_pntr_);
}


void LinearBasisSetExpansion::updateTargetDistribution() {
  plumed_massert(targetdist_pntr_!=NULL,"the target distribution hasn't been setup!");
  plumed_massert(targetdist_pntr_->isDynamic(),"this should only be used for dynamically updated target distributions!");
  if(targetdist_pntr_->biasGridNeeded()) {updateBiasGrid();}
  if(biasCutoffActive()) {updateBiasWithoutCutoffGrid();}
  if(targetdist_pntr_->fesGridNeeded()) {updateFesGrid();}
  targetdist_pntr_->updateTargetDist();
  calculateTargetDistAveragesFromGrid(targetdist_grid_pntr_);
}


void LinearBasisSetExpansion::readInRestartTargetDistribution(const std::string& grid_fname) {
  targetdist_pntr_->readInRestartTargetDistGrid(grid_fname);
  if(biasCutoffActive()) {targetdist_pntr_->clearLogTargetDistGrid();}
}


void LinearBasisSetExpansion::restartTargetDistribution() {
  plumed_massert(targetdist_pntr_!=NULL,"the target distribution hasn't been setup!");
  plumed_massert(targetdist_pntr_->isDynamic(),"this should only be used for dynamically updated target distributions!");
  if(biasCutoffActive()) {updateBiasWithoutCutoffGrid();}
  calculateTargetDistAveragesFromGrid(targetdist_grid_pntr_);
}


void LinearBasisSetExpansion::calculateTargetDistAveragesFromGrid(const Grid* targetdist_grid_pntr) {
  plumed_assert(targetdist_grid_pntr!=NULL);
  std::vector<double> targetdist_averages(ncoeffs_,0.0);
  std::vector<double> integration_weights = GridIntegrationWeights::getIntegrationWeights(targetdist_grid_pntr);
  Grid::index_t stride=mycomm_.Get_size();
  Grid::index_t rank=mycomm_.Get_rank();
  for(Grid::index_t l=rank; l<targetdist_grid_pntr->getSize(); l+=stride) {
    std::vector<double> args_values = targetdist_grid_pntr->getPoint(l);
    std::vector<double> basisset_values(ncoeffs_);
    getBasisSetValues(args_values,basisset_values,false);
    double weight = integration_weights[l]*targetdist_grid_pntr->getValue(l);
    for(unsigned int i=0; i<ncoeffs_; i++) {
      targetdist_averages[i] += weight*basisset_values[i];
    }
  }
  mycomm_.Sum(targetdist_averages);
  // the overall constant;
  targetdist_averages[0] = 1.0;
  TargetDistAverages() = targetdist_averages;
}


void LinearBasisSetExpansion::setBiasMinimumToZero() {
  plumed_massert(bias_grid_pntr_!=NULL,"setBiasMinimumToZero can only be used if the bias grid is defined");
  updateBiasGrid();
  BiasCoeffs()[0]-=bias_grid_pntr_->getMinValue();
}


void LinearBasisSetExpansion::setBiasMaximumToZero() {
  plumed_massert(bias_grid_pntr_!=NULL,"setBiasMaximumToZero can only be used if the bias grid is defined");
  updateBiasGrid();
  BiasCoeffs()[0]-=bias_grid_pntr_->getMaxValue();
}


bool LinearBasisSetExpansion::biasCutoffActive() const {
  if(vesbias_pntr_!=NULL) {return vesbias_pntr_->biasCutoffActive();}
  else {return false;}
}


double LinearBasisSetExpansion::calculateReweightFactor() const {
  plumed_massert(targetdist_grid_pntr_!=NULL,"calculateReweightFactor only be used if the target distribution grid is defined");
  plumed_massert(bias_grid_pntr_!=NULL,"calculateReweightFactor only be used if the bias grid is defined");
  double sum = 0.0;
  std::vector<double> integration_weights = GridIntegrationWeights::getIntegrationWeights(targetdist_grid_pntr_);
  //
  for(Grid::index_t l=0; l<targetdist_grid_pntr_->getSize(); l++) {
    sum += integration_weights[l] * targetdist_grid_pntr_->getValue(l) * exp(+beta_*bias_grid_pntr_->getValue(l));
  }
  return (1.0/beta_)*std::log(sum);
}




}

}
