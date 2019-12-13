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

#include "VesBias.h"
#include "BasisFunctions.h"
#include "CoeffsVector.h"
#include "CoeffsMatrix.h"
#include "Optimizer.h"
#include "FermiSwitchingFunction.h"
#include "VesTools.h"
#include "TargetDistribution.h"

#include "tools/Communicator.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/File.h"


namespace PLMD {
namespace ves {

VesBias::VesBias(const ActionOptions&ao):
  Action(ao),
  Bias(ao),
  ncoeffssets_(0),
  coeffs_pntrs_(0),
  targetdist_averages_pntrs_(0),
  gradient_pntrs_(0),
  hessian_pntrs_(0),
  sampled_averages(0),
  sampled_cross_averages(0),
  use_multiple_coeffssets_(false),
  coeffs_fnames(0),
  ncoeffs_total_(0),
  optimizer_pntr_(NULL),
  optimize_coeffs_(false),
  compute_hessian_(false),
  diagonal_hessian_(true),
  aver_counters(0),
  kbt_(0.0),
  targetdist_pntrs_(0),
  dynamic_targetdist_(false),
  grid_bins_(0),
  grid_min_(0),
  grid_max_(0),
  bias_filename_(""),
  fes_filename_(""),
  targetdist_filename_(""),
  coeffs_id_prefix_("c-"),
  bias_file_fmt_("14.9f"),
  fes_file_fmt_("14.9f"),
  targetdist_file_fmt_("14.9f"),
  targetdist_restart_file_fmt_("30.16e"),
  filenames_have_iteration_number_(false),
  bias_fileoutput_active_(false),
  fes_fileoutput_active_(false),
  fesproj_fileoutput_active_(false),
  dynamic_targetdist_fileoutput_active_(false),
  static_targetdist_fileoutput_active_(true),
  bias_cutoff_active_(false),
  bias_cutoff_value_(0.0),
  bias_current_max_value(0.0),
  bias_cutoff_swfunc_pntr_(NULL),
  calc_reweightfactor_(false)
{
  log.printf("  VES bias, please read and cite ");
  log << plumed.cite("Valsson and Parrinello, Phys. Rev. Lett. 113, 090601 (2014)");
  log.printf("\n");

  double temp=0.0;
  parse("TEMP",temp);
  if(temp>0.0) {
    kbt_=plumed.getAtoms().getKBoltzmann()*temp;
  }
  else {
    kbt_=plumed.getAtoms().getKbT();
  }
  if(kbt_>0.0) {
    log.printf("  KbT: %f\n",kbt_);
  }
  // NOTE: the check for that the temperature is given is done when linking the optimizer later on.

  if(keywords.exists("COEFFS")) {
    parseVector("COEFFS",coeffs_fnames);
  }

  if(keywords.exists("GRID_BINS")) {
    parseMultipleValues<unsigned int>("GRID_BINS",grid_bins_,getNumberOfArguments(),100);
  }

  if(keywords.exists("GRID_MIN") && keywords.exists("GRID_MAX")) {
    parseMultipleValues("GRID_MIN",grid_min_,getNumberOfArguments());
    parseMultipleValues("GRID_MAX",grid_max_,getNumberOfArguments());
  }

  std::vector<std::string> targetdist_labels;
  if(keywords.exists("TARGET_DISTRIBUTION")) {
    parseVector("TARGET_DISTRIBUTION",targetdist_labels);
    if(targetdist_labels.size()>1) {
      plumed_merror(getName()+" with label "+getLabel()+": multiple target distribution labels not allowed");
    }
  }
  else if(keywords.exists("TARGET_DISTRIBUTIONS")) {
    parseVector("TARGET_DISTRIBUTIONS",targetdist_labels);
  }

  std::string error_msg = "";
  targetdist_pntrs_ = VesTools::getPointersFromLabels<TargetDistribution*>(targetdist_labels,plumed.getActionSet(),error_msg);
  if(error_msg.size()>0) {plumed_merror("Problem with target distribution in "+getName()+": "+error_msg);}

  for(unsigned int i=0; i<targetdist_pntrs_.size(); i++) {
    targetdist_pntrs_[i]->linkVesBias(this);
  }


  if(getNumberOfArguments()>2) {
    disableStaticTargetDistFileOutput();
  }


  if(keywords.exists("BIAS_FILE")) {
    parse("BIAS_FILE",bias_filename_);
    if(bias_filename_.size()==0) {
      bias_filename_ = "bias." + getLabel() + ".data";
    }
  }
  if(keywords.exists("FES_FILE")) {
    parse("FES_FILE",fes_filename_);
    if(fes_filename_.size()==0) {
      fes_filename_ = "fes." + getLabel() + ".data";
    }
  }
  if(keywords.exists("TARGETDIST_FILE")) {
    parse("TARGETDIST_FILE",targetdist_filename_);
    if(targetdist_filename_.size()==0) {
      targetdist_filename_ = "targetdist." + getLabel() + ".data";
    }
  }
  //
  if(keywords.exists("BIAS_FILE_FMT")) {
    parse("BIAS_FILE_FMT",bias_file_fmt_);
  }
  if(keywords.exists("FES_FILE_FMT")) {
    parse("FES_FILE_FMT",fes_file_fmt_);
  }
  if(keywords.exists("TARGETDIST_FILE_FMT")) {
    parse("TARGETDIST_FILE_FMT",targetdist_file_fmt_);
  }
  if(keywords.exists("TARGETDIST_RESTART_FILE_FMT")) {
    parse("TARGETDIST_RESTART_FILE_FMT",targetdist_restart_file_fmt_);
  }

  //
  if(keywords.exists("BIAS_CUTOFF")) {
    double cutoff_value=0.0;
    parse("BIAS_CUTOFF",cutoff_value);
    if(cutoff_value<0.0) {
      plumed_merror("the value given in BIAS_CUTOFF doesn't make sense, it should be larger than 0.0");
    }
    //
    if(cutoff_value>0.0) {
      double fermi_lambda=10.0;
      parse("BIAS_CUTOFF_FERMI_LAMBDA",fermi_lambda);
      setupBiasCutoff(cutoff_value,fermi_lambda);
      log.printf("  Employing a bias cutoff of %f (the lambda value for the Fermi switching function is %f), see and cite ",cutoff_value,fermi_lambda);
      log << plumed.cite("McCarty, Valsson, Tiwary, and Parrinello, Phys. Rev. Lett. 115, 070601 (2015)");
      log.printf("\n");
    }
  }


  if(keywords.exists("PROJ_ARG")) {
    std::vector<std::string> proj_arg;
    for(int i=1;; i++) {
      if(!parseNumberedVector("PROJ_ARG",i,proj_arg)) {break;}
      // checks
      if(proj_arg.size() > (getNumberOfArguments()-1) ) {
        plumed_merror("PROJ_ARG must be a subset of ARG");
      }
      //
      for(unsigned int k=0; k<proj_arg.size(); k++) {
        bool found = false;
        for(unsigned int l=0; l<getNumberOfArguments(); l++) {
          if(proj_arg[k]==getPntrToArgument(l)->getName()) {
            found = true;
            break;
          }
        }
        if(!found) {
          std::string s1; Tools::convert(i,s1);
          std::string error = "PROJ_ARG" + s1 + ": label " + proj_arg[k] + " is not among the arguments given in ARG";
          plumed_merror(error);
        }
      }
      //
      projection_args_.push_back(proj_arg);
    }
  }

  if(keywords.exists("CALC_REWEIGHT_FACTOR")) {
    parseFlag("CALC_REWEIGHT_FACTOR",calc_reweightfactor_);
    if(calc_reweightfactor_) {
      addComponent("rct"); componentIsNotPeriodic("rct");
      updateReweightFactor();
    }
  }


}


VesBias::~VesBias() {
  for(unsigned int i=0; i<coeffs_pntrs_.size(); i++) {
    delete coeffs_pntrs_[i];
  }
  for(unsigned int i=0; i<targetdist_averages_pntrs_.size(); i++) {
    delete targetdist_averages_pntrs_[i];
  }
  for(unsigned int i=0; i<gradient_pntrs_.size(); i++) {
    delete gradient_pntrs_[i];
  }
  for(unsigned int i=0; i<hessian_pntrs_.size(); i++) {
    delete hessian_pntrs_[i];
  }
  if(bias_cutoff_swfunc_pntr_!=NULL) {
    delete bias_cutoff_swfunc_pntr_;
  }
}


void VesBias::registerKeywords( Keywords& keys ) {
  Bias::registerKeywords(keys);
  keys.add("optional","TEMP","the system temperature - this is needed if the MD code does not pass the temperature to PLUMED.");
  //
  keys.reserve("optional","COEFFS","read in the coefficients from files.");
  //
  keys.reserve("optional","TARGET_DISTRIBUTION","the label of the target distribution to be used.");
  keys.reserve("optional","TARGET_DISTRIBUTIONS","the label of the target distribution to be used. Here you are allows to use multiple labels.");
  //
  keys.reserve("optional","GRID_BINS","the number of bins used for the grid. The default value is 100 bins per dimension.");
  keys.reserve("optional","GRID_MIN","the lower bounds used for the grid.");
  keys.reserve("optional","GRID_MAX","the upper bounds used for the grid.");
  //
  keys.add("optional","BIAS_FILE","filename of the file on which the bias should be written out. By default it is bias.LABEL.data. Note that suffixes indicating the iteration number (iter-#) are added to the filename when optimizing coefficients.");
  keys.add("optional","FES_FILE","filename of the file on which the FES should be written out. By default it is fes.LABEL.data. Note that suffixes indicating the iteration number (iter-#) are added to the filename when optimizing coefficients.");
  keys.add("optional","TARGETDIST_FILE","filename of the file on which the target distribution should be written out. By default it is targetdist.LABEL.data. Note that suffixes indicating the iteration number (iter-#) are added to the filename when optimizing coefficients and the target distribution is dynamic.");
  //
  // keys.add("optional","BIAS_FILE_FMT","the format of the bias files, by default it is %14.9f.");
  // keys.add("optional","FES_FILE_FMT","the format of the FES files, by default it is %14.9f.");
  // keys.add("optional","TARGETDIST_FILE_FMT","the format of the target distribution files, by default it is %14.9f.");
  // keys.add("hidden","TARGETDIST_RESTART_FILE_FMT","the format of the target distribution files that are used for restarting, by default it is %30.16e.");
  //
  keys.reserve("optional","BIAS_CUTOFF","cutoff the bias such that it only fills the free energy surface up to certain level F_cutoff, here you should give the value of the F_cutoff.");
  keys.reserve("optional","BIAS_CUTOFF_FERMI_LAMBDA","the lambda value used in the Fermi switching function for the bias cutoff (BIAS_CUTOFF), the default value is 10.0.");
  //
  keys.reserve("numbered","PROJ_ARG","arguments for doing projections of the FES or the target distribution.");
  //
  keys.reserveFlag("CALC_REWEIGHT_FACTOR",false,"enable the calculation of the reweight factor c(t). You should also give a stride for updating the reweight factor in the optimizer by using the REWEIGHT_FACTOR_STRIDE keyword if the coefficients are updated.");

}


void VesBias::useInitialCoeffsKeywords(Keywords& keys) {
  keys.use("COEFFS");
}


void VesBias::useTargetDistributionKeywords(Keywords& keys) {
  plumed_massert(!keys.exists("TARGET_DISTRIBUTIONS"),"you cannot use both useTargetDistributionKeywords and useMultipleTargetDistributionKeywords");
  keys.use("TARGET_DISTRIBUTION");
}


void VesBias::useMultipleTargetDistributionKeywords(Keywords& keys) {
  plumed_massert(!keys.exists("TARGET_DISTRIBUTION"),"you cannot use both useTargetDistributionKeywords and useMultipleTargetDistributionKeywords");
  keys.use("TARGET_DISTRIBUTIONS");
}


void VesBias::useGridBinKeywords(Keywords& keys) {
  keys.use("GRID_BINS");
}


void VesBias::useGridLimitsKeywords(Keywords& keys) {
  keys.use("GRID_MIN");
  keys.use("GRID_MAX");
}


void VesBias::useBiasCutoffKeywords(Keywords& keys) {
  keys.use("BIAS_CUTOFF");
  keys.use("BIAS_CUTOFF_FERMI_LAMBDA");
}


void VesBias::useProjectionArgKeywords(Keywords& keys) {
  keys.use("PROJ_ARG");
}


void VesBias::useReweightFactorKeywords(Keywords& keys) {
  keys.use("CALC_REWEIGHT_FACTOR");
  keys.addOutputComponent("rct","CALC_REWEIGHT_FACTOR","the reweight factor c(t).");
}


void VesBias::addCoeffsSet(const std::vector<std::string>& dimension_labels,const std::vector<unsigned int>& indices_shape) {
  CoeffsVector* coeffs_pntr_tmp = new CoeffsVector("coeffs",dimension_labels,indices_shape,comm,true);
  initializeCoeffs(coeffs_pntr_tmp);
}


void VesBias::addCoeffsSet(std::vector<Value*>& args,std::vector<BasisFunctions*>& basisf) {
  CoeffsVector* coeffs_pntr_tmp = new CoeffsVector("coeffs",args,basisf,comm,true);
  initializeCoeffs(coeffs_pntr_tmp);
}


void VesBias::addCoeffsSet(CoeffsVector* coeffs_pntr_in) {
  initializeCoeffs(coeffs_pntr_in);
}


void VesBias::initializeCoeffs(CoeffsVector* coeffs_pntr_in) {
  //
  coeffs_pntr_in->linkVesBias(this);
  //
  std::string label;
  if(!use_multiple_coeffssets_ && ncoeffssets_==1) {
    plumed_merror("you are not allowed to use multiple coefficient sets");
  }
  //
  label = getCoeffsSetLabelString("coeffs",ncoeffssets_);
  coeffs_pntr_in->setLabels(label);

  coeffs_pntrs_.push_back(coeffs_pntr_in);
  CoeffsVector* aver_ps_tmp = new CoeffsVector(*coeffs_pntr_in);
  label = getCoeffsSetLabelString("targetdist_averages",ncoeffssets_);
  aver_ps_tmp->setLabels(label);
  aver_ps_tmp->setValues(0.0);
  targetdist_averages_pntrs_.push_back(aver_ps_tmp);
  //
  CoeffsVector* gradient_tmp = new CoeffsVector(*coeffs_pntr_in);
  label = getCoeffsSetLabelString("gradient",ncoeffssets_);
  gradient_tmp->setLabels(label);
  gradient_pntrs_.push_back(gradient_tmp);
  //
  label = getCoeffsSetLabelString("hessian",ncoeffssets_);
  CoeffsMatrix* hessian_tmp = new CoeffsMatrix(label,coeffs_pntr_in,comm,diagonal_hessian_);
  hessian_pntrs_.push_back(hessian_tmp);
  //
  std::vector<double> aver_sampled_tmp;
  aver_sampled_tmp.assign(coeffs_pntr_in->numberOfCoeffs(),0.0);
  sampled_averages.push_back(aver_sampled_tmp);
  //
  std::vector<double> cross_aver_sampled_tmp;
  cross_aver_sampled_tmp.assign(hessian_tmp->getSize(),0.0);
  sampled_cross_averages.push_back(cross_aver_sampled_tmp);
  //
  aver_counters.push_back(0);
  //
  ncoeffssets_++;
}


bool VesBias::readCoeffsFromFiles() {
  plumed_assert(ncoeffssets_>0);
  plumed_massert(keywords.exists("COEFFS"),"you are not allowed to use this function as the COEFFS keyword is not enabled");
  bool read_coeffs = false;
  if(coeffs_fnames.size()>0) {
    plumed_massert(coeffs_fnames.size()==ncoeffssets_,"COEFFS keyword is of the wrong size");
    if(ncoeffssets_==1) {
      log.printf("  Read in coefficients from file ");
    }
    else {
      log.printf("  Read in coefficients from files:\n");
    }
    for(unsigned int i=0; i<ncoeffssets_; i++) {
      IFile ifile;
      ifile.link(*this);
      ifile.open(coeffs_fnames[i]);
      if(!ifile.FieldExist(coeffs_pntrs_[i]->getDataLabel())) {
        std::string error_msg = "Problem with reading coefficients from file " + ifile.getPath() + ": no field with name " + coeffs_pntrs_[i]->getDataLabel() + "\n";
        plumed_merror(error_msg);
      }
      size_t ncoeffs_read = coeffs_pntrs_[i]->readFromFile(ifile,false,false);
      coeffs_pntrs_[i]->setIterationCounterAndTime(0,getTime());
      if(ncoeffssets_==1) {
        log.printf("%s (read %zu of %zu values)\n", ifile.getPath().c_str(),ncoeffs_read,coeffs_pntrs_[i]->numberOfCoeffs());
      }
      else {
        log.printf("   coefficient %u: %s (read %zu of %zu values)\n",i,ifile.getPath().c_str(),ncoeffs_read,coeffs_pntrs_[i]->numberOfCoeffs());
      }
      ifile.close();
    }
    read_coeffs = true;
  }
  return read_coeffs;
}


void VesBias::updateGradientAndHessian(const bool use_mwalkers_mpi) {
  for(unsigned int k=0; k<ncoeffssets_; k++) {
    //
    comm.Sum(sampled_averages[k]);
    comm.Sum(sampled_cross_averages[k]);
    if(use_mwalkers_mpi) {
      double walker_weight=1.0;
      if(aver_counters[k]==0) {walker_weight=0.0;}
      multiSimSumAverages(k,walker_weight);
    }
    // NOTE: this assumes that all walkers have the same TargetDist, might change later on!!
    Gradient(k).setValues( TargetDistAverages(k) - sampled_averages[k] );
    Hessian(k) = computeCovarianceFromAverages(k);
    Hessian(k) *= getBeta();
    //
    Gradient(k).activate();
    Hessian(k).activate();
    //
    // Check the total number of samples (from all walkers) and deactivate the Gradient and Hessian if it
    // is zero
    unsigned int total_samples = aver_counters[k];
    if(use_mwalkers_mpi) {
      if(comm.Get_rank()==0) {multi_sim_comm.Sum(total_samples);}
      comm.Bcast(total_samples,0);
    }
    if(total_samples==0) {
      Gradient(k).deactivate();
      Gradient(k).clear();
      Hessian(k).deactivate();
      Hessian(k).clear();
    }
    //
    std::fill(sampled_averages[k].begin(), sampled_averages[k].end(), 0.0);
    std::fill(sampled_cross_averages[k].begin(), sampled_cross_averages[k].end(), 0.0);
    aver_counters[k]=0;
  }
}


void VesBias::multiSimSumAverages(const unsigned int c_id, const double walker_weight) {
  plumed_massert(walker_weight>=0.0,"the weight of the walker cannot be negative!");
  if(walker_weight!=1.0) {
    for(size_t i=0; i<sampled_averages[c_id].size(); i++) {
      sampled_averages[c_id][i] *= walker_weight;
    }
    for(size_t i=0; i<sampled_cross_averages[c_id].size(); i++) {
      sampled_cross_averages[c_id][i] *= walker_weight;
    }
  }
  //
  if(comm.Get_rank()==0) {
    multi_sim_comm.Sum(sampled_averages[c_id]);
    multi_sim_comm.Sum(sampled_cross_averages[c_id]);
    double norm_weights = walker_weight;
    multi_sim_comm.Sum(norm_weights);
    if(norm_weights>0.0) {norm_weights=1.0/norm_weights;}
    for(size_t i=0; i<sampled_averages[c_id].size(); i++) {
      sampled_averages[c_id][i] *= norm_weights;
    }
    for(size_t i=0; i<sampled_cross_averages[c_id].size(); i++) {
      sampled_cross_averages[c_id][i] *= norm_weights;
    }
  }
  comm.Bcast(sampled_averages[c_id],0);
  comm.Bcast(sampled_cross_averages[c_id],0);
}


void VesBias::addToSampledAverages(const std::vector<double>& values, const unsigned int c_id) {
  /*
  use the following online equation to calculate the average and covariance
  (see https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Covariance)
      xm[n+1] = xm[n] + (x[n+1]-xm[n])/(n+1)
  */
  double counter_dbl = static_cast<double>(aver_counters[c_id]);
  size_t ncoeffs = numberOfCoeffs(c_id);
  std::vector<double> deltas(ncoeffs,0.0);
  size_t stride = comm.Get_size();
  size_t rank = comm.Get_rank();
  // update average and diagonal part of Hessian
  for(size_t i=rank; i<ncoeffs; i+=stride) {
    size_t midx = getHessianIndex(i,i,c_id);
    deltas[i] = (values[i]-sampled_averages[c_id][i])/(counter_dbl+1); // (x[n+1]-xm[n])/(n+1)
    sampled_averages[c_id][i] += deltas[i];
    sampled_cross_averages[c_id][midx] += (values[i]*values[i]-sampled_cross_averages[c_id][midx])/(counter_dbl+1);
  }
  comm.Sum(deltas);
  // update off-diagonal part of the Hessian
  if(!diagonal_hessian_) {
    for(size_t i=rank; i<ncoeffs; i+=stride) {
      for(size_t j=(i+1); j<ncoeffs; j++) {
        size_t midx = getHessianIndex(i,j,c_id);
        sampled_cross_averages[c_id][midx] += (values[i]*values[j]-sampled_cross_averages[c_id][midx])/(counter_dbl+1);
      }
    }
  }
  // NOTE: the MPI sum for sampled_averages and sampled_cross_averages is done later
  aver_counters[c_id] += 1;
}


void VesBias::setTargetDistAverages(const std::vector<double>& coeffderivs_aver_ps, const unsigned int coeffs_id) {
  TargetDistAverages(coeffs_id) = coeffderivs_aver_ps;
  TargetDistAverages(coeffs_id).setIterationCounterAndTime(this->getIterationCounter(),this->getTime());
}


void VesBias::setTargetDistAverages(const CoeffsVector& coeffderivs_aver_ps, const unsigned int coeffs_id) {
  TargetDistAverages(coeffs_id).setValues( coeffderivs_aver_ps );
  TargetDistAverages(coeffs_id).setIterationCounterAndTime(this->getIterationCounter(),this->getTime());
}


void VesBias::setTargetDistAveragesToZero(const unsigned int coeffs_id) {
  TargetDistAverages(coeffs_id).setAllValuesToZero();
  TargetDistAverages(coeffs_id).setIterationCounterAndTime(this->getIterationCounter(),this->getTime());
}


void VesBias::checkThatTemperatureIsGiven() {
  if(kbt_==0.0) {
    std::string err_msg = "VES bias " + getLabel() + " of type " + getName() + ": the temperature is needed so you need to give it using the TEMP keyword as the MD engine does not pass it to PLUMED.";
    plumed_merror(err_msg);
  }
}


unsigned int VesBias::getIterationCounter() const {
  unsigned int iteration = 0;
  if(optimizeCoeffs()) {
    iteration = getOptimizerPntr()->getIterationCounter();
  }
  else {
    iteration = getCoeffsPntrs()[0]->getIterationCounter();
  }
  return iteration;
}


void VesBias::linkOptimizer(Optimizer* optimizer_pntr_in) {
  //
  if(optimizer_pntr_==NULL) {
    optimizer_pntr_ = optimizer_pntr_in;
  }
  else {
    std::string err_msg = "VES bias " + getLabel() + " of type " + getName() + " has already been linked with optimizer " + optimizer_pntr_->getLabel() + " of type " + optimizer_pntr_->getName() + ". You cannot link two optimizer to the same VES bias.";
    plumed_merror(err_msg);
  }
  checkThatTemperatureIsGiven();
  optimize_coeffs_ = true;
  filenames_have_iteration_number_ = true;
}


void VesBias::enableHessian(const bool diagonal_hessian) {
  compute_hessian_=true;
  diagonal_hessian_=diagonal_hessian;
  sampled_cross_averages.clear();
  for (unsigned int i=0; i<ncoeffssets_; i++) {
    delete hessian_pntrs_[i];
    std::string label = getCoeffsSetLabelString("hessian",i);
    hessian_pntrs_[i] = new CoeffsMatrix(label,coeffs_pntrs_[i],comm,diagonal_hessian_);
    //
    std::vector<double> cross_aver_sampled_tmp;
    cross_aver_sampled_tmp.assign(hessian_pntrs_[i]->getSize(),0.0);
    sampled_cross_averages.push_back(cross_aver_sampled_tmp);
  }
}


void VesBias::disableHessian() {
  compute_hessian_=false;
  diagonal_hessian_=true;
  sampled_cross_averages.clear();
  for (unsigned int i=0; i<ncoeffssets_; i++) {
    delete hessian_pntrs_[i];
    std::string label = getCoeffsSetLabelString("hessian",i);
    hessian_pntrs_[i] = new CoeffsMatrix(label,coeffs_pntrs_[i],comm,diagonal_hessian_);
    //
    std::vector<double> cross_aver_sampled_tmp;
    cross_aver_sampled_tmp.assign(hessian_pntrs_[i]->getSize(),0.0);
    sampled_cross_averages.push_back(cross_aver_sampled_tmp);
  }
}


std::string VesBias::getCoeffsSetLabelString(const std::string& type, const unsigned int coeffs_id) const {
  std::string label_prefix = getLabel() + ".";
  std::string label_postfix = "";
  if(use_multiple_coeffssets_) {
    Tools::convert(coeffs_id,label_postfix);
    label_postfix = "-" + label_postfix;
  }
  return label_prefix+type+label_postfix;
}


OFile* VesBias::getOFile(const std::string& filepath, const bool multi_sim_single_file, const bool enforce_backup) {
  OFile* ofile_pntr = new OFile();
  std::string fp = filepath;
  ofile_pntr->link(*static_cast<Action*>(this));
  if(enforce_backup) {ofile_pntr->enforceBackup();}
  if(multi_sim_single_file) {
    unsigned int r=0;
    if(comm.Get_rank()==0) {r=multi_sim_comm.Get_rank();}
    comm.Bcast(r,0);
    if(r>0) {fp="/dev/null";}
    ofile_pntr->enforceSuffix("");
  }
  ofile_pntr->open(fp);
  return ofile_pntr;
}


void VesBias::setGridBins(const std::vector<unsigned int>& grid_bins_in) {
  plumed_massert(grid_bins_in.size()==getNumberOfArguments(),"the number of grid bins given doesn't match the number of arguments");
  grid_bins_=grid_bins_in;
}


void VesBias::setGridBins(const unsigned int nbins) {
  std::vector<unsigned int> grid_bins_in(getNumberOfArguments(),nbins);
  grid_bins_=grid_bins_in;
}


void VesBias::setGridMin(const std::vector<double>& grid_min_in) {
  plumed_massert(grid_min_in.size()==getNumberOfArguments(),"the number of lower bounds given for the grid doesn't match the number of arguments");
  grid_min_=grid_min_in;
}


void VesBias::setGridMax(const std::vector<double>& grid_max_in) {
  plumed_massert(grid_max_in.size()==getNumberOfArguments(),"the number of upper bounds given for the grid doesn't match the number of arguments");
  grid_max_=grid_max_in;
}


std::string VesBias::getCurrentOutputFilename(const std::string& base_filename, const std::string& suffix) const {
  std::string filename = base_filename;
  if(suffix.size()>0) {
    filename = FileBase::appendSuffix(filename,"."+suffix);
  }
  if(filenamesIncludeIterationNumber()) {
    filename = FileBase::appendSuffix(filename,"."+getIterationFilenameSuffix());
  }
  return filename;
}


std::string VesBias::getCurrentTargetDistOutputFilename(const std::string& suffix) const {
  std::string filename = targetdist_filename_;
  if(suffix.size()>0) {
    filename = FileBase::appendSuffix(filename,"."+suffix);
  }
  if(filenamesIncludeIterationNumber() && dynamicTargetDistribution()) {
    filename = FileBase::appendSuffix(filename,"."+getIterationFilenameSuffix());
  }
  return filename;
}


std::string VesBias::getIterationFilenameSuffix() const {
  std::string iter_str;
  Tools::convert(getIterationCounter(),iter_str);
  iter_str = "iter-" + iter_str;
  return iter_str;
}


std::string VesBias::getCoeffsSetFilenameSuffix(const unsigned int coeffs_id) const {
  std::string suffix = "";
  if(use_multiple_coeffssets_) {
    Tools::convert(coeffs_id,suffix);
    suffix = coeffs_id_prefix_ + suffix;
  }
  return suffix;
}


void VesBias::setupBiasCutoff(const double bias_cutoff_value, const double fermi_lambda) {
  //
  double fermi_exp_max = 100.0;
  //
  std::string str_bias_cutoff_value; VesTools::convertDbl2Str(bias_cutoff_value,str_bias_cutoff_value);
  std::string str_fermi_lambda; VesTools::convertDbl2Str(fermi_lambda,str_fermi_lambda);
  std::string str_fermi_exp_max; VesTools::convertDbl2Str(fermi_exp_max,str_fermi_exp_max);
  std::string swfunc_keywords = "FERMI R_0=" + str_bias_cutoff_value + " FERMI_LAMBDA=" + str_fermi_lambda + " FERMI_EXP_MAX=" + str_fermi_exp_max;
  //
  if(bias_cutoff_swfunc_pntr_!=NULL) {delete bias_cutoff_swfunc_pntr_;}
  std::string swfunc_errors="";
  bias_cutoff_swfunc_pntr_ = new FermiSwitchingFunction();
  bias_cutoff_swfunc_pntr_->set(swfunc_keywords,swfunc_errors);
  if(swfunc_errors.size()>0) {
    plumed_merror("problem with setting up Fermi switching function: " + swfunc_errors);
  }
  //
  bias_cutoff_value_=bias_cutoff_value;
  bias_cutoff_active_=true;
  enableDynamicTargetDistribution();
}


double VesBias::getBiasCutoffSwitchingFunction(const double bias, double& deriv_factor) const {
  plumed_massert(bias_cutoff_active_,"The bias cutoff is not active so you cannot call this function");
  double arg = -(bias-bias_current_max_value);
  double deriv=0.0;
  double value = bias_cutoff_swfunc_pntr_->calculate(arg,deriv);
  // as FermiSwitchingFunction class has different behavior from normal SwitchingFunction class
  // I was having problems with NaN as it was dividing with zero
  // deriv *= arg;
  deriv_factor = value-bias*deriv;
  return value;
}


bool VesBias::useMultipleWalkers() const {
  bool use_mwalkers_mpi=false;
  if(optimizeCoeffs() && getOptimizerPntr()->useMultipleWalkers()) {
    use_mwalkers_mpi=true;
  }
  return use_mwalkers_mpi;
}


void VesBias::updateReweightFactor() {
  if(calc_reweightfactor_) {
    double value = calculateReweightFactor();
    getPntrToComponent("rct")->set(value);
  }
}


double VesBias::calculateReweightFactor() const {
  plumed_merror(getName()+" with label "+getLabel()+": calculation of the reweight factor c(t) has not been implemented for this type of VES bias");
  return 0.0;
}


}
}
