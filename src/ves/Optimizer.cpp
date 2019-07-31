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

#include "Optimizer.h"
#include "CoeffsVector.h"
#include "CoeffsMatrix.h"
#include "VesBias.h"
#include "VesTools.h"

#include "tools/Exception.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "tools/Communicator.h"
#include "tools/File.h"
#include "tools/FileBase.h"

namespace PLMD {
namespace ves {

Optimizer::Optimizer(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionWithValue(ao),
  description_("Undefined"),
  type_("Undefined"),
  stepsizes_(0),
  current_stepsizes(0),
  fixed_stepsize_(true),
  iter_counter(0),
  use_hessian_(false),
  diagonal_hessian_(true),
  monitor_instantaneous_gradient_(false),
  use_mwalkers_mpi_(false),
  mwalkers_mpi_single_files_(true),
  dynamic_targetdists_(0),
  ustride_targetdist_(0),
  ustride_reweightfactor_(0),
  coeffssetid_prefix_(""),
  coeffs_wstride_(100),
  coeffsOFiles_(0),
  coeffs_output_fmt_(""),
  gradient_wstride_(100),
  gradientOFiles_(0),
  gradient_output_fmt_(""),
  hessian_wstride_(100),
  hessianOFiles_(0),
  hessian_output_fmt_(""),
  targetdist_averages_wstride_(0),
  targetdist_averagesOFiles_(0),
  targetdist_averages_output_fmt_(""),
  nbiases_(0),
  bias_pntrs_(0),
  ncoeffssets_(0),
  coeffs_pntrs_(0),
  aux_coeffs_pntrs_(0),
  gradient_pntrs_(0),
  aver_gradient_pntrs_(0),
  hessian_pntrs_(0),
  coeffs_mask_pntrs_(0),
  targetdist_averages_pntrs_(0),
  identical_coeffs_shape_(true),
  bias_output_active_(false),
  bias_output_stride_(0),
  fes_output_active_(false),
  fes_output_stride_(0),
  fesproj_output_active_(false),
  fesproj_output_stride_(0),
  targetdist_output_active_(false),
  targetdist_output_stride_(0),
  targetdist_proj_output_active_(false),
  targetdist_proj_output_stride_(0),
  isFirstStep(true)
{
  std::vector<std::string> bias_labels(0);
  parseVector("BIAS",bias_labels);
  plumed_massert(bias_labels.size()>0,"problem with BIAS keyword");
  nbiases_ = bias_labels.size();
  //
  std::string error_msg = "";
  bias_pntrs_ = VesTools::getPointersFromLabels<VesBias*>(bias_labels,plumed.getActionSet(),error_msg);
  if(error_msg.size()>0) {plumed_merror("Error in keyword BIAS of "+getName()+": "+error_msg);}

  for(unsigned int i=0; i<bias_pntrs_.size(); i++) {
    bias_pntrs_[i]->linkOptimizer(this);
    //
    std::vector<CoeffsVector*> pntrs_coeffs = bias_pntrs_[i]->getCoeffsPntrs();
    std::vector<CoeffsVector*> pntrs_gradient = bias_pntrs_[i]->getGradientPntrs();
    std::vector<CoeffsVector*> pntrs_targetdist_averages = bias_pntrs_[i]->getTargetDistAveragesPntrs();
    plumed_massert(pntrs_coeffs.size()==pntrs_gradient.size(),"something wrong in the coefficients and gradient passed from VES bias");
    plumed_massert(pntrs_coeffs.size()==pntrs_targetdist_averages.size(),"something wrong in the coefficients and target distribution averages passed from VES bias");
    for(unsigned int k=0; k<pntrs_coeffs.size(); k++) {
      plumed_massert(pntrs_coeffs[k] != NULL,"some coefficient is not linked correctly");
      plumed_massert(pntrs_gradient[k] != NULL,"some gradient is not linked correctly");
      plumed_massert(pntrs_targetdist_averages[k] != NULL,"some target distribution average is not linked correctly");
      pntrs_coeffs[k]->turnOnIterationCounter();
      coeffs_pntrs_.push_back(pntrs_coeffs[k]);
      pntrs_gradient[k]->turnOnIterationCounter();
      gradient_pntrs_.push_back(pntrs_gradient[k]);
      pntrs_targetdist_averages[k]->turnOnIterationCounter();
      targetdist_averages_pntrs_.push_back(pntrs_targetdist_averages[k]);
      //
      CoeffsVector* aux_coeffs_tmp = new CoeffsVector(*pntrs_coeffs[k]);
      std::string aux_label = pntrs_coeffs[k]->getLabel();
      if(aux_label.find("coeffs")!=std::string::npos) {
        aux_label.replace(aux_label.find("coeffs"), std::string("coeffs").length(), "aux_coeffs");
      }
      else {
        aux_label += "_aux";
      }
      aux_coeffs_tmp->setLabels(aux_label);
      aux_coeffs_pntrs_.push_back(aux_coeffs_tmp);
      AuxCoeffs(i).setValues( Coeffs(i) );
    }
  }
  ncoeffssets_ = coeffs_pntrs_.size();
  plumed_massert(aux_coeffs_pntrs_.size()==ncoeffssets_,"problems in linking aux coefficients");
  plumed_massert(gradient_pntrs_.size()==ncoeffssets_,"problems in linking gradients");
  setAllCoeffsSetIterationCounters();


  //
  identical_coeffs_shape_ = true;
  for(unsigned int i=1; i<ncoeffssets_; i++) {
    if(!coeffs_pntrs_[0]->sameShape(*coeffs_pntrs_[i])) {
      identical_coeffs_shape_ = false;
      break;
    }
  }
  //
  if(keywords.exists("STEPSIZE")) {
    plumed_assert(!keywords.exists("INITIAL_STEPSIZE"));
    fixed_stepsize_=true;
    parseMultipleValues("STEPSIZE",stepsizes_);
    setCurrentStepSizes(stepsizes_);
  }
  if(keywords.exists("INITIAL_STEPSIZE")) {
    plumed_assert(!keywords.exists("STEPSIZE"));
    fixed_stepsize_=false;
    parseMultipleValues("INITIAL_STEPSIZE",stepsizes_);
    setCurrentStepSizes(stepsizes_);
  }
  //
  if(ncoeffssets_==1) {
    log.printf("  optimizing VES bias %s with label %s: \n",bias_pntrs_[0]->getName().c_str(),bias_pntrs_[0]->getLabel().c_str());
    log.printf("   KbT: %f\n",bias_pntrs_[0]->getKbT());
    log.printf("  number of coefficients: %zu\n",coeffs_pntrs_[0]->numberOfCoeffs());
    if(stepsizes_.size()>0) {
      if(fixed_stepsize_) {log.printf("  using a constant step size of %f\n",stepsizes_[0]);}
      else {log.printf("  using an initial step size of %f\n",stepsizes_[0]);}
    }
  }
  else {
    log.printf("  optimizing %u coefficient sets from following %u VES biases:\n",ncoeffssets_,nbiases_);
    for(unsigned int i=0; i<nbiases_; i++) {
      log.printf("   %s of type %s (KbT: %f) \n",bias_pntrs_[i]->getLabel().c_str(),bias_pntrs_[i]->getName().c_str(),bias_pntrs_[i]->getKbT());
    }
    size_t tot_ncoeffs = 0;
    for(unsigned int i=0; i<ncoeffssets_; i++) {
      log.printf("  coefficient set %u: \n",i);
      log.printf("   used in bias %s (type %s)\n",coeffs_pntrs_[i]->getPntrToAction()->getLabel().c_str(),coeffs_pntrs_[i]->getPntrToAction()->getName().c_str());
      log.printf("   number of coefficients: %zu\n",coeffs_pntrs_[i]->numberOfCoeffs());
      if(stepsizes_.size()>0) {
        if(fixed_stepsize_) {log.printf("   using a constant step size of %f\n",stepsizes_[i]);}
        else {log.printf("   using an initial step size of %f\n",stepsizes_[i]);}
      }
      tot_ncoeffs += coeffs_pntrs_[i]->numberOfCoeffs();
    }
    log.printf("  total number of coefficients: %zu\n",tot_ncoeffs);
    if(identical_coeffs_shape_) {
      log.printf("  the indices shape is identical for all coefficient sets\n");
    }
    else {
      log.printf("  the indices shape differs between coefficient sets\n");
    }
  }

  //
  if(keywords.exists("FULL_HESSIAN")) {
    bool full_hessian=false;
    parseFlag("FULL_HESSIAN",full_hessian);
    diagonal_hessian_ = !full_hessian;
  }
  //
  bool mw_single_files = false;
  if(keywords.exists("MULTIPLE_WALKERS")) {
    parseFlag("MULTIPLE_WALKERS",use_mwalkers_mpi_);
    if(use_mwalkers_mpi_) {
      mw_single_files=true;
    }
  }

  int numwalkers=1;
  int walker_rank=0;
  if(comm.Get_rank()==0) {
    numwalkers = multi_sim_comm.Get_size();
    walker_rank = multi_sim_comm.Get_rank();
  }
  comm.Bcast(numwalkers,0);
  comm.Bcast(walker_rank,0);
  if(use_mwalkers_mpi_ && numwalkers==1) {
    plumed_merror("using the MULTIPLE_WALKERS keyword does not make sense if running the MD code with a single replica");
  }
  if(use_mwalkers_mpi_) {
    log.printf("  optimization performed using multiple walkers connected via MPI:\n");
    log.printf("   number of walkers: %d\n",numwalkers);
    log.printf("   walker number: %d\n",walker_rank);
    log.printf("   please see and cite ");
    log << plumed.cite("Raiteri, Laio, Gervasio, Micheletti, and Parrinello, J. Phys. Chem. B 110, 3533 (2006)");
    log.printf("\n");
  }

  dynamic_targetdists_.resize(nbiases_,false);
  if(keywords.exists("TARGETDIST_STRIDE")) {
    bool need_stride = false;
    for(unsigned int i=0; i<nbiases_; i++) {
      dynamic_targetdists_[i] = bias_pntrs_[i]->dynamicTargetDistribution();
      if(dynamic_targetdists_[i]) {need_stride = true;}
    }
    parse("TARGETDIST_STRIDE",ustride_targetdist_);
    if(need_stride && ustride_targetdist_==0) {
      plumed_merror("one of the biases has a dynamic target distribution so you need to give stride for updating it by using the TARGETDIST_STRIDE keyword");
    }
    if(!need_stride && ustride_targetdist_!=0) {
      plumed_merror("using the TARGETDIST_STRIDE keyword doesn't make sense as there is no dynamic target distribution to update");
    }
    if(ustride_targetdist_>0) {
      if(nbiases_==1) {
        log.printf("  the target distribution will be updated very %u coefficient iterations\n",ustride_targetdist_);
      }
      else {
        log.printf("  the target distribution will be updated very %u coefficient iterations for the following biases\n   ",ustride_targetdist_);
        for(unsigned int i=0; i<nbiases_; i++) {
          log.printf("%s ",bias_pntrs_[i]->getLabel().c_str());
        }
        log.printf("\n");
      }
      log.printf("  See and cite ");
      log << plumed.cite("Valsson and Parrinello, J. Chem. Theory Comput. 11, 1996-2002 (2015)");
      log.printf("\n");
    }
  }

  if(keywords.exists("REWEIGHT_FACTOR_STRIDE")) {
    bool reweightfactor_calculated = false;
    for(unsigned int i=0; i<nbiases_; i++) {
      reweightfactor_calculated = bias_pntrs_[i]->isReweightFactorCalculated();
    }
    parse("REWEIGHT_FACTOR_STRIDE",ustride_reweightfactor_);
    if(ustride_reweightfactor_==0 && reweightfactor_calculated) {
      plumed_merror("the calculation of the reweight factor is enabled, You need to use the REWEIGHT_FACTOR_STRIDE keyword to specfiy how often it should be updated.");
    }
    if(ustride_reweightfactor_>0) {
      if(!reweightfactor_calculated) {
        plumed_merror("In order to use the REWEIGHT_FACTOR_STRIDE keyword you need to enable the calculation of the reweight factor in the VES bias by using the CALC_REWEIGHT_FACTOR flag.");
      }
      log.printf("  the reweight factor c(t) will be updated very %u coefficient iterations\n",ustride_reweightfactor_);
    }
  }

  if(keywords.exists("MONITOR_INSTANTANEOUS_GRADIENT")) {
    parseFlag("MONITOR_INSTANTANEOUS_GRADIENT",monitor_instantaneous_gradient_);
  }

  if(keywords.exists("MONITOR_AVERAGE_GRADIENT")) {
    bool monitor_aver_gradient = false;
    parseFlag("MONITOR_AVERAGE_GRADIENT",monitor_aver_gradient);
    if(monitor_aver_gradient) {
      unsigned int averaging_exp_decay=0;
      parse("MONITOR_AVERAGES_GRADIENT_EXP_DECAY",averaging_exp_decay);
      aver_gradient_pntrs_.clear();
      for(unsigned int i=0; i<ncoeffssets_; i++) {
        CoeffsVector* aver_gradient_tmp = new CoeffsVector(*gradient_pntrs_[i]);
        aver_gradient_tmp->clear();
        std::string aver_grad_label = aver_gradient_tmp->getLabel();
        if(aver_grad_label.find("gradient")!=std::string::npos) {
          aver_grad_label.replace(aver_grad_label.find("gradient"), std::string("gradient").length(), "aver_gradient");
        }
        else {
          aver_grad_label += "_aver";
        }
        aver_gradient_tmp->setLabels(aver_grad_label);
        if(averaging_exp_decay>0) {
          aver_gradient_tmp->setupExponentiallyDecayingAveraging(averaging_exp_decay);
        }
        aver_gradient_pntrs_.push_back(aver_gradient_tmp);
      }
    }
  }


  if(ncoeffssets_>1) {
    coeffssetid_prefix_="c-";
    if(keywords.exists("COEFFS_SET_ID_PREFIX")) {
      parse("COEFFS_SET_ID_PREFIX",coeffssetid_prefix_);
    }
  }
  else {
    coeffssetid_prefix_="";
    if(keywords.exists("COEFFS_SET_ID_PREFIX")) {
      parse("COEFFS_SET_ID_PREFIX",coeffssetid_prefix_);
    }
    if(coeffssetid_prefix_.size()>0) {
      plumed_merror("COEFFS_SET_ID_PREFIX should only be given if optimizing multiple coefficient sets");
    }
  }

  if(keywords.exists("INITIAL_COEFFS")) {
    std::vector<std::string> initial_coeffs_fnames;
    parseFilenames("INITIAL_COEFFS",initial_coeffs_fnames);
    if(initial_coeffs_fnames.size()>0) {
      readCoeffsFromFiles(initial_coeffs_fnames,false);
      comm.Barrier();
      if(comm.Get_rank()==0 && use_mwalkers_mpi_) {
        multi_sim_comm.Barrier();
      }
      setAllCoeffsSetIterationCounters();
    }
  }
  //

  std::vector<std::string> coeffs_fnames;
  if(keywords.exists("COEFFS_FILE")) {
    parseFilenames("COEFFS_FILE",coeffs_fnames,"coeffs.data");
    bool start_opt_afresh=false;
    if(keywords.exists("START_OPTIMIZATION_AFRESH")) {
      parseFlag("START_OPTIMIZATION_AFRESH",start_opt_afresh);
      if(start_opt_afresh && !getRestart()) {
        plumed_merror("the START_OPTIMIZATION_AFRESH keyword should only be used when a restart has been triggered by the RESTART keyword or the MD code");
      }
    }
    if(getRestart()) {
      for(unsigned int i=0; i<coeffs_fnames.size(); i++) {
        IFile ifile;
        ifile.link(*this);
        if(use_mwalkers_mpi_) {ifile.enforceSuffix("");}
        bool file_exist = ifile.FileExist(coeffs_fnames[i]);
        if(!file_exist) {
          std::string fname = FileBase::appendSuffix(coeffs_fnames[i],ifile.getSuffix());
          plumed_merror("Cannot find coefficient file " + fname + " when trying to restart an optimzation. If you don't want to restart the optimzation please remove the RESTART keyword or use the RESTART=NO within the "+getName()+" action to locally disable the restart.");
        }
      }
      readCoeffsFromFiles(coeffs_fnames,true);
      comm.Barrier();
      if(comm.Get_rank()==0 && use_mwalkers_mpi_) {
        multi_sim_comm.Barrier();
      }
      unsigned int iter_opt_tmp = coeffs_pntrs_[0]->getIterationCounter();
      for(unsigned int i=1; i<ncoeffssets_; i++) {
        plumed_massert(coeffs_pntrs_[i]->getIterationCounter()==iter_opt_tmp,"the iteraton counter should be the same for all files when restarting from previous coefficient files\n");
      }
      if(start_opt_afresh) {
        setIterationCounter(0);
        log.printf("  Optimization started afresh at iteration %u\n",getIterationCounter());
        for(unsigned int i=0; i<ncoeffssets_; i++) {
          AuxCoeffs(i).setValues( Coeffs(i) );
        }
      }
      else {
        setIterationCounter(coeffs_pntrs_[0]->getIterationCounter());
        log.printf("  Optimization restarted at iteration %u\n",getIterationCounter());
      }
      setAllCoeffsSetIterationCounters();
    }

    std::string coeffs_wstride_tmpstr="";
    parse("COEFFS_OUTPUT",coeffs_wstride_tmpstr);
    if(coeffs_wstride_tmpstr!="OFF" && coeffs_wstride_tmpstr.size()>0) {
      Tools::convert(coeffs_wstride_tmpstr,coeffs_wstride_);
    }
    if(coeffs_wstride_tmpstr=="OFF") {
      coeffs_fnames.clear();
    }
    setupOFiles(coeffs_fnames,coeffsOFiles_,mw_single_files);
    parse("COEFFS_FMT",coeffs_output_fmt_);
    if(coeffs_output_fmt_.size()>0) {
      for(unsigned int i=0; i<ncoeffssets_; i++) {
        coeffs_pntrs_[i]->setOutputFmt(coeffs_output_fmt_);
      }
    }
    if(!getRestart()) {
      for(unsigned int i=0; i<coeffsOFiles_.size(); i++) {
        coeffs_pntrs_[i]->writeToFile(*coeffsOFiles_[i],aux_coeffs_pntrs_[i],false);
      }
    }
    if(coeffs_fnames.size()>0) {
      if(ncoeffssets_==1) {
        log.printf("  Coefficients will be written out to file %s every %u iterations\n",coeffsOFiles_[0]->getPath().c_str(),coeffs_wstride_);
      }
      else {
        log.printf("  Coefficients will be written out to the following files every %u iterations:\n",coeffs_wstride_);
        for(unsigned int i=0; i<coeffs_fnames.size(); i++) {
          log.printf("   coefficient set %u: %s\n",i,coeffsOFiles_[i]->getPath().c_str());
        }
      }
    }
    else {
      log.printf("  Output of coefficients to file has been disabled\n");
    }
  }

  std::vector<std::string> gradient_fnames;
  if(keywords.exists("GRADIENT_FILE")) {
    parseFilenames("GRADIENT_FILE",gradient_fnames);
    parse("GRADIENT_OUTPUT",gradient_wstride_);

    if(coeffs_fnames.size()>0) {
      for(unsigned int i=0; i<gradient_fnames.size(); i++) {
        plumed_massert(gradient_fnames[i]!=coeffs_fnames[i],"COEFFS_FILE and GRADIENT_FILE cannot be the same");
      }
    }
    setupOFiles(gradient_fnames,gradientOFiles_,mw_single_files);
    parse("GRADIENT_FMT",gradient_output_fmt_);
    if(gradient_output_fmt_.size()>0) {
      for(unsigned int i=0; i<ncoeffssets_; i++) {
        gradient_pntrs_[i]->setOutputFmt(gradient_output_fmt_);
      }
    }

    if(gradient_fnames.size()>0) {
      if(ncoeffssets_==1) {
        log.printf("  Gradient will be written out to file %s every %u iterations\n",gradientOFiles_[0]->getPath().c_str(),gradient_wstride_);
      }
      else {
        log.printf("  Gradient will be written out to the following files every %u iterations:\n",gradient_wstride_);
        for(unsigned int i=0; i<gradient_fnames.size(); i++) {
          log.printf("   coefficient set %u: %s\n",i,gradientOFiles_[i]->getPath().c_str());
        }
      }
    }
  }

  std::vector<std::string> hessian_fnames;
  if(keywords.exists("HESSIAN_FILE")) {
    parseFilenames("HESSIAN_FILE",hessian_fnames);
    parse("HESSIAN_OUTPUT",hessian_wstride_);

    if(coeffs_fnames.size()>0) {
      for(unsigned int i=0; i<hessian_fnames.size(); i++) {
        plumed_massert(hessian_fnames[i]!=coeffs_fnames[i],"COEFFS_FILE and HESSIAN_FILE cannot be the same");
      }
    }
    if(gradient_fnames.size()>0) {
      for(unsigned int i=0; i<hessian_fnames.size(); i++) {
        plumed_massert(hessian_fnames[i]!=gradient_fnames[i],"GRADIENT_FILE and HESSIAN_FILE cannot be the same");
      }
    }
    setupOFiles(hessian_fnames,hessianOFiles_,mw_single_files);
    parse("HESSIAN_FMT",hessian_output_fmt_);

    if(hessian_fnames.size()>0) {
      if(ncoeffssets_==1) {
        log.printf("  Hessian will be written out to file %s every %u iterations\n",hessianOFiles_[0]->getPath().c_str(),hessian_wstride_);
      }
      else {
        log.printf("  Hessian will be written out to the following files every %u iterations:\n",hessian_wstride_);
        for(unsigned int i=0; i<hessian_fnames.size(); i++) {
          log.printf("   coefficient set %u: %s\n",i,hessianOFiles_[i]->getPath().c_str());
        }
      }
    }
  }


  //
  if(keywords.exists("MASK_FILE")) {
    std::vector<std::string> mask_fnames_in;
    parseVector("MASK_FILE",mask_fnames_in);
    if(mask_fnames_in.size()==1 && ncoeffssets_>1) {
      if(identical_coeffs_shape_) {mask_fnames_in.resize(ncoeffssets_,mask_fnames_in[0]);}
      else {plumed_merror("the coefficients indices shape differs between biases so you need to give a separate file for each coefficient set\n");}
    }
    if(mask_fnames_in.size()>0 && mask_fnames_in.size()!=ncoeffssets_) {
      plumed_merror("Error in MASK_FILE keyword: either give one value for all biases or a separate value for each coefficient set");
    }

    coeffs_mask_pntrs_.resize(ncoeffssets_);
    for(unsigned int i=0; i<ncoeffssets_; i++) {
      coeffs_mask_pntrs_[i] = new CoeffsVector(*coeffs_pntrs_[i]);
      coeffs_mask_pntrs_[i]->setLabels("mask");
      coeffs_mask_pntrs_[i]->setValues(1.0);
      coeffs_mask_pntrs_[i]->setOutputFmt("%f");
    }

    if(mask_fnames_in.size()>0) {
      if(ncoeffssets_==1) {
        size_t nread = coeffs_mask_pntrs_[0]->readFromFile(mask_fnames_in[0],true,true);
        log.printf("  read %zu values from mask file %s\n",nread,mask_fnames_in[0].c_str());
        size_t ndeactived = coeffs_mask_pntrs_[0]->countValues(0.0);
        log.printf("  deactived optimization of %zu coefficients\n",ndeactived);
      }
      else {
        for(unsigned int i=0; i<ncoeffssets_; i++) {
          size_t nread = coeffs_mask_pntrs_[i]->readFromFile(mask_fnames_in[i],true,true);
          log.printf("  mask for coefficient set %u:\n",i);
          log.printf("   read %zu values from file %s\n",nread,mask_fnames_in[i].c_str());
          size_t ndeactived = coeffs_mask_pntrs_[0]->countValues(0.0);
          log.printf("   deactived optimization of %zu coefficients\n",ndeactived);
        }
      }
    }

    std::vector<std::string> mask_fnames_out;
    parseFilenames("OUTPUT_MASK_FILE",mask_fnames_out);

    for(unsigned int i=0; i<mask_fnames_out.size(); i++) {
      if(mask_fnames_in.size()>0) {
        plumed_massert(mask_fnames_out[i]!=mask_fnames_in[i],"MASK_FILE and OUTPUT_MASK_FILE cannot be the same");
      }
      OFile maskOFile;
      maskOFile.link(*this);
      maskOFile.enforceBackup();
      if(use_mwalkers_mpi_ && mwalkers_mpi_single_files_) {
        unsigned int r=0;
        if(comm.Get_rank()==0) {r=multi_sim_comm.Get_rank();}
        comm.Bcast(r,0);
        if(r>0) {mask_fnames_out[i]="/dev/null";}
        maskOFile.enforceSuffix("");
      }
      maskOFile.open(mask_fnames_out[i]);
      coeffs_mask_pntrs_[i]->writeToFile(maskOFile,true);
      maskOFile.close();
    }
  }

  if(getRestart() && ustride_targetdist_>0) {
    for(unsigned int i=0; i<nbiases_; i++) {
      if(dynamic_targetdists_[i]) {
        bias_pntrs_[i]->restartTargetDistributions();
      }
    }
  }


  std::vector<std::string> targetdist_averages_fnames;
  if(keywords.exists("TARGETDIST_AVERAGES_FILE")) {
    parseFilenames("TARGETDIST_AVERAGES_FILE",targetdist_averages_fnames,"targetdist-averages.data");
    parse("TARGETDIST_AVERAGES_OUTPUT",targetdist_averages_wstride_);

    if(coeffs_fnames.size()>0) {
      for(unsigned int i=0; i<targetdist_averages_fnames.size(); i++) {
        plumed_massert(targetdist_averages_fnames[i]!=coeffs_fnames[i],"COEFFS_FILE and TARGETDIST_AVERAGES_FILE cannot be the same");
      }
    }
    if(gradient_fnames.size()>0) {
      for(unsigned int i=0; i<targetdist_averages_fnames.size(); i++) {
        plumed_massert(targetdist_averages_fnames[i]!=gradient_fnames[i],"GRADIENT_FILE and TARGETDIST_AVERAGES_FILE cannot be the same");
      }
    }
    if(hessian_fnames.size()>0) {
      for(unsigned int i=0; i<targetdist_averages_fnames.size(); i++) {
        plumed_massert(targetdist_averages_fnames[i]!=hessian_fnames[i],"HESSIAN_FILE and TARGETDIST_AVERAGES_FILE cannot be the same");
      }
    }
    setupOFiles(targetdist_averages_fnames,targetdist_averagesOFiles_,mw_single_files);
    parse("TARGETDIST_AVERAGES_FMT",targetdist_averages_output_fmt_);
    if(targetdist_averages_output_fmt_.size()>0) {
      for(unsigned int i=0; i<ncoeffssets_; i++) {
        targetdist_averages_pntrs_[i]->setOutputFmt(targetdist_averages_output_fmt_);
      }
    }

    for(unsigned int i=0; i<targetdist_averagesOFiles_.size(); i++) {
      targetdist_averages_pntrs_[i]->writeToFile(*targetdist_averagesOFiles_[i]);
    }

    if(targetdist_averages_wstride_==0) {
      for(unsigned int i=0; i<targetdist_averagesOFiles_.size(); i++) {
        targetdist_averagesOFiles_[i]->close();
        delete targetdist_averagesOFiles_[i];
      }
      targetdist_averagesOFiles_.clear();
    }

    if(targetdist_averages_fnames.size()>0 && targetdist_averages_wstride_ > 0) {
      if(ncoeffssets_==1) {
        log.printf("  Target distribution averages will be written out to file %s every %u iterations\n",targetdist_averagesOFiles_[0]->getPath().c_str(),targetdist_averages_wstride_);
      }
      else {
        log.printf("  Target distribution averages will be written out to the following files every %u iterations:\n",targetdist_averages_wstride_);
        for(unsigned int i=0; i<targetdist_averages_fnames.size(); i++) {
          log.printf("   coefficient set %u: %s\n",i,targetdist_averagesOFiles_[i]->getPath().c_str());
        }
      }
    }
  }


  if(keywords.exists("BIAS_OUTPUT")) {
    parse("BIAS_OUTPUT",bias_output_stride_);
    if(bias_output_stride_>0) {
      bias_output_active_=true;
      for(unsigned int i=0; i<nbiases_; i++) {
        bias_pntrs_[i]->enableBiasFileOutput();
        bias_pntrs_[i]->setupBiasFileOutput();
        bias_pntrs_[i]->writeBiasToFile();
      }
    }
    else {
      bias_output_active_=false;
      bias_output_stride_=1000;
    }
  }

  if(keywords.exists("FES_OUTPUT")) {
    parse("FES_OUTPUT",fes_output_stride_);
    if(fes_output_stride_>0) {
      fes_output_active_=true;
      for(unsigned int i=0; i<nbiases_; i++) {
        bias_pntrs_[i]->enableFesFileOutput();
        bias_pntrs_[i]->setupFesFileOutput();
        bias_pntrs_[i]->writeFesToFile();
      }
    }
    else {
      fes_output_active_=false;
      fes_output_stride_=1000;
    }
  }

  if(keywords.exists("FES_PROJ_OUTPUT")) {
    parse("FES_PROJ_OUTPUT",fesproj_output_stride_);
    if(fesproj_output_stride_>0) {
      fesproj_output_active_=true;
      for(unsigned int i=0; i<nbiases_; i++) {
        bias_pntrs_[i]->enableFesProjFileOutput();
        bias_pntrs_[i]->setupFesProjFileOutput();
        bias_pntrs_[i]->writeFesProjToFile();
      }
    }
    else {
      fesproj_output_active_=false;
      fesproj_output_stride_=1000;
    }
  }

  for(unsigned int i=0; i<nbiases_; i++) {
    if(!dynamic_targetdists_[i] && bias_pntrs_[i]->isStaticTargetDistFileOutputActive()) {
      bias_pntrs_[i]->setupTargetDistFileOutput();
      bias_pntrs_[i]->writeTargetDistToFile();
      bias_pntrs_[i]->setupTargetDistProjFileOutput();
      bias_pntrs_[i]->writeTargetDistProjToFile();
    }
  }

  if(keywords.exists("TARGETDIST_OUTPUT")) {
    parse("TARGETDIST_OUTPUT",targetdist_output_stride_);
    if(targetdist_output_stride_>0) {
      if(ustride_targetdist_==0) {
        plumed_merror("it doesn't make sense to use the TARGETDIST_OUTPUT keyword if you don't have a target distribution that needs to be updated");
      }
      if(targetdist_output_stride_%ustride_targetdist_!=0) {
        plumed_merror("the value given in TARGETDIST_OUTPUT doesn't make sense, it should be multiple of TARGETDIST_STRIDE");
      }

      targetdist_output_active_=true;
      for(unsigned int i=0; i<nbiases_; i++) {
        if(dynamic_targetdists_[i]) {
          bias_pntrs_[i]->enableDynamicTargetDistFileOutput();
          bias_pntrs_[i]->setupTargetDistFileOutput();
          bias_pntrs_[i]->writeTargetDistToFile();
        }
      }
    }
    else {
      targetdist_output_active_=false;
      targetdist_output_stride_=1000;
    }
  }

  if(keywords.exists("TARGETDIST_PROJ_OUTPUT")) {
    parse("TARGETDIST_PROJ_OUTPUT",targetdist_proj_output_stride_);
    if(targetdist_proj_output_stride_>0) {
      if(ustride_targetdist_==0) {
        plumed_merror("it doesn't make sense to use the TARGETDIST_PROJ_OUTPUT keyword if you don't have a target distribution that needs to be updated");
      }
      if(targetdist_proj_output_stride_%ustride_targetdist_!=0) {
        plumed_merror("the value given in TARGETDIST_PROJ_OUTPUT doesn't make sense, it should be multiple of TARGETDIST_STRIDE");
      }

      targetdist_proj_output_active_=true;
      for(unsigned int i=0; i<nbiases_; i++) {
        if(dynamic_targetdists_[i]) {
          bias_pntrs_[i]->enableDynamicTargetDistFileOutput();
          bias_pntrs_[i]->setupTargetDistProjFileOutput();
          bias_pntrs_[i]->writeTargetDistProjToFile();
        }
      }
    }
    else {
      targetdist_proj_output_active_=false;
      targetdist_proj_output_stride_=1000;
    }
  }

  if(ncoeffssets_==1) {
    log.printf("  Output Components:\n");
    log.printf(" ");
    if(monitor_instantaneous_gradient_) {
      addComponent("gradrms"); componentIsNotPeriodic("gradrms");
      log.printf(" ");
      addComponent("gradmax"); componentIsNotPeriodic("gradmax");
    }
    if(aver_gradient_pntrs_.size()>0) {
      log.printf(" ");
      addComponent("avergradrms"); componentIsNotPeriodic("avergradrms");
      log.printf(" ");
      addComponent("avergradmax"); componentIsNotPeriodic("avergradmax");
    }
    if(!fixed_stepsize_) {
      log.printf(" ");
      addComponent("stepsize"); componentIsNotPeriodic("stepsize");
      getPntrToComponent("stepsize")->set( getCurrentStepSize(0) );
    }
  }
  else {
    for(unsigned int i=0; i<ncoeffssets_; i++) {
      log.printf("  Output Components for coefficient set %u:\n",i);
      std::string is=""; Tools::convert(i,is); is = "-" + coeffssetid_prefix_ + is;
      log.printf(" ");
      if(monitor_instantaneous_gradient_) {
        addComponent("gradrms"+is); componentIsNotPeriodic("gradrms"+is);
        log.printf(" ");
        addComponent("gradmax"+is); componentIsNotPeriodic("gradmax"+is);
      }
      if(aver_gradient_pntrs_.size()>0) {
        log.printf(" ");
        addComponent("avergradrms"+is); componentIsNotPeriodic("avergradrms"+is);
        log.printf(" ");
        addComponent("avergradmax"+is); componentIsNotPeriodic("avergradmax"+is);
      }
      if(!fixed_stepsize_) {
        log.printf(" ");
        addComponent("stepsize"+is); componentIsNotPeriodic("stepsize"+is);
        getPntrToComponent("stepsize"+is)->set( getCurrentStepSize(i) );
      }
    }
  }

}


Optimizer::~Optimizer() {
  //
  if(!isTargetDistOutputActive()) {
    for(unsigned int i=0; i<nbiases_; i++) {
      if(dynamic_targetdists_[i]) {
        bias_pntrs_[i]->enableDynamicTargetDistFileOutput();
        bias_pntrs_[i]->setupTargetDistFileOutput();
        bias_pntrs_[i]->writeTargetDistToFile();
      }
    }
  }
  else if(isTargetDistOutputActive() && getIterationCounter()%getTargetDistOutputStride()!=0) {
    writeTargetDistOutputFiles();
  }
  //
  for(unsigned int i=0; i<aux_coeffs_pntrs_.size(); i++) {
    delete aux_coeffs_pntrs_[i];
  }
  aux_coeffs_pntrs_.clear();
  //
  for(unsigned int i=0; i<aver_gradient_pntrs_.size(); i++) {
    delete aver_gradient_pntrs_[i];
  }
  aver_gradient_pntrs_.clear();
  //
  for(unsigned int i=0; i<coeffs_mask_pntrs_.size(); i++) {
    delete coeffs_mask_pntrs_[i];
  }
  coeffs_mask_pntrs_.clear();
  //
  for(unsigned int i=0; i<coeffsOFiles_.size(); i++) {
    coeffsOFiles_[i]->close();
    delete coeffsOFiles_[i];
  }
  coeffsOFiles_.clear();
  for(unsigned int i=0; i<gradientOFiles_.size(); i++) {
    gradientOFiles_[i]->close();
    delete gradientOFiles_[i];
  }
  gradientOFiles_.clear();
  for(unsigned int i=0; i<hessianOFiles_.size(); i++) {
    hessianOFiles_[i]->close();
    delete hessianOFiles_[i];
  }
  hessianOFiles_.clear();
  for(unsigned int i=0; i<targetdist_averagesOFiles_.size(); i++) {
    targetdist_averagesOFiles_[i]->close();
    delete targetdist_averagesOFiles_[i];
  }
  targetdist_averagesOFiles_.clear();
}


void Optimizer::registerKeywords( Keywords& keys ) {
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  //
  keys.remove("NUMERICAL_DERIVATIVES");
  // Default always active keywords
  keys.add("compulsory","BIAS","the label of the VES bias to be optimized");
  keys.add("compulsory","STRIDE","the frequency of updating the coefficients given in the number of MD steps.");
  keys.add("compulsory","COEFFS_FILE","coeffs.data","the name of output file for the coefficients");
  keys.add("compulsory","COEFFS_OUTPUT","100","how often the coefficients should be written to file. This parameter is given as the number of iterations.");
  keys.add("optional","COEFFS_FMT","specify format for coefficient file(s) (useful for decrease the number of digits in regtests)");
  keys.add("optional","COEFFS_SET_ID_PREFIX","suffix to add to the filename given in FILE to identify the bias, should only be given if a single filename is given in FILE when optimizing multiple biases.");
  //
  keys.add("optional","INITIAL_COEFFS","the name(s) of file(s) with the initial coefficients");
  // Hidden keywords to output the gradient to a file.
  keys.add("hidden","GRADIENT_FILE","the name of output file for the gradient");
  keys.add("hidden","GRADIENT_OUTPUT","how often the gradient should be written to file. This parameter is given as the number of bias iterations. It is by default 100 if GRADIENT_FILE is specficed");
  keys.add("hidden","GRADIENT_FMT","specify format for gradient file(s) (useful for decrease the number of digits in regtests)");
  // Either use a fixed stepsize (useFixedStepSizeKeywords) or changing stepsize (useDynamicsStepSizeKeywords)
  keys.reserve("compulsory","STEPSIZE","the step size used for the optimization");
  keys.reserve("compulsory","INITIAL_STEPSIZE","the initial step size used for the optimization");
  // Keywords related to the Hessian, actived with the useHessianKeywords function
  keys.reserveFlag("FULL_HESSIAN",false,"if the full Hessian matrix should be used for the optimization, otherwise only the diagonal part of the Hessian is used");
  keys.reserve("hidden","HESSIAN_FILE","the name of output file for the Hessian");
  keys.reserve("hidden","HESSIAN_OUTPUT","how often the Hessian should be written to file. This parameter is given as the number of bias iterations. It is by default 100 if HESSIAN_FILE is specficed");
  keys.reserve("hidden","HESSIAN_FMT","specify format for hessian file(s) (useful for decrease the number of digits in regtests)");
  // Keywords related to the multiple walkers, actived with the useMultipleWalkersKeywords function
  keys.reserveFlag("MULTIPLE_WALKERS",false,"if optimization is to be performed using multiple walkers connected via MPI");
  // Keywords related to the mask file, actived with the useMaskKeywords function
  keys.reserve("optional","MASK_FILE","read in a mask file which allows one to employ different step sizes for different coefficients and/or deactivate the optimization of certain coefficients (by putting values of 0.0). One can write out the resulting mask by using the OUTPUT_MASK_FILE keyword.");
  keys.reserve("optional","OUTPUT_MASK_FILE","Name of the file to write out the mask resulting from using the MASK_FILE keyword. Can also be used to generate a template mask file.");
  //
  keys.reserveFlag("START_OPTIMIZATION_AFRESH",false,"if the iterations should be started afresh when a restart has been triggered by the RESTART keyword or the MD code.");
  //
  keys.addFlag("MONITOR_INSTANTANEOUS_GRADIENT",false,"if quantities related to the instantaneous gradient should be outputted.");
  //
  keys.reserveFlag("MONITOR_AVERAGE_GRADIENT",false,"if the averaged gradient should be monitored and quantities related to it should be outputted.");
  keys.reserve("optional","MONITOR_AVERAGES_GRADIENT_EXP_DECAY","use an exponentially decaying averaging with a given time constant when monitoring the averaged gradient");
  //
  keys.reserve("optional","TARGETDIST_STRIDE","stride for updating a target distribution that is iteratively updated during the optimization. Note that the value is given in terms of coefficient iterations.");
  keys.reserve("optional","TARGETDIST_OUTPUT","how often the dynamic target distribution(s) should be written out to file. Note that the value is given in terms of coefficient iterations.");
  keys.reserve("optional","TARGETDIST_PROJ_OUTPUT","how often the projections of the dynamic target distribution(s) should be written out to file. Note that the value is given in terms of coefficient iterations.");
  //
  keys.add("optional","TARGETDIST_AVERAGES_FILE","the name of output file for the target distribution averages. By default it is targetdist-averages.data.");
  keys.add("optional","TARGETDIST_AVERAGES_OUTPUT","how often the target distribution averages should be written out to file. Note that the value is given in terms of coefficient iterations. If no value is given are the averages only written at the beginning of the optimization");
  keys.add("hidden","TARGETDIST_AVERAGES_FMT","specify format for target distribution averages file(s) (useful for decrease the number of digits in regtests)");
  //
  keys.add("optional","BIAS_OUTPUT","how often the bias(es) should be written out to file. Note that the value is given in terms of coefficient iterations.");
  keys.add("optional","FES_OUTPUT","how often the FES(s) should be written out to file. Note that the value is given in terms of coefficient iterations.");
  keys.add("optional","FES_PROJ_OUTPUT","how often the projections of the FES(s) should be written out to file. Note that the value is given in terms of coefficient iterations.");
  //
  keys.reserve("optional","REWEIGHT_FACTOR_STRIDE","stride for updating the reweighting factor c(t). Note that the value is given in terms of coefficient iterations.");
  //
  keys.use("RESTART");
  //
  keys.use("UPDATE_FROM");
  keys.use("UPDATE_UNTIL");
  // Components that are always active
  keys.addOutputComponent("gradrms","MONITOR_INSTANTANEOUS_GRADIENT","the root mean square value of the coefficient gradient. For multiple biases this component is labeled using the number of the bias as gradrms-#.");
  keys.addOutputComponent("gradmax","MONITOR_INSTANTANEOUS_GRADIENT","the largest absolute value of the coefficient gradient. For multiple biases this component is labeled using the number of the bias as gradmax-#.");
  // keys.addOutputComponent("gradmaxidx","default","the index of the maximum absolute value of the gradient");

}


void Optimizer::useHessianKeywords(Keywords& keys) {
  // keys.use("FULL_HESSIAN");
  keys.use("HESSIAN_FILE");
  keys.use("HESSIAN_OUTPUT");
  keys.use("HESSIAN_FMT");
}


void Optimizer::useMultipleWalkersKeywords(Keywords& keys) {
  keys.use("MULTIPLE_WALKERS");
}


void Optimizer::useFixedStepSizeKeywords(Keywords& keys) {
  keys.use("STEPSIZE");
}


void Optimizer::useDynamicStepSizeKeywords(Keywords& keys) {
  keys.use("INITIAL_STEPSIZE");
  keys.addOutputComponent("stepsize","default","the current value of step size used to update the coefficients. For multiple biases this component is labeled using the number of the bias as stepsize-#.");
}


void Optimizer::useMaskKeywords(Keywords& keys) {
  keys.use("MASK_FILE");
  keys.use("OUTPUT_MASK_FILE");
}


void Optimizer::useRestartKeywords(Keywords& keys) {
  keys.use("START_OPTIMIZATION_AFRESH");
}


void Optimizer::useMonitorAverageGradientKeywords(Keywords& keys) {
  keys.use("MONITOR_AVERAGE_GRADIENT");
  keys.use("MONITOR_AVERAGES_GRADIENT_EXP_DECAY");
  keys.addOutputComponent("avergradrms","MONITOR_AVERAGE_GRADIENT","the root mean square value of the averaged coefficient gradient. For multiple biases this component is labeled using the number of the bias as gradrms-#.");
  keys.addOutputComponent("avergradmax","MONITOR_AVERAGE_GRADIENT","the largest absolute value of the averaged coefficient gradient. For multiple biases this component is labeled using the number of the bias as gradmax-#.");
}


void Optimizer::useDynamicTargetDistributionKeywords(Keywords& keys) {
  keys.use("TARGETDIST_STRIDE");
  keys.use("TARGETDIST_OUTPUT");
  keys.use("TARGETDIST_PROJ_OUTPUT");
}


void Optimizer::useReweightFactorKeywords(Keywords& keys) {
  keys.use("REWEIGHT_FACTOR_STRIDE");
}


void Optimizer::turnOnHessian() {
  plumed_massert(hessian_pntrs_.size()==0,"turnOnHessian() should only be run during initialization");
  use_hessian_=true;
  hessian_pntrs_.clear();
  for(unsigned int i=0; i<nbiases_; i++) {
    std::vector<CoeffsMatrix*> pntrs_hessian = enableHessian(bias_pntrs_[i],diagonal_hessian_);
    for(unsigned int k=0; k<pntrs_hessian.size(); k++) {
      pntrs_hessian[k]->turnOnIterationCounter();
      pntrs_hessian[k]->setIterationCounterAndTime(getIterationCounter(),getTime());
      hessian_pntrs_.push_back(pntrs_hessian[k]);
    }
  }
  plumed_massert(hessian_pntrs_.size()==ncoeffssets_,"problems in linking Hessians");
  if(diagonal_hessian_) {
    log.printf("  Optimization performed using diagonal Hessian matrix\n");
  }
  else {
    log.printf("  Optimization performed using full Hessian matrix\n");
  }
  //
  if(hessian_output_fmt_.size()>0) {
    for(unsigned int i=0; i<ncoeffssets_; i++) {
      hessian_pntrs_[i]->setOutputFmt(hessian_output_fmt_);
    }
  }

}


void Optimizer::turnOffHessian() {
  use_hessian_=false;
  for(unsigned int i=0; i<nbiases_; i++) {
    bias_pntrs_[i]->disableHessian();
  }
  hessian_pntrs_.clear();
  for(unsigned int i=0; i<hessianOFiles_.size(); i++) {
    hessianOFiles_[i]->close();
    delete hessianOFiles_[i];
  }
  hessianOFiles_.clear();
}


std::vector<CoeffsMatrix*> Optimizer::enableHessian(VesBias* bias_pntr_in, const bool diagonal_hessian) {
  plumed_massert(use_hessian_,"the Hessian should not be used");
  bias_pntr_in->enableHessian(diagonal_hessian);
  std::vector<CoeffsMatrix*> hessian_pntrs_out = bias_pntr_in->getHessianPntrs();
  for(unsigned int k=0; k<hessian_pntrs_out.size(); k++) {
    plumed_massert(hessian_pntrs_out[k] != NULL,"Hessian is needed but not linked correctly");
  }
  return hessian_pntrs_out;
}


// CoeffsMatrix* Optimizer::switchToDiagonalHessian(VesBias* bias_pntr_in) {
//   plumed_massert(use_hessian_,"it does not make sense to switch to diagonal Hessian if it Hessian is not used");
//   diagonal_hessian_=true;
//   bias_pntr_in->enableHessian(diagonal_hessian_);
//   CoeffsMatrix* hessian_pntr_out = bias_pntr_in->getHessianPntr();
//   plumed_massert(hessian_pntr_out != NULL,"Hessian is needed but not linked correctly");
//   //
//   log.printf("  %s (with label %s): switching to a diagonal Hessian for VES bias %s (with label %s) at time  %f\n",getName().c_str(),getLabel().c_str(),bias_pntr_in->getName().c_str(),bias_pntr_in->getLabel().c_str(),getTime());
//   return hessian_pntr_out;
// }


// CoeffsMatrix* Optimizer::switchToFullHessian(VesBias* bias_pntr_in) {
//   plumed_massert(use_hessian_,"it does not make sense to switch to diagonal Hessian if it Hessian is not used");
//   diagonal_hessian_=false;
//   bias_pntr_in->enableHessian(diagonal_hessian_);
//   CoeffsMatrix* hessian_pntr_out = bias_pntr_in->getHessianPntr();
//   plumed_massert(hessian_pntr_out != NULL,"Hessian is needed but not linked correctly");
//   //
//   log.printf("  %s (with label %s): switching to a diagonal Hessian for VES bias %s (with label %s) at time  %f\n",getName().c_str(),getLabel().c_str(),bias_pntr_in->getName().c_str(),bias_pntr_in->getLabel().c_str(),getTime());
//   return hessian_pntr_out;
// }


void Optimizer::update() {
  if(onStep() && !isFirstStep) {
    for(unsigned int i=0; i<nbiases_; i++) {
      bias_pntrs_[i]->updateGradientAndHessian(use_mwalkers_mpi_);
    }
    for(unsigned int i=0; i<ncoeffssets_; i++) {
      if(gradient_pntrs_[i]->isActive()) {coeffsUpdate(i);}
      else {
        std::string msg = "iteration " + getIterationCounterStr(+1) +
                          " for " + bias_pntrs_[i]->getLabel() +
                          " - the coefficients are not updated as CV values are outside the bias intervals";
        warning(msg);
      }

      // +1 as this is done before increaseIterationCounter() is used
      unsigned int curr_iter = getIterationCounter()+1;
      double curr_time = getTime();
      coeffs_pntrs_[i]->setIterationCounterAndTime(curr_iter,curr_time);
      aux_coeffs_pntrs_[i]->setIterationCounterAndTime(curr_iter,curr_time);
      gradient_pntrs_[i]->setIterationCounterAndTime(curr_iter,curr_time);
      targetdist_averages_pntrs_[i]->setIterationCounterAndTime(curr_iter,curr_time);
      if(use_hessian_) {
        hessian_pntrs_[i]->setIterationCounterAndTime(curr_iter,curr_time);
      }
      if(aver_gradient_pntrs_.size()>0) {
        aver_gradient_pntrs_[i]->setIterationCounterAndTime(curr_iter,curr_time);
        aver_gradient_pntrs_[i]->addToAverage(*gradient_pntrs_[i]);
      }
    }
    increaseIterationCounter();
    updateOutputComponents();
    for(unsigned int i=0; i<ncoeffssets_; i++) {
      writeOutputFiles(i);
    }
    if(ustride_targetdist_>0 && getIterationCounter()%ustride_targetdist_==0) {
      for(unsigned int i=0; i<nbiases_; i++) {
        if(dynamic_targetdists_[i]) {
          bias_pntrs_[i]->updateTargetDistributions();
        }
      }
    }
    if(ustride_reweightfactor_>0 && getIterationCounter()%ustride_reweightfactor_==0) {
      for(unsigned int i=0; i<nbiases_; i++) {
        bias_pntrs_[i]->updateReweightFactor();
      }
    }


    //
    if(isBiasOutputActive() && getIterationCounter()%getBiasOutputStride()==0) {
      writeBiasOutputFiles();
    }
    if(isFesOutputActive() && getIterationCounter()%getFesOutputStride()==0) {
      writeFesOutputFiles();
    }
    if(isFesProjOutputActive() && getIterationCounter()%getFesProjOutputStride()==0) {
      writeFesProjOutputFiles();
    }
    if(isTargetDistOutputActive() && getIterationCounter()%getTargetDistOutputStride()==0) {
      writeTargetDistOutputFiles();
    }
    if(isTargetDistProjOutputActive() && getIterationCounter()%getTargetDistProjOutputStride()==0) {
      writeTargetDistProjOutputFiles();
    }
  }
  else {
    isFirstStep=false;
  }
}


void Optimizer::updateOutputComponents() {
  if(ncoeffssets_==1) {
    if(!fixed_stepsize_) {
      getPntrToComponent("stepsize")->set( getCurrentStepSize(0) );
    }
    if(monitor_instantaneous_gradient_) {
      getPntrToComponent("gradrms")->set( gradient_pntrs_[0]->getRMS() );
      size_t gradient_maxabs_idx=0;
      getPntrToComponent("gradmax")->set( gradient_pntrs_[0]->getMaxAbsValue(gradient_maxabs_idx) );
    }
    if(aver_gradient_pntrs_.size()>0) {
      getPntrToComponent("avergradrms")->set( aver_gradient_pntrs_[0]->getRMS() );
      size_t avergradient_maxabs_idx=0;
      getPntrToComponent("avergradmax")->set( aver_gradient_pntrs_[0]->getMaxAbsValue(avergradient_maxabs_idx) );
    }
  }
  else {
    for(unsigned int i=0; i<ncoeffssets_; i++) {
      std::string is=""; Tools::convert(i,is); is = "-" + coeffssetid_prefix_ + is;
      if(!fixed_stepsize_) {
        getPntrToComponent("stepsize"+is)->set( getCurrentStepSize(i) );
      }
      if(monitor_instantaneous_gradient_) {
        getPntrToComponent("gradrms"+is)->set( gradient_pntrs_[i]->getRMS() );
        size_t gradient_maxabs_idx=0;
        getPntrToComponent("gradmax"+is)->set( gradient_pntrs_[i]->getMaxAbsValue(gradient_maxabs_idx) );
      }
      if(aver_gradient_pntrs_.size()>0) {
        getPntrToComponent("avergradrms"+is)->set( aver_gradient_pntrs_[0]->getRMS() );
        size_t avergradient_maxabs_idx=0;
        getPntrToComponent("avergradmax"+is)->set( aver_gradient_pntrs_[0]->getMaxAbsValue(avergradient_maxabs_idx) );
      }
    }
  }
}


void Optimizer::turnOffCoeffsOutputFiles() {
  for(unsigned int i=0; i<coeffsOFiles_.size(); i++) {
    coeffsOFiles_[i]->close();
    delete coeffsOFiles_[i];
  }
  coeffsOFiles_.clear();
}


void Optimizer::writeOutputFiles(const unsigned int coeffs_id) {
  if(coeffsOFiles_.size()>0 && iter_counter%coeffs_wstride_==0) {
    coeffs_pntrs_[coeffs_id]->writeToFile(*coeffsOFiles_[coeffs_id],aux_coeffs_pntrs_[coeffs_id],false);
  }
  if(gradientOFiles_.size()>0 && iter_counter%gradient_wstride_==0) {
    if(aver_gradient_pntrs_.size()==0) {
      gradient_pntrs_[coeffs_id]->writeToFile(*gradientOFiles_[coeffs_id],false);
    }
    else {
      gradient_pntrs_[coeffs_id]->writeToFile(*gradientOFiles_[coeffs_id],aver_gradient_pntrs_[coeffs_id],false);
    }
  }
  if(hessianOFiles_.size()>0 && iter_counter%hessian_wstride_==0) {
    hessian_pntrs_[coeffs_id]->writeToFile(*hessianOFiles_[coeffs_id]);
  }
  if(targetdist_averagesOFiles_.size()>0 && iter_counter%targetdist_averages_wstride_==0) {
    targetdist_averages_pntrs_[coeffs_id]->writeToFile(*targetdist_averagesOFiles_[coeffs_id]);
  }
}


void Optimizer::setupOFiles(std::vector<std::string>& fnames, std::vector<OFile*>& OFiles, const bool multi_sim_single_files) {
  plumed_assert(ncoeffssets_>0);
  OFiles.resize(fnames.size(),NULL);
  for(unsigned int i=0; i<fnames.size(); i++) {
    OFiles[i] = new OFile();
    OFiles[i]->link(*this);
    if(multi_sim_single_files) {
      unsigned int r=0;
      if(comm.Get_rank()==0) {r=multi_sim_comm.Get_rank();}
      comm.Bcast(r,0);
      if(r>0) {fnames[i]="/dev/null";}
      OFiles[i]->enforceSuffix("");
    }
    OFiles[i]->open(fnames[i]);
    OFiles[i]->setHeavyFlush();
  }
}


void Optimizer::readCoeffsFromFiles(const std::vector<std::string>& fnames, const bool read_aux_coeffs) {
  plumed_assert(ncoeffssets_>0);
  plumed_assert(fnames.size()==ncoeffssets_);
  if(ncoeffssets_==1) {
    log.printf("  Read in coefficients from file ");
  }
  else {
    log.printf("  Read in coefficients from files:\n");
  }
  for(unsigned int i=0; i<ncoeffssets_; i++) {
    IFile ifile;
    ifile.link(*this);
    if(use_mwalkers_mpi_ && mwalkers_mpi_single_files_) {
      ifile.enforceSuffix("");
    }
    ifile.open(fnames[i]);
    if(!ifile.FieldExist(coeffs_pntrs_[i]->getDataLabel())) {
      std::string error_msg = "Problem with reading coefficients from file " + ifile.getPath() + ": no field with name " + coeffs_pntrs_[i]->getDataLabel() + "\n";
      plumed_merror(error_msg);
    }
    size_t ncoeffs_read = coeffs_pntrs_[i]->readFromFile(ifile,false,false);
    if(ncoeffssets_==1) {
      log.printf("%s (read %zu of %zu values)\n", ifile.getPath().c_str(),ncoeffs_read,coeffs_pntrs_[i]->numberOfCoeffs());
    }
    else {
      log.printf("   coefficient set %u: %s (read %zu of %zu values)\n",i,ifile.getPath().c_str(),ncoeffs_read,coeffs_pntrs_[i]->numberOfCoeffs());
    }
    ifile.close();
    if(read_aux_coeffs) {
      ifile.open(fnames[i]);
      if(!ifile.FieldExist(aux_coeffs_pntrs_[i]->getDataLabel())) {
        std::string error_msg = "Problem with reading coefficients from file " + ifile.getPath() + ": no field with name " + aux_coeffs_pntrs_[i]->getDataLabel() + "\n";
        plumed_merror(error_msg);
      }
      aux_coeffs_pntrs_[i]->readFromFile(ifile,false,false);
      ifile.close();
    }
    else {
      AuxCoeffs(i).setValues( Coeffs(i) );
    }
  }
}


void Optimizer::addCoeffsSetIDsToFilenames(std::vector<std::string>& fnames, std::string& coeffssetid_prefix) {
  if(ncoeffssets_==1) {return;}
  //
  if(fnames.size()==1) {
    fnames.resize(ncoeffssets_,fnames[0]);
  }
  plumed_assert(fnames.size()==ncoeffssets_);
  //
  for(unsigned int i=0; i<ncoeffssets_; i++) {
    std::string is=""; Tools::convert(i,is);
    fnames[i] = FileBase::appendSuffix(fnames[i],"."+coeffssetid_prefix_+is);
  }
}


void Optimizer::setAllCoeffsSetIterationCounters() {
  for(unsigned int i=0; i<ncoeffssets_; i++) {
    coeffs_pntrs_[i]->setIterationCounterAndTime(getIterationCounter(),getTime());
    aux_coeffs_pntrs_[i]->setIterationCounterAndTime(getIterationCounter(),getTime());
    gradient_pntrs_[i]->setIterationCounterAndTime(getIterationCounter(),getTime());
    targetdist_averages_pntrs_[i]->setIterationCounterAndTime(getIterationCounter(),getTime());
    if(use_hessian_) {
      hessian_pntrs_[i]->setIterationCounterAndTime(getIterationCounter(),getTime());
    }
  }
}


std::string Optimizer::getIterationCounterStr(const int offset) const {
  std::string s;
  Tools::convert(getIterationCounter()+offset,s);
  return s;
}


void Optimizer::writeBiasOutputFiles() const {
  for(unsigned int i=0; i<nbiases_; i++) {
    bias_pntrs_[i]->writeBiasToFile();
  }
}


void Optimizer::writeFesOutputFiles() const {
  for(unsigned int i=0; i<nbiases_; i++) {
    bias_pntrs_[i]->writeFesToFile();
  }
}


void Optimizer::writeFesProjOutputFiles() const {
  for(unsigned int i=0; i<nbiases_; i++) {
    bias_pntrs_[i]->writeFesProjToFile();
  }
}


void Optimizer::writeTargetDistOutputFiles() const {
  for(unsigned int i=0; i<nbiases_; i++) {
    if(dynamic_targetdists_[i]) {
      bias_pntrs_[i]->writeTargetDistToFile();
    }
  }
}


void Optimizer::writeTargetDistProjOutputFiles() const {
  for(unsigned int i=0; i<nbiases_; i++) {
    if(dynamic_targetdists_[i]) {
      bias_pntrs_[i]->writeTargetDistProjToFile();
    }
  }
}



}
}
