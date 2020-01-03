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
#include "LinearBasisSetExpansion.h"
#include "CoeffsVector.h"
#include "CoeffsMatrix.h"
#include "BasisFunctions.h"
#include "Optimizer.h"
#include "TargetDistribution.h"
#include "VesTools.h"

#include "bias/Bias.h"
#include "core/ActionRegister.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"


namespace PLMD {
namespace ves {

//+PLUMEDOC VES_BIAS VES_LINEAR_EXPANSION
/*
Linear basis set expansion bias.

This VES bias action takes the bias potential to be a linear expansion
in some basis set that is written as a product of one-dimensional basis functions.
For example, for one CV the bias would be written as
\f[
V(s_{1};\boldsymbol{\alpha}) = \sum_{i_{1}} \alpha_{i_{1}} \, f_{i_{1}}(s_{1}),
\f]
while for two CVs it is written as
\f[
V(s_{1},s_{2};\boldsymbol{\alpha}) = \sum_{i_{1},i_{2}} \alpha_{i_{1},i_{2}} \, f_{i_{1}}(s_{1}) \, f_{i_{2}}(s_{2})
\f]
where \f$\boldsymbol{\alpha}\f$ is the set of expansion coefficients that
are optimized within VES. With an appropriate choice of the basis functions
it is possible to represent any generic free energy surface.
The relationship between the bias and the free energy surface is given by
\f[
V(\mathbf{s}) = - F(\mathbf{s}) - \frac{1}{\beta} \log p(\mathbf{s}).
\f]
where \f$p(\mathbf{s})\f$ is the target distribution that is employed in the VES simulation.

\par Basis Functions

Various one-dimensional basis functions are available in the VES code,
see the complete list \ref ves_basisf "here".
At the current moment we recommend to use Legendre polynomials (\ref BF_LEGENDRE)
for non-periodic CVs and Fourier basis functions (\ref BF_FOURIER)
for periodic CV (e.g. dihedral angles).

To use basis functions within VES_LINEAR_EXPANSION you first need to
define them in the input file before the VES_LINEAR_EXPANSION action and
then give their labels using the BASIS_FUNCTIONS keyword.


\par Target Distributions

Various target distributions \f$p(\mathbf{s})\f$ are available in the VES code,
see the complete list \ref ves_targetdist "here".

To use a target distribution within VES_LINEAR_EXPANSION you first need to
define it in the input file before the VES_LINEAR_EXPANSION action and
then give its label using the TARGET_DISTRIBUTION keyword.
The default behavior if no TARGET_DISTRIBUTION is given is to
employ a uniform target distribution.

Some target distribution, like the well-tempered one (\ref TD_WELLTEMPERED),
are dynamic and need to be iteratively updated during the optimization.

\par Optimizer

In order to optimize the coefficients you will need to use VES_LINEAR_EXPANSION
in combination with an optimizer, see the list of optimizers available in the
VES code \ref ves_optimizer "here". At the current moment we recommend to
use the averaged stochastic gradient decent optimizer (\ref OPT_AVERAGED_SGD).

The optimizer should be defined after the VES_LINEAR_EXPANSION action.

\par Grid

Internally the code uses grids to calculate the basis set averages
over the target distribution that is needed for the gradient. The same grid is
also used for the output files (see next section).
The size of the grid is determined by the GRID_BINS keyword. By default it has
100 grid points in each dimension, and generally this value should be sufficient.

\par Outputting Free Energy Surfaces and Other Files

It is possible to output on-the-fly during the simulation the free energy surface
estimated from the bias potential. How often this is done is specified within
the \ref ves_optimizer "optimizer" by using the FES_OUTPUT keyword. The filename
is specified by the FES_FILE keyword, but by default is it fes.LABEL.data,
with an added suffix indicating
the iteration number (iter-#).

For multi-dimensional case is it possible to also output projections of the
free energy surfaces. The arguments for which to do these projections is
specified using the numbered PROJ_ARG keywords. For these files a suffix
indicating the projection (proj-#) will be added to the filenames.
You will also need to specify the frequency of the output by using the
FES_PROJ_OUTPUT keyword within the optimizer.

It is also possible to output the bias potential itself, for this the relevant
keyword is BIAS_OUTPUT within the optimizer. The filename
is specified by the BIAS_FILE keyword, but by default is it bias.LABEL.data,
with an added suffix indicating the iteration number (iter-#).

Furthermore is it possible to output the target distribution, and its projections
(i.e. marginal distributions). The filenames of these files are specified with
the TARGETDIST_FILE, but by default is it targetdist.LABEL.data. The
logarithm of the target distribution will also be outputted to file that has the
added suffix log. For static target distribution these files will be outputted in
the beginning of the
simulation while for dynamic ones you will need to specify the frequency
of the output by using the TARGETDIST_OUTPUT and TARGETDIST_PROJ_OUTPUT
keywords within the optimizer.

It is also possible to output free energy surfaces and bias in post processing
by using the \ref VES_OUTPUT_FES action. However, be aware that this action
does does not support dynamic target distribution (e.g. well-tempered).

\par Static Bias

It is also possible to use VES_LINEAR_EXPANSION as a static bias that uses
previously obtained coefficients. In this case the coefficients should be
read in from the coefficient file given in the COEFFS keyword.

\par Bias Cutoff

It is possible to impose a cutoff on the bias potential using the procedure
introduced in \cite McCarty-PRL-2015 such that the free energy surface
is only flooded up to a certain value. The bias that results from this procedure
can then be used as a static bias for obtaining kinetic rates.
The value of the cutoff is given by the BIAS_CUTOFF keyword.
To impose the cutoff the code uses a Fermi switching function \f$1/(1+e^{\lambda x})\f$
where the parameter \f$\lambda\f$ controls how sharply the switchingfunction goes to zero.
The default value is \f$\lambda=10\f$ but this can be changed by using the
BIAS_CUTOFF_FERMI_LAMBDA keyword.

\par Examples

In the following example we run a VES_LINEAR_EXPANSION for one CV using
a Legendre basis functions (\ref BF_LEGENDRE) and a uniform target
distribution as no target distribution is specified. The coefficients
are optimized using averaged stochastic gradient descent optimizer
(\ref OPT_AVERAGED_SGD). Within the optimizer we specify that the
FES should be outputted to file every 500 coefficients iterations (the
FES_OUTPUT keyword).
Parameters that are very specific to the problem at hand, like the
order of the basis functions, the interval on which the
basis functions are defined, and the step size used
in the optimizer, are left unfilled.
\plumedfile
bf1: BF_LEGENDRE ORDER=__ MINIMUM=__ MAXIMUM=__

VES_LINEAR_EXPANSION ...
 ARG=d1
 BASIS_FUNCTIONS=bf1
 TEMP=__
 GRID_BINS=200
 LABEL=b1
... VES_LINEAR_EXPANSION

OPT_AVERAGED_SGD ...
 BIAS=b1
 STRIDE=1000
 LABEL=o1
 STEPSIZE=__
 FES_OUTPUT=500
 COEFFS_OUTPUT=10
... OPT_AVERAGED_SGD
\endplumedfile

In the following example we employ VES_LINEAR_EXPANSION for two CVs,
The first CV is periodic and therefore we employ a Fourier basis functions
(\ref BF_LEGENDRE) while the second CV is non-periodic so we employ a
Legendre polynomials as in the previous example. For the target distribution
we employ a well-tempered target distribution (\ref TD_WELLTEMPERED), which is
dynamic and needs to be iteratively updated with a stride that is given
using the TARGETDIST_STRIDE within the optimizer.

\plumedfile
bf1: BF_FOURIER  ORDER=__ MINIMUM=__ MAXIMUM=__
bf2: BF_LEGENDRE ORDER=__ MINIMUM=__ MAXIMUM=__

td_wt: TD_WELLTEMPERED BIASFACTOR=10.0

VES_LINEAR_EXPANSION ...
 ARG=cv1,cv2
 BASIS_FUNCTIONS=bf1,bf2
 TEMP=__
 GRID_BINS=100
 LABEL=b1
 TARGET_DISTRIBUTION=td_wt
... VES_LINEAR_EXPANSION

OPT_AVERAGED_SGD ...
 BIAS=b1
 STRIDE=1000
 LABEL=o1
 STEPSIZE=__
 FES_OUTPUT=500
 COEFFS_OUTPUT=10
 TARGETDIST_STRIDE=500
... OPT_AVERAGED_SGD
\endplumedfile


In the following example we employ a bias cutoff such that the bias
only fills the free energy landscape up a certain level. In this case
the target distribution is also dynamic and needs to iteratively updated.

\plumedfile
bf1: BF_LEGENDRE ORDER=__ MINIMUM=__ MAXIMUM=__
bf2: BF_LEGENDRE ORDER=__ MINIMUM=__ MAXIMUM=__

VES_LINEAR_EXPANSION ...
 ARG=cv1,cv2
 BASIS_FUNCTIONS=bf1,bf2
 TEMP=__
 GRID_BINS=100
 LABEL=b1
 BIAS_CUTOFF=20.0
... VES_LINEAR_EXPANSION

OPT_AVERAGED_SGD ...
 BIAS=b1
 STRIDE=1000
 LABEL=o1
 STEPSIZE=__
 FES_OUTPUT=500
 COEFFS_OUTPUT=10
 TARGETDIST_STRIDE=500
... OPT_AVERAGED_SGD
\endplumedfile

The optimized bias potential can then be used as a static bias for obtaining
kinetics. For this you need read in the final coefficients from file
(e.g. coeffs_final.data in this case) by using the
COEFFS keyword (also, no optimizer should be defined in the input)
\plumedfile
bf1: BF_LEGENDRE ORDER=__ MINIMUM=__ MAXIMUM=__
bf2: BF_LEGENDRE ORDER=__ MINIMUM=__ MAXIMUM=__

VES_LINEAR_EXPANSION ...
 ARG=cv1,cv2
 BASIS_FUNCTIONS=bf1,bf2
 TEMP=__
 GRID_BINS=100
 LABEL=b1
 BIAS_CUTOFF=20.0
 COEFFS=coeffs_final.data
... VES_LINEAR_EXPANSION
\endplumedfile



*/
//+ENDPLUMEDOC


class VesLinearExpansion : public VesBias {
private:
  unsigned int nargs_;
  std::vector<BasisFunctions*> basisf_pntrs_;
  LinearBasisSetExpansion* bias_expansion_pntr_;
  size_t ncoeffs_;
  Value* valueForce2_;
  bool all_values_inside;
  std::vector<double> bf_values;
  bool bf_values_set;
public:
  explicit VesLinearExpansion(const ActionOptions&);
  ~VesLinearExpansion();
  void calculate() override;
  void update() override;
  void updateTargetDistributions() override;
  void restartTargetDistributions() override;
  //
  void setupBiasFileOutput() override;
  void writeBiasToFile() override;
  void resetBiasFileOutput() override;
  //
  void setupFesFileOutput() override;
  void writeFesToFile() override;
  void resetFesFileOutput() override;
  //
  void setupFesProjFileOutput() override;
  void writeFesProjToFile() override;
  //
  void writeTargetDistToFile() override;
  void writeTargetDistProjToFile() override;
  //
  double calculateReweightFactor() const override;
  //
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(VesLinearExpansion,"VES_LINEAR_EXPANSION")

void VesLinearExpansion::registerKeywords( Keywords& keys ) {
  VesBias::registerKeywords(keys);
  //
  VesBias::useInitialCoeffsKeywords(keys);
  VesBias::useTargetDistributionKeywords(keys);
  VesBias::useBiasCutoffKeywords(keys);
  VesBias::useGridBinKeywords(keys);
  VesBias::useProjectionArgKeywords(keys);
  //
  keys.use("ARG");
  keys.add("compulsory","BASIS_FUNCTIONS","the label of the one dimensional basis functions that should be used.");
  keys.addOutputComponent("force2","default","the instantaneous value of the squared force due to this bias potential.");
}

VesLinearExpansion::VesLinearExpansion(const ActionOptions&ao):
  PLUMED_VES_VESBIAS_INIT(ao),
  nargs_(getNumberOfArguments()),
  basisf_pntrs_(0),
  bias_expansion_pntr_(NULL),
  valueForce2_(NULL),
  all_values_inside(true),
  bf_values(0),
  bf_values_set(false)
{
  std::vector<std::string> basisf_labels;
  parseMultipleValues("BASIS_FUNCTIONS",basisf_labels,nargs_);
  checkRead();

  std::string error_msg = "";
  basisf_pntrs_ = VesTools::getPointersFromLabels<BasisFunctions*>(basisf_labels,plumed.getActionSet(),error_msg);
  if(error_msg.size()>0) {plumed_merror("Error in keyword BASIS_FUNCTIONS of "+getName()+": "+error_msg);}
  //

  std::vector<Value*> args_pntrs = getArguments();
  // check arguments and basis functions
  // this is done to avoid some issues with integration of target distribution
  // and periodic CVs, needs to be fixed later on.
  for(unsigned int i=0; i<args_pntrs.size(); i++) {
    if(args_pntrs[i]->isPeriodic() && !(basisf_pntrs_[i]->arePeriodic()) ) {
      plumed_merror("argument "+args_pntrs[i]->getName()+" is periodic while the basis functions " + basisf_pntrs_[i]->getLabel()+ " are not. You need to use the COMBINE action to remove the periodicity of the argument if you want to use these basis functions");
    }
    else if(!(args_pntrs[i]->isPeriodic()) && basisf_pntrs_[i]->arePeriodic() ) {
      log.printf("  warning: argument %s is not periodic while the basis functions %s used for it are periodic\n",args_pntrs[i]->getName().c_str(),basisf_pntrs_[i]->getLabel().c_str());
    }
  }

  addCoeffsSet(args_pntrs,basisf_pntrs_);
  ncoeffs_ = numberOfCoeffs();
  bool coeffs_read = readCoeffsFromFiles();

  checkThatTemperatureIsGiven();
  bias_expansion_pntr_ = new LinearBasisSetExpansion(getLabel(),getBeta(),comm,args_pntrs,basisf_pntrs_,getCoeffsPntr());
  bias_expansion_pntr_->linkVesBias(this);
  bias_expansion_pntr_->setGridBins(this->getGridBins());
  //
  bf_values.assign(ncoeffs_,0.0);



  if(getNumberOfTargetDistributionPntrs()==0) {
    log.printf("  using an uniform target distribution: \n");
    bias_expansion_pntr_->setupUniformTargetDistribution();
    disableStaticTargetDistFileOutput();
  }
  else if(getNumberOfTargetDistributionPntrs()==1) {
    if(biasCutoffActive()) {getTargetDistributionPntrs()[0]->setupBiasCutoff();}
    bias_expansion_pntr_->setupTargetDistribution(getTargetDistributionPntrs()[0]);
    log.printf("  using target distribution of type %s with label %s \n",getTargetDistributionPntrs()[0]->getName().c_str(),getTargetDistributionPntrs()[0]->getLabel().c_str());
  }
  else {
    plumed_merror("problem with the TARGET_DISTRIBUTION keyword, either give no label or just one label.");
  }
  setTargetDistAverages(bias_expansion_pntr_->TargetDistAverages());
  //
  if(coeffs_read && biasCutoffActive()) {
    VesLinearExpansion::updateTargetDistributions();
  }
  //
  if(coeffs_read) {
    VesLinearExpansion::setupBiasFileOutput();
    VesLinearExpansion::writeBiasToFile();
  }

  addComponent("force2"); componentIsNotPeriodic("force2");
  valueForce2_=getPntrToComponent("force2");
}


VesLinearExpansion::~VesLinearExpansion() {
  if(bias_expansion_pntr_!=NULL) {
    delete bias_expansion_pntr_;
  }
}


void VesLinearExpansion::calculate() {

  std::vector<double> cv_values(nargs_);
  std::vector<double> forces(nargs_);

  for(unsigned int k=0; k<nargs_; k++) {
    cv_values[k]=getArgument(k);
  }

  all_values_inside = true;
  double bias = bias_expansion_pntr_->getBiasAndForces(cv_values,all_values_inside,forces,bf_values);
  if(biasCutoffActive()) {
    applyBiasCutoff(bias,forces,bf_values);
    bf_values[0]=1.0;
  }
  double totalForce2 = 0.0;
  for(unsigned int k=0; k<nargs_; k++) {
    setOutputForce(k,forces[k]);
    totalForce2 += forces[k]*forces[k];
  }

  setBias(bias);
  valueForce2_->set(totalForce2);

  bf_values_set = true;
}


void VesLinearExpansion::update() {
  if(!bf_values_set) {
    warning("VesLinearExpansion::update() is being called without calling VesLinearExpansion::calculate() first to calculate the basis function values. This can lead to incorrect behavior.");
  }
  if(all_values_inside && bf_values_set) {
    addToSampledAverages(bf_values);
  }
  std::fill(bf_values.begin(), bf_values.end(), 0.0);
  bf_values_set = false;
}






void VesLinearExpansion::updateTargetDistributions() {
  bias_expansion_pntr_->updateTargetDistribution();
  setTargetDistAverages(bias_expansion_pntr_->TargetDistAverages());
}


void VesLinearExpansion::restartTargetDistributions() {
  bias_expansion_pntr_->readInRestartTargetDistribution(getCurrentTargetDistOutputFilename());
  bias_expansion_pntr_->restartTargetDistribution();
  setTargetDistAverages(bias_expansion_pntr_->TargetDistAverages());
}


void VesLinearExpansion::setupBiasFileOutput() {
  bias_expansion_pntr_->setupBiasGrid(true);
}


void VesLinearExpansion::writeBiasToFile() {
  bias_expansion_pntr_->updateBiasGrid();
  OFile* ofile_pntr = getOFile(getCurrentBiasOutputFilename(),useMultipleWalkers());
  bias_expansion_pntr_->writeBiasGridToFile(*ofile_pntr);
  ofile_pntr->close(); delete ofile_pntr;
  if(biasCutoffActive()) {
    bias_expansion_pntr_->updateBiasWithoutCutoffGrid();
    OFile* ofile_pntr2 = getOFile(getCurrentBiasOutputFilename("without-cutoff"),useMultipleWalkers());
    bias_expansion_pntr_->writeBiasWithoutCutoffGridToFile(*ofile_pntr2);
    ofile_pntr2->close(); delete ofile_pntr2;
  }
}


void VesLinearExpansion::resetBiasFileOutput() {
  bias_expansion_pntr_->resetStepOfLastBiasGridUpdate();
}


void VesLinearExpansion::setupFesFileOutput() {
  bias_expansion_pntr_->setupFesGrid();
}


void VesLinearExpansion::writeFesToFile() {
  bias_expansion_pntr_->updateFesGrid();
  OFile* ofile_pntr = getOFile(getCurrentFesOutputFilename(),useMultipleWalkers());
  bias_expansion_pntr_->writeFesGridToFile(*ofile_pntr);
  ofile_pntr->close(); delete ofile_pntr;
}


void VesLinearExpansion::resetFesFileOutput() {
  bias_expansion_pntr_->resetStepOfLastFesGridUpdate();
}


void VesLinearExpansion::setupFesProjFileOutput() {
  if(getNumberOfProjectionArguments()>0) {
    bias_expansion_pntr_->setupFesProjGrid();
  }
}


void VesLinearExpansion::writeFesProjToFile() {
  bias_expansion_pntr_->updateFesGrid();
  for(unsigned int i=0; i<getNumberOfProjectionArguments(); i++) {
    std::string suffix;
    Tools::convert(i+1,suffix);
    suffix = "proj-" + suffix;
    OFile* ofile_pntr = getOFile(getCurrentFesOutputFilename(suffix),useMultipleWalkers());
    std::vector<std::string> args = getProjectionArgument(i);
    bias_expansion_pntr_->writeFesProjGridToFile(args,*ofile_pntr);
    ofile_pntr->close(); delete ofile_pntr;
  }
}


void VesLinearExpansion::writeTargetDistToFile() {
  OFile* ofile1_pntr = getOFile(getCurrentTargetDistOutputFilename(),useMultipleWalkers());
  OFile* ofile2_pntr = getOFile(getCurrentTargetDistOutputFilename("log"),useMultipleWalkers());
  bias_expansion_pntr_->writeTargetDistGridToFile(*ofile1_pntr);
  bias_expansion_pntr_->writeLogTargetDistGridToFile(*ofile2_pntr);
  ofile1_pntr->close(); delete ofile1_pntr;
  ofile2_pntr->close(); delete ofile2_pntr;
}


void VesLinearExpansion::writeTargetDistProjToFile() {
  for(unsigned int i=0; i<getNumberOfProjectionArguments(); i++) {
    std::string suffix;
    Tools::convert(i+1,suffix);
    suffix = "proj-" + suffix;
    OFile* ofile_pntr = getOFile(getCurrentTargetDistOutputFilename(suffix),useMultipleWalkers());
    std::vector<std::string> args = getProjectionArgument(i);
    bias_expansion_pntr_->writeTargetDistProjGridToFile(args,*ofile_pntr);
    ofile_pntr->close(); delete ofile_pntr;
  }
}


double VesLinearExpansion::calculateReweightFactor() const {
  return bias_expansion_pntr_->calculateReweightFactor();
}


}
}
