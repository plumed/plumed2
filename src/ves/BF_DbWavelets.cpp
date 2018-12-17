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


#include "BasisFunctions.h"
#include "DbWaveletGrid.h"
#include "tools/Grid.h"
#include "VesTools.h"
#include "core/ActionRegister.h"
#include "tools/Exception.h"

#include <algorithm>


namespace PLMD {
namespace ves {


//+PLUMEDOC VES_BASISF BF_DB_WAVELETS
/*
Daubechies Wavelets as basis functions

Note: at the moment only the scaling function is working as intended as multiscale is not yet implemented.

This basis set uses Daubechies Wavelets to construct a complete and orthogonal basis.
Because no analytic formula for these wavelets exist, they are instead constructed iteratively on a grid.
The method of construction is close to the "Vector cascade algorithm" described by Strang, Nguyen. (It is sometimes also called Daubechies-Lagarias method)

\par Input parameters

order N: number of vanishing moments

Intrinsic support of the function is then [0, 2*N-1).

(Give formula!)

Each basis function is a translate by an integer value k.

\par Number of basis functions

If the support is scaled to match the desired range of the CV exactly there would be 4*N -3 basis functions whose support is at least partially in this region: k = {(-2*N)+2, … , 0, … 2*N-1}
As some of these translates will not have significant contributions in this area because of their function values being close to zero over the full range, these 'tail wavelets' are being omitted.
The formula for cutting of was found empirically. In short it omits roughly all translates that have only function values with less than 1 % of the maximum value in the CV range.
It is not perfect but sometimes keeps a few more translates (depending on the actually chosen order).
Giving the range of k in dependency of the chosen order is therefore not possible.
Instead it will be given in the logfile if wanted.
As a rule of thumb there are about 3*order basis functions.

(Note for future: The number of needed BFs could also be lowered by scaling the wavelets less so that their support is larger than the desired CV range.)

Additionally a constant basis function is included.


\par Examples

\par Test

*/
//+ENDPLUMEDOC


class BF_DbWavelets : public BasisFunctions {
  // Grid that holds the Wavelet values and its derivative
  std::unique_ptr<Grid> waveletGrid_;
  bool use_scaling_function_;
  double scale_; // scale factor of the individual BFs to match specified length
  std::vector<double> shifts_; // shift of the individual BFs
  void setupLabels() override;
protected:
  std::vector<double> getCutoffPoints(const double& threshold);

public:
  static void registerKeywords( Keywords&);
  explicit BF_DbWavelets(const ActionOptions&);
  void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const override;
};


PLUMED_REGISTER_ACTION(BF_DbWavelets,"BF_DB_WAVELETS")


void BF_DbWavelets::registerKeywords(Keywords& keys) {
  BasisFunctions::registerKeywords(keys);
  keys.add("optional","GRID_SIZE","The number of grid bins of the Wavelet function. Because of the used construction algorithm this value will be used as guiding value only, while the true number will be \"(ORDER*2 - 1) * 2**n\" with the smallest n such that the grid is at least as large as the specified number. Defaults to 1000"); // Change the commentary a bit?
  keys.add("optional","FUNCTION_LENGTH","The length of the support of the scaled basis functions. This can be used to alter the scaling of the basis functions. Is by default set to the total size of the interval. This also influences the number of actually used basis functions, as all shifted functions that are partially supported in the CV space are used.");
  keys.add("optional","TAILS_THRESHOLD","The threshold for cutting off tail wavelets with respect to the maximum value. All shifted wavelet functions that will only have values lower below the threshold in the CV space will be excluded from the basis set. Defaults to 0 (include all), generally only small values (e.g. 0.01) should be used.");
  keys.addFlag("SCALING_FUNCTION", false, "If this flag is set the scaling function (mother wavelet) will be used instead of the \"true\" wavelet function (father wavelet).");
  keys.addFlag("DUMP_WAVELET_GRID", false, "If this flag is set the grid with the wavelet values will be written to a file called \"wavelet_grid.data\".");
  // why is this removed?
  keys.remove("NUMERICAL_INTEGRALS");
}


BF_DbWavelets::BF_DbWavelets(const ActionOptions& ao):
  PLUMED_VES_BASISFUNCTIONS_INIT(ao),
  use_scaling_function_(false)
{

  // parse grid properties and set it up
  parseFlag("SCALING_FUNCTION", use_scaling_function_);
  unsigned gridsize = 1000;
  parse("GRID_SIZE", gridsize);
  waveletGrid_ = DbWaveletGrid::setupGrid(getOrder(), gridsize, !use_scaling_function_);
  unsigned true_gridsize = waveletGrid_->getNbin()[0];
  if(true_gridsize != 1000) {addKeywordToList("GRID_SIZE",true_gridsize);}
  bool dump_wavelet_grid=false;
  parseFlag("DUMP_WAVELET_GRID", dump_wavelet_grid);
  if (dump_wavelet_grid) {
    OFile wavelet_gridfile;
    wavelet_gridfile.link(*this);
    wavelet_gridfile.enforceBackup();
    wavelet_gridfile.open(getLabel()+".wavelet_grid.data");
    waveletGrid_->writeToFile(wavelet_gridfile);
  }

  // calculate the number of basis functions from the specified length
  unsigned intrinsic_length = 2*getOrder() - 1;
  double length = intervalMax() - intervalMin();
  parse("FUNCTION_LENGTH",length);
  if(length != intervalMax() - intervalMin()) {addKeywordToList("FUNCTION_LENGTH",length);}
  scale_ = intrinsic_length / length;

  // parse threshold for tail wavelets and get respective cutoff points
  double threshold = 0.0;
  std::vector<double> cutoffpoints (2);
  parse("TAILS_THRESHOLD",threshold);
  if(threshold >= 1) {plumed_merror("TAILS_THRESHOLD should be significantly smaller than 1.");}
  if(threshold == 0.0) {
    cutoffpoints = {0.0, static_cast<double>(intrinsic_length)};
  }
  else {
    addKeywordToList("TAILS_THRESHOLD",threshold);
    cutoffpoints = getCutoffPoints(threshold);
  };

  // calculate number of Basis functions and the shifts
  unsigned int num_BFs = 1; // constant one
  num_BFs += static_cast<unsigned>(ceil(cutoffpoints[1])); // left shifts including 0
  num_BFs += static_cast<unsigned>(ceil((intervalMax()-intervalMin())*scale_ - cutoffpoints[0] - 1)); // right shifts
  setNumberOfBasisFunctions(num_BFs);
  shifts_.push_back(0.0); // constant BF – never used, just for clearer notation
  for(unsigned int i = 1; i < getNumberOfBasisFunctions(); ++i) {
    shifts_.push_back(-intervalMin()*scale_ + ceil(cutoffpoints[1]) - i);
  }

  // set some properties
  setIntrinsicInterval(0.0,intrinsic_length);
  setNonPeriodic();
  setIntervalBounded();
  setType("daubechies_wavelets");
  setDescription("Daubechies Wavelets (maximum phase type)");
  setLabelPrefix("k");
  setupBF();
  checkRead();
}


void BF_DbWavelets::getAllValues(const double arg, double& argT, bool& inside_range, std::vector<double>& values, std::vector<double>& derivs) const {
  argT=checkIfArgumentInsideInterval(arg,inside_range);
  //
  values[0]=1.0;
  derivs[0]=0.0;
  for(unsigned int i = 1; i < getNumberOfBasisFunctions(); ++i) {
    // scale and shift argument to match current wavelet
    double x = shifts_[i] + arg*scale_;

    if (x < 0 || x >= intrinsicIntervalMax()) { // Wavelets are 0 outside the defined range
      values[i] = 0.0; derivs[i] = 0.0;
    }
    else {
      // declare vectors and fill them with value
      std::vector<double> temp_deriv (1);
      std::vector<double> x_vec {x};
      values[i] = waveletGrid_->getValueAndDerivatives(x_vec, temp_deriv);
      derivs[i] = temp_deriv[0] * intervalDerivf(); // scale derivative
    }
  }
  if(!inside_range) {for(auto& deriv : derivs) {deriv=0.0;}}
}


// returns left and right cutoff point of Wavelet
// threshold is a percent value of maximum
std::vector<double> BF_DbWavelets::getCutoffPoints(const double& threshold) {
  double threshold_value = threshold * waveletGrid_->getMaxValue();
  std::vector<double> cutoffpoints;

  for (unsigned i = 0; i < waveletGrid_->getSize(); ++i) {
    if (fabs(waveletGrid_->getValue(i)) >= threshold_value) {
      cutoffpoints.push_back(waveletGrid_->getPoint(i)[0]);
      break;
    }
  }

  for (int i = waveletGrid_->getSize() - 1; i >= 0; --i) {
    if (fabs(waveletGrid_->getValue(i)) >= threshold_value) {
      cutoffpoints.push_back(waveletGrid_->getPoint(i)[0]);
      break;
    }
  }

  return cutoffpoints;
}


// labels according to minimum position in CV space
void BF_DbWavelets::setupLabels() {
  setLabel(0,"const");
  for(unsigned int i=1; i < getNumberOfBasisFunctions(); i++) {
    double pos = -shifts_[i]/scale_;
    std::string is; Tools::convert(pos, is);
    setLabel(i,"i="+is);
  }
}


}
}
