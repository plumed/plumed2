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


namespace PLMD {
namespace ves {


//+PLUMEDOC VES_BASISF BF_DB_WAVELETS
/*
Daubechies Wavelets as basis functions.

Note: at the moment only the scaling function is working as intended as multiscale is not yet implemented.

This basis set uses Daubechies Wavelets \cite daubechies_orthonormal_1988 to construct a complete and orthogonal basis.
It is based on using a pair of functions, the scaling function (or father wavelet) \f$\phi\f$ and the wavelet function (or mother wavelet) \f$\psi\f$.
They are defined via the two-scale relations for scale \f$j\f$ and shift \f$k\f$:

\f{align*}{
  \phi_k^j \left(x\right) = 2^{-j/2} \phi \left( 2^{-j} x - k\right)\\
  \psi_k^j \left(x\right) = 2^{-j/2} \psi \left( 2^{-j} x - k\right)
\f}

The exact properties are set by choosing filter coefficients, e.g. choosing \f$h_k\f$ for the father wavelet:

\f[
  \phi\left(x\right) = \sqrt{2} \sum_k h_k\, \phi \left( 2 x - k\right)
\f]

The filter coefficients by Daubechies result in an orthonormal basis of all integer shifted functions:
\f[
  \int \phi(x+i) \phi(x+j) \mathop{}\!\mathrm{d}x = \delta_{ij} \quad \text{for} \quad i,j \in \mathbb{Z}
\f]

Because no analytic formula for these wavelets exist, they are instead constructed iteratively on a grid.
The method of construction is close to the "Vector cascade algorithm" described in \cite strang_wavelets_1997 .
The needed filter coefficients of the scaling function are hardcoded, and were previously generated via a python script.
Currently only the "maximum phase" type is implemented, but the "least asymmetric" type can be added easily.

As an example the two Db8 wavelet functions can be seen below

\image html ves-basisf-db8.png


\par Some details on the input parameters

The specified order of the basis set defines the coefficients and the corresponding wavelet used.
Order N results in DbN wavelets, which is equal to the number of vanishing moments of the wavelet basis and double the number of filter coefficients.

The intrinsic support of the wavelets is then \f$ \left[0, N*2 -1 \right) \f$.
Using the cascade algorithm results in a doubling of the grid values per integer for each iteration.
This means that the grid size will always be a power of two multiplied by the intrinsic support length of the wavelets.
The used grid size is calculated by \f$n_{\text{bins}} = (N*2 - 1) * 2^m\f$ with the smallest \f$ m \f$ such that the grid is at least as large as the specified number.

By default the basis functions are scaled to match the specified size of the CV space (MINIMUM and MAXIMUM keywords), which often is a good initial choice.
The "FUNCTION_LENGTH" keyword can be used to alter this and use a different scaling.
A smaller value will use more and smaller basis functions which results in a more localized optimization, while a larger one will use less and more globally defined functions.

\par Number of basis functions

The resulting basis set consists of integer shifts of the wavelet function at scale \f$j\f$ ,
\f[
  \phi_i (x) = \phi(\frac{x+i}{j})
\f]

where \f$i\f$ in principal would span all positive and negative integers.
Because of the compact support of the wavelets clearly not all shifts are needed.

If the wavelets are scaled to match the CV range exactly there would be \f$4*N -3\f$ basis functions whose support is at least partially in this region.
This number is adjusted automatically if a different FUNCTION_LENGTH is specified.
Additionally, some of the shifted basis functions will not have significant contributions because of their function values being close to zero over the full range.
These 'tail wavelets' can be omitted by using the TAILS_THRESHOLD keyword.
By default all are included but a value of e.g. 0.01 will already reduce the number of basis functions significantly.
The number of basis functions is not easily determinable a priori but will be given in the logfile.
Additionally the starting point (leftmost defined point) of the individual basis functions is printed.

Additionally a constant basis function is included.

\par Examples


First a very simple example that relies on the default values.
We want to bias some CV in the range of 0 to 4.
The wavelets will therefore be scaled to match that range.
Using Db8 wavelets this results in 30 basis functions (including the constant one), with their starting points given by \f$ -14 \frac{4}{15}, -13 \frac{4}{15}, \cdots , 0 , \cdots, 13 \frac{4}{15}, 14 \frac{4}{15} \f$.
\plumedfile
BF_DB_WAVELETS ...
 ORDER=8
 MINIMUM=0.0
 MAXIMUM=4.0
 LABEL=bf
... BF_DB_WAVELETS
\endplumedfile


By omitting wavelets with only insignificant parts, we can reduce the number of basis functions, e.g. a threshold of 0.01 will remove the 8 leftmost shifts.
\plumedfile
BF_DB_WAVELETS ...
 ORDER=8
 MINIMUM=0.0
 MAXIMUM=4.0
 TAILS_THRESHOLD=0.01
 LABEL=bf
... BF_DB_WAVELETS
\endplumedfile


The length of the individual basis functions can also be adjusted to fit the specific problem.
If for example the wavelets are instead scaled to length 3, there will be 35 basis functions, with leftmost points at \f$ -14 \frac{3}{15}, -13 \frac{3}{15}, \cdots, 0, \cdots, 18 \frac{3}{15}, 19 \frac{3}{15} \f$.

\plumedfile
BF_DB_WAVELETS ...
 ORDER=8
 MINIMUM=0.0
 MAXIMUM=4.0
 FUNCTION_LENGTH=3
 LABEL=bf
... BF_DB_WAVELETS
\endplumedfile

*/
//+ENDPLUMEDOC


class BF_DbWavelets : public BasisFunctions {
  // Grid that holds the Wavelet values and its derivative
  std::unique_ptr<Grid> waveletGrid_;
  bool use_mother_wavelet_;
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
  keys.add("optional","GRID_SIZE","The number of grid bins of the Wavelet function. Because of the used construction algorithm this value definess the minimum number, while the true number will probably be larger. Defaults to 1000.");
  keys.add("optional","FUNCTION_LENGTH","The length of the support of the scaled basis functions. This can be used to alter the scaling of the basis functions. Is by default set to the total size of the interval. This also influences the number of actually used basis functions, as all shifted functions that are partially supported in the CV space are used.");
  keys.add("optional","TAILS_THRESHOLD","The threshold for cutting off tail wavelets with respect to the maximum value. All shifted wavelet functions that will only have values lower below the threshold in the CV space will be excluded from the basis set. Defaults to 0 (include all).");
  keys.addFlag("MOTHER_WAVELET", false, "If this flag is set the \"true\" wavelet function (mother wavelet) will be used instead of the scaling function (father wavelet). Makes only sense for multiresolution, which is at the moment not implemented.");
  keys.addFlag("DUMP_WAVELET_GRID", false, "If this flag is set the grid with the wavelet values will be written to a file called \"wavelet_grid.data\".");
  // why is this removed?
  keys.remove("NUMERICAL_INTEGRALS");
}


BF_DbWavelets::BF_DbWavelets(const ActionOptions& ao):
  PLUMED_VES_BASISFUNCTIONS_INIT(ao),
  use_mother_wavelet_(false)
{

  // parse grid properties and set it up
  parseFlag("MOTHER_WAVELET", use_mother_wavelet_);
  unsigned gridsize = 1000;
  parse("GRID_SIZE", gridsize);
  waveletGrid_ = DbWaveletGrid::setupGrid(getOrder(), gridsize, use_mother_wavelet_);
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
  std::string type; use_mother_wavelet_ ? type = "Mother Wavelet" : type = "Father Wavelet";
  setType(type);
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
