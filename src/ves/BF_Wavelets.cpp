/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2021 The VES code team
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
#include "GridLinearInterpolation.h"
#include "tools/Grid.h"
#include "VesTools.h"
#include "WaveletGrid.h"
#include "core/ActionRegister.h"
#include "tools/Exception.h"
#include "core/PlumedMain.h"


namespace PLMD {
namespace ves {


//+PLUMEDOC VES_BASISF BF_WAVELETS
/*
Daubechies Wavelets basis functions.

Note: at the moment only bases with a single level of scaling functions are usable, as multiscale optimization is not yet implemented.

This basis set uses Daubechies Wavelets \cite daubechies_ten_1992 to construct a complete and orthogonal basis. See \cite ValssonPampel_Wavelets_2022 for full details.

The basis set is based on using a pair of functions, the scaling function (or father wavelet) \f$\phi\f$ and the wavelet function (or mother wavelet) \f$\psi\f$.
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
Currently the "maximum phase" type (Db) and the "least asymmetric" (Sym) type are implemented.
We recommend to use Symlets.

As an example two adjacent basis functions of both Sym8 (ORDER=8, TYPE=SYMLET) and Db8 (ORDER=8, TYPE=DAUBECHIES) is shown in the figure.
The full basis consists of shifted wavelets in the full specified interval.

\image html ves_basisf-wavelets.png


\par Specify the wavelet type

The TYPE keyword sets the type of Wavelet, at the moment "DAUBECHIES" and "SYMLETS" are available.
The specified ORDER of the basis corresponds to the number of vanishing moments of the wavelet, i.e. if TYPE was specified as "DAUBECHIES" an order of 8 results in Db8 wavelets.


\par Specify the number of functions

The resulting basis set consists of integer shifts of the wavelet with some scaling \f$j\f$,
\f[
  V(x) = \sum_i \alpha_i * \phi_i (x) = \sum_i \alpha_i * \phi(\frac{x+i}{j})
\f]
with the variational parameters \f$ \alpha \f$.
Additionally a constant basis function is included.

There are two different ways to specify the number of used basis functions implemented.
You can either specify the scale or alternatively a fixed number of basis function.

Coming from the multiresolution aspect of wavelets, you can set the scale of the father wavelets, i.e. the largest scale used for approximation.
This can be done with the FUNCTION_LENGTH keyword.
It should be given in the same units as the used CV and specifies the length (of the domain interval) of the individual father wavelet functions.

Alternatively a fixed number of basis functions for the bias expansion can be specified with the NUM_BF keyword, which will set the scale automatically to match the desired number of functions.
Note that this also includes the constant function.

If you do not specify anything, it is assumed that the range of the bias should match the scale of the wavelet functions.
More precise, the basis functions are scaled to match the specified size of the CV space (MINIMUM and MAXIMUM keywords).
This has so far been a good initial choice.

If the wavelets are scaled to match the CV range exactly there would be \f$4*\text{ORDER} -3\f$ basis functions whose domain is at least partially in this region.
This number is adjusted if FUNCTION_LENGTH or NUM_BF is specified.
Additionally, some of the shifted basis functions will not have significant contributions because of their function values being close to zero over the full range of the bias.
These 'tail wavelets' can be omitted by using the TAILS_THRESHOLD keyword.
This omits all shifted functions that have only function values smaller than a fraction of their maximum value inside the bias range.
Using a value of e.g. 0.01 will already reduce the number of basis functions significantly.
The default setting will not omit any tail wavelets (i.e. TAILS_THRESHOLD=0).

The number of basis functions is then not easily determinable a priori but will be given in the logfile.
Additionally the starting point (leftmost defined point) of the individual basis functions is printed.


With the PERIODIC keyword the basis set can also be used to bias periodic CVs.
Then the shift between the functions will be chosen such that the function at the left border and right border coincide.
If the FUNCTION_LENGTH keyword is used together with PERIODIC, a smaller length might be chosen to satisfy this requirement.


\par Grid

The values of the wavelet function are generated on a grid.
Using the cascade algorithm results in doubling the grid values for each iteration.
This means that the grid size will always be a power of two multiplied by the number of coefficients (\f$ 2*\text{ORDER} -1\f$) for the specified wavelet.
Using the MIN_GRID_SIZE keyword a lower bound for the number of grid points can be specified.
By default at least 1,000 grid points are used.
Function values in between grid points are calculated by linear interpolation.

\par Optimization notes

To avoid 'blind' optimization of the basis functions outside the currently sampled area, it is often beneficial to use the OPTIMIZATION_THRESHOLD keyword of the \ref VES_LINEAR_EXPANSION (set it to a small value, e.g. 1e-6)

\par Examples


First a very simple example that relies on the default values.
We want to bias some CV in the range of 0 to 4.
The wavelets will therefore be scaled to match that range.
Using Db8 wavelets this results in 30 basis functions (including the constant one), with their starting points given by \f$ -14*\frac{4}{15}, -13*\frac{4}{15}, \cdots , 0 , \cdots, 13*\frac{4}{15}, 14*\frac{4}{15} \f$.
\plumedfile
BF_WAVELETS ...
 ORDER=8
 TYPE=DAUBECHIES
 MINIMUM=0.0
 MAXIMUM=4.0
 LABEL=bf
... BF_WAVELETS
\endplumedfile


By omitting wavelets with only insignificant parts, we can reduce the number of basis functions. Using a threshold of 0.01 will in this example remove the 8 leftmost shifts, which we can check in the logfile.
\plumedfile
BF_WAVELETS ...
 ORDER=8
 TYPE=DAUBECHIES
 MINIMUM=0.0
 MAXIMUM=4.0
 TAILS_THRESHOLD=0.01
 LABEL=bf
... BF_WAVELETS
\endplumedfile


The length of the individual basis functions can also be adjusted to fit the specific problem.
If for example the wavelets are instead scaled to length 3, there will be 35 basis functions, with leftmost points at \f$ -14*\frac{3}{15}, -13*\frac{3}{15}, \cdots, 0, \cdots, 18*\frac{3}{15}, 19*\frac{3}{15} \f$.
\plumedfile
BF_WAVELETS ...
 ORDER=8
 TYPE=DAUBECHIES
 MINIMUM=0.0
 MAXIMUM=4.0
 FUNCTION_LENGTH=3
 LABEL=bf
... BF_WAVELETS
\endplumedfile


Alternatively you can also specify the number of basis functions. Here we specify the usage of 40 Sym10 wavelet functions. We also used a custom minimum size for the grid and want it to be printed to a file with a specific numerical format.
\plumedfile
BF_WAVELETS ...
 ORDER=10
 TYPE=SYMLETS
 MINIMUM=0.0
 MAXIMUM=4.0
 NUM_BF=40
 MIN_GRID_SIZE=500
 DUMP_WAVELET_GRID
 WAVELET_FILE_FMT=%11.4f
 LABEL=bf
... BF_WAVELETS
\endplumedfile

*/
//+ENDPLUMEDOC


class BF_Wavelets : public BasisFunctions {
private:
  void setupLabels() override;
  /// ptr to Grid that holds the Wavelet values and its derivative
  std::unique_ptr<Grid> waveletGrid_;
  /// calculate threshold for omitted tail wavelets
  std::vector<double> getCutoffPoints(const double& threshold);
  /// scale factor of the individual BFs to match specified length
  double scale_;
  /// shift of the individual BFs
  std::vector<double> shifts_;
public:
  static void registerKeywords( Keywords&);
  explicit BF_Wavelets(const ActionOptions&);
  void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const override;
};


PLUMED_REGISTER_ACTION(BF_Wavelets,"BF_WAVELETS")


void BF_Wavelets::registerKeywords(Keywords& keys) {
  BasisFunctions::registerKeywords(keys);
  keys.add("compulsory","TYPE","Specify the wavelet type. Currently available are DAUBECHIES Wavelets with minimum phase and the more symmetric SYMLETS");
  keys.add("optional","FUNCTION_LENGTH","The domain size of the individual basis functions. (length) This is used to alter the scaling of the basis functions. By default it is set to the total size of the interval. This also influences the number of actually used basis functions, as all shifted functions that are partially supported in the CV space are used.");
  keys.add("optional","NUM_BF","The number of basis functions that should be used. Includes the constant one and N-1 shifted wavelets within the specified range. Cannot be used together with FUNCTION_LENGTH.");
  keys.add("optional","TAILS_THRESHOLD","The threshold for cutting off tail wavelets as a fraction of the maximum value. All shifted wavelet functions that only have values smaller than the threshold in the bias range will be excluded from the basis set. Defaults to 0 (include all).");
  keys.addFlag("MOTHER_WAVELET", false, "If this flag is set mother wavelets will be used instead of the scaling function (father wavelet). Makes only sense for multiresolution, which is at the moment not usable.");
  keys.add("optional","MIN_GRID_SIZE","The minimal number of grid bins of the Wavelet function. The true number depends also on the used wavelet type and will probably be larger. Defaults to 1000.");
  keys.addFlag("DUMP_WAVELET_GRID", false, "If this flag is set the grid with the wavelet values will be written to a file.  This file is called wavelet_grid.data.");
  keys.add("optional","WAVELET_FILE_FMT","The number format of the wavelet grid values and derivatives written to file. By default it is %15.8f.\n");
  keys.addFlag("PERIODIC", false, "Use periodic version of basis set.");
  keys.remove("NUMERICAL_INTEGRALS");
}


BF_Wavelets::BF_Wavelets(const ActionOptions& ao):
  PLUMED_VES_BASISFUNCTIONS_INIT(ao),
  waveletGrid_(nullptr),
  scale_(0.0)
{
  log.printf("  Wavelet basis functions, see and cite ");
  log << plumed.cite("Pampel and Valsson, J. Chem. Theory Comput. 18, 4127-4141 (2022) - DOI:10.1021/acs.jctc.2c00197");

  // parse properties for waveletGrid and set it up
  bool use_mother_wavelet;
  parseFlag("MOTHER_WAVELET", use_mother_wavelet);

  std::string wavelet_type_str;
  parse("TYPE", wavelet_type_str);
  addKeywordToList("TYPE", wavelet_type_str);

  unsigned min_grid_size = 1000;
  parse("MIN_GRID_SIZE", min_grid_size);
  if(min_grid_size != 1000) {addKeywordToList("MIN_GRID_SIZE",min_grid_size);}

  waveletGrid_ = WaveletGrid::setupGrid(getOrder(), min_grid_size, use_mother_wavelet, WaveletGrid::stringToType(wavelet_type_str));
  bool dump_wavelet_grid=false;
  parseFlag("DUMP_WAVELET_GRID", dump_wavelet_grid);
  if (dump_wavelet_grid) {
    OFile wavelet_gridfile;
    std::string fmt = "%13.6f";
    parse("WAVELET_FILE_FMT",fmt);
    waveletGrid_->setOutputFmt(fmt); // property of grid not OFile determines fmt
    wavelet_gridfile.link(*this);
    wavelet_gridfile.enforceBackup();
    wavelet_gridfile.open(getLabel()+".wavelet_grid.data");
    waveletGrid_->writeToFile(wavelet_gridfile);
  }

  bool periodic = false;
  parseFlag("PERIODIC",periodic);
  if (periodic) {addKeywordToList("PERIODIC",periodic);}

  // now set up properties of basis set
  unsigned intrinsic_length = 2*getOrder() - 1; // length of unscaled wavelet
  double bias_length = intervalMax() - intervalMin(); // intervalRange() is not yet set

  // parse threshold for tail wavelets and get respective cutoff points
  double threshold = 0.0;
  std::vector<double> cutoffpoints (2);
  parse("TAILS_THRESHOLD",threshold);
  plumed_massert(threshold < 1, "TAILS_THRESHOLD should be significantly smaller than 1.");
  if(threshold == 0.0) {
    cutoffpoints = {0.0, static_cast<double>(intrinsic_length)};
  }
  else {
    plumed_massert(!periodic, "TAILS_THRESHOLD can't be used with the periodic wavelet variant");
    addKeywordToList("TAILS_THRESHOLD",threshold);
    cutoffpoints = getCutoffPoints(threshold);
  };

  double function_length = bias_length;
  parse("FUNCTION_LENGTH",function_length);
  if(function_length != bias_length) {
    if (periodic) {  // shifted functions need to fit into interval exactly -> reduce size if not
      unsigned num_shifts = ceil(bias_length * intrinsic_length / function_length);
      function_length = bias_length * intrinsic_length / num_shifts;
    }
    addKeywordToList("FUNCTION_LENGTH",function_length);
  }

  // determine number of BFs and needed scaling
  unsigned num_BFs = 0;
  parse("NUM_BF",num_BFs);
  if(num_BFs == 0) { // get from function length
    scale_ = intrinsic_length / function_length;
    if (periodic) {
      // this is the same value as num_shifts above + constant
      num_BFs = static_cast<unsigned>(bias_length * scale_) + 1;
    }
    else {
      num_BFs = 1; // constant one
      // left shifts (w/o left cutoff) + right shifts - right cutoff - 1
      num_BFs += static_cast<unsigned>(ceil(cutoffpoints[1] + (bias_length)*scale_ - cutoffpoints[0]) - 1);
    }
  }
  else {
    plumed_massert(num_BFs > 0, "The number of basis functions has to be positive (NUM_BF > 0)");
    // check does not work if function length was given as intrinsic length, but can't check for keyword use directly
    plumed_massert(function_length==bias_length,"The keywords \"NUM_BF\" and \"FUNCTION_LENGTH\" cannot be used at the same time");
    addKeywordToList("NUM_BF",num_BFs);

    if (periodic) {  // inverted num_BFs calculation from where FUNCTION_LENGTH is specified
      scale_ = (num_BFs  - 1) / bias_length ;
    }
    else {
      double cutoff_length = cutoffpoints[1] - cutoffpoints [0];
      double intrinsic_bias_length = num_BFs - cutoff_length + 1; // length of bias in intrinsic scale of wavelets
      scale_ = intrinsic_bias_length / bias_length;
    }
  }

  setNumberOfBasisFunctions(num_BFs);

  // now set up the starting points of the basis functions
  shifts_.push_back(0.0); // constant BF – never used, just for clearer notation
  for(unsigned int i = 1; i < getNumberOfBasisFunctions(); ++i) {
    shifts_.push_back(-intervalMin()*scale_ + cutoffpoints[1] - i);
  }

  // set some properties
  setIntrinsicInterval(0.0,intrinsic_length);
  periodic ? setPeriodic() : setNonPeriodic();
  setIntervalBounded();
  setType(wavelet_type_str);
  setDescription("Wavelets as localized basis functions");
  setupBF();
  checkRead();

  log.printf("  Each basisfunction spans %f in CV space\n", intrinsic_length/scale_);
}


void BF_Wavelets::getAllValues(const double arg, double& argT, bool& inside_range, std::vector<double>& values, std::vector<double>& derivs) const {
  argT=checkIfArgumentInsideInterval(arg,inside_range);
  //
  values[0]=1.0;
  derivs[0]=0.0;
  for(unsigned int i = 1; i < getNumberOfBasisFunctions(); ++i) {
    // scale and shift argument to match current wavelet
    double x = shifts_[i] + argT*scale_;
    if (arePeriodic()) { // periodic interval [0,intervalRange*scale]
      x = x - floor(x/(intervalRange()*scale_))*intervalRange()*scale_;
    }

    if (x < 0 || x >= intrinsicIntervalMax()) { // Wavelets are 0 outside the defined range
      values[i] = 0.0; derivs[i] = 0.0;
    }
    else {
      std::vector<double> temp_deriv (1);
      values[i] = GridLinearInterpolation::getGridValueAndDerivativesWithLinearInterpolation(waveletGrid_.get(), {x}, temp_deriv);
      derivs[i] = temp_deriv[0] * scale_; // scale derivative
    }
  }
  if(!inside_range) {for(auto& deriv : derivs) {deriv=0.0;}}
}


// returns left and right cutoff point of Wavelet
// threshold is a percent value of maximum
std::vector<double> BF_Wavelets::getCutoffPoints(const double& threshold) {
  double threshold_value = threshold * waveletGrid_->getMaxValue();
  std::vector<double> cutoffpoints;

  for (size_t i = 0; i < waveletGrid_->getSize(); ++i) {
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
void BF_Wavelets::setupLabels() {
  setLabel(0,"const");
  for(unsigned int i=1; i < getNumberOfBasisFunctions(); i++) {
    double pos = -shifts_[i]/scale_;
    if (arePeriodic()) {
      pos = pos - floor((pos-intervalMin())/intervalRange())*intervalRange();
    }
    std::string is; Tools::convert(pos, is);
    setLabel(i,"i="+is);
  }
}


}
}
