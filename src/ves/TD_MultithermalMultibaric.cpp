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

#include "TargetDistribution.h"
#include "GridIntegrationWeights.h"
#include "core/ActionRegister.h"
#include "tools/Grid.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include <cfloat>


namespace PLMD {
namespace ves {

//+PLUMEDOC VES_TARGETDIST TD_MULTITHERMAL_MULTIBARIC
/*
Multithermal-multibaric target distribution (dynamic).

Use the target distribution to sample the multithermal-multibaric ensemble \cite Piaggi-PRL-2019 \cite Okumura-CPL-2004.
In this way, in a single molecular dynamics simulation one can obtain information
about the simulated system in a range of temperatures and pressures.
This range is determined through the keywords MIN_TEMP, MAX_TEMP, MIN_PRESSURE, and MAX_PRESSURE.
One should also specified the target pressure of the barostat with the keyword PRESSURE.

The collective variables (CVs) used to construct the bias potential must be:
  1. the potential energy and the volume or,
  2. the potential energy, the volume, and an order parameter.

Other choices of CVs or a different order of the above mentioned CVs are nonsensical.
The third CV, the order parameter, must be used when the region of the phase diagram under study is crossed by a first order phase transition \cite Piaggi-arXiv-2019 .

The algorithm will explore the free energy at each temperature and pressure up to a predefined free
 energy threshold \f$\epsilon\f$ specified through the keyword THRESHOLD (in kT units).
If only the energy and the volume are being biased, i.e. no phase transition is considered, then THRESHOLD can be set to 1.
If also an order parameter is used then the THRESHOLD should be greater than the barrier for the transformation in kT.
For small systems undergoing a freezing transition THRESHOLD is typically between 20 and 50.

It is also important to specify the number of intermediate temperatures and pressures to consider.
This is done through the keywords STEPS_TEMP and STEPS_PRESSURE.
If the number of intermediate temperature and pressures is too small, then holes might appear in the target distribution.
If it is too large, the performance will deteriorate with no additional advantage.

We now describe the algorithm more rigurously.
The target distribution is given by
\f[
p(E,\mathcal{V},s)=
  \begin{cases}
    1/\Omega_{E,\mathcal{V},s} & \text{if there is at least one } \beta',P' \text{ such} \\
             & \text{that } \beta' F_{\beta',P'}(E,\mathcal{V},s)<\epsilon \text{ with}  \\
             & \beta_1>\beta'>\beta_2 \text{ and } P_1<P'<P_2 \\
    0 & \text{otherwise}
  \end{cases}
\f]
with \f$F_{\beta',P'}(E,\mathcal{V},s)\f$ the free energy as a function of energy \f$E\f$ and volume \f$\mathcal{V}\f$ (and optionally the order parameter \f$s\f$) at temperature \f$\beta'\f$ and pressure \f$P'\f$, \f$\Omega_{E,\mathcal{V},s}\f$ is a normalization constant, and \f$\epsilon\f$ is the THRESHOLD.
In practice the condition \f$\beta' F_{\beta',P'}(E,\mathcal{V},s)<\epsilon\f$  is checked in equally spaced points in each dimension \f$\beta'\f$ and \f$P'\f$.
The number of points is determined with the keywords STEPS_TEMP and STEPS_PRESSURE.

Much like in the Wang-Landau algorithm \cite wanglandau or in the multicanonical ensemble \cite Berg-PRL-1992 , a flat histogram is targeted.
The idea behind this choice of target distribution is that all regions of potential energy and volume (and optionally order parameter) that are relevant at all temperatures \f$\beta_1<\beta'<\beta_2\f$ and pressure \f$P_1<P'<P_2\f$ are included in the distribution.

The free energy at temperature \f$\beta'\f$ and pressure \f$P'\f$ is calculated from the free energy at \f$\beta\f$ and \f$P\f$ using:
\f[
\beta' F_{\beta',P'}(E,\mathcal{V},s) = \beta F_{\beta,P}(E,\mathcal{V},s) + (\beta' - \beta) E + (\beta' P' - \beta P ) \mathcal{V} + C
\f]
with \f$C\f$ such that \f$F_{\beta',P'}(E_m,\mathcal{V}_m,s_m)=0\f$ with \f$E_{m},\mathcal{V}_m,s_m\f$ the position of the free energy minimum.
\f$ \beta F_{\beta,P}(E,\mathcal{V},s) \f$ is not know from the start and is instead found during the simulation.
Therefore \f$ p(E,\mathcal{V},s) \f$ is determined iteratively as done in the well tempered target distribution \cite Valsson-JCTC-2015.

The output of these simulations can be reweighted in order to obtain information at all temperatures and pressures in the targeted region of TP plane.
The reweighting can be performed using the action \ref REWEIGHT_TEMP_PRESS.

The multicanonical ensemble (fixed volume) can be targeted using \ref TD_MULTICANONICAL.

\par Examples

The following input can be used to run a simulation in the multithermal-multibaric ensemble.
The region of the temperature-pressure plane that will be explored is 260-350 K and 1 bar- 300 MPa.
The energy and the volume are used as collective variables.
Legendre polynomials are used to construct the two dimensional bias potential.
The averaged stochastic gradient descent algorithm is chosen to optimize the VES functional.
The target distribution is updated every 100 optimization steps (200 ps here) using the last estimation of the free energy.

\plumedfile
# Use energy and volume as CVs
energy: ENERGY
vol: VOLUME

# Basis functions
bf1: BF_LEGENDRE ORDER=10 MINIMUM=-14750 MAXIMUM=-12250
bf2: BF_LEGENDRE ORDER=10 MINIMUM=6.5 MAXIMUM=8.25

# Target distribution
TD_MULTITHERMAL_MULTIBARIC ...
 MIN_TEMP=260
 MAX_TEMP=350
 MAX_PRESSURE=180.66422571 # 300 MPa
 MIN_PRESSURE=0.06022140857 # 1 bar
 PRESSURE=0.06022140857 # 1 bar
 STEPS_PRESSURE=20
 STEPS_TEMP=20
 SIGMA=50.,0.05
 THRESHOLD=1
 LABEL=td_multi
... TD_MULTITHERMAL_MULTIBARIC

# Bias expansion
VES_LINEAR_EXPANSION ...
 ARG=energy,vol
 BASIS_FUNCTIONS=bf1,bf2
 TEMP=300.0
 GRID_BINS=200,200
 TARGET_DISTRIBUTION=td_multi
 LABEL=b1
... VES_LINEAR_EXPANSION

# Optimization algorithm
OPT_AVERAGED_SGD ...
  BIAS=b1
  STRIDE=500
  LABEL=o1
  STEPSIZE=1.0
  FES_OUTPUT=500
  BIAS_OUTPUT=500
  TARGETDIST_OUTPUT=500
  COEFFS_OUTPUT=100
  TARGETDIST_STRIDE=100
... OPT_AVERAGED_SGD

\endplumedfile


The multithermal-multibaric target distribution can also be used to explore regions of the phase diagram crossed by first order phase transitions.
Consider a system of 250 atoms that crystallizes in the fcc crystal structure.
The region of the temperature-pressure plane that will be explored is 350-450 K and 1bar-1GPa.
We assume that inside this region we can find the liquid-fcc coexistence line that we would like to obtain.
In this case in addition to the energy and volume, an order parameter must also be biased.
The energy, volume, and an order parameter are used as collective variables to construct the bias potential.
We choose as order parameter the \ref FCCUBIC.
Legendre polynomials are used to construct the three dimensional bias potential.
The averaged stochastic gradient descent algorithm is chosen to optimize the VES functional.
The target distribution is updated every 100 optimization steps (200 ps here) using the last estimation of the free energy.

\plumedfile
# Use energy, volume and FCCUBIC as CVs
energy: ENERGY
vol: VOLUME
fcc: FCCUBIC SPECIES=1-256 SWITCH={CUBIC D_0=0.4 D_MAX=0.5} MORE_THAN={RATIONAL R_0=0.45 NN=12 MM=24}

# Basis functions
bf1: BF_LEGENDRE ORDER=8 MINIMUM=-26500 MAXIMUM=-23500
bf2: BF_LEGENDRE ORDER=8 MINIMUM=8.0 MAXIMUM=11.5
bf3: BF_LEGENDRE ORDER=8 MINIMUM=0.0 MAXIMUM=256.0

# Target distribution
TD_MULTITHERMAL_MULTIBARIC ...
 LABEL=td_multitp
 MIN_TEMP=350.0
 MAX_TEMP=450.0
 MIN_PRESSURE=0.06022140857
 MAX_PRESSURE=602.2140857
 PRESSURE=301.10704285
 SIGMA=250.0,0.1,10.0
 THRESHOLD=15
 STEPS_TEMP=20
 STEPS_PRESSURE=20
... TD_MULTITHERMAL_MULTIBARIC

# Expansion
VES_LINEAR_EXPANSION ...
 ARG=energy,vol,fcc.morethan
 BASIS_FUNCTIONS=bf1,bf2,bf3
 TEMP=400.0
 GRID_BINS=40,40,40
 TARGET_DISTRIBUTION=td_multitp
 LABEL=b1
... VES_LINEAR_EXPANSION

# Optimization algorithm
OPT_AVERAGED_SGD ...
  BIAS=b1
  STRIDE=500
  LABEL=o1
  STEPSIZE=1.0
  FES_OUTPUT=500
  BIAS_OUTPUT=500
  TARGETDIST_OUTPUT=500
  COEFFS_OUTPUT=100
  TARGETDIST_STRIDE=500
... OPT_AVERAGED_SGD

\endplumedfile

*/
//+ENDPLUMEDOC

class TD_MultithermalMultibaric: public TargetDistribution {
private:
  double threshold_, min_temp_, max_temp_;
  double min_press_, max_press_, press_;
  std::vector<double> sigma_;
  unsigned steps_temp_, steps_pressure_;
public:
  static void registerKeywords(Keywords&);
  explicit TD_MultithermalMultibaric(const ActionOptions& ao);
  void updateGrid();
  double getValue(const std::vector<double>&) const;
  ~TD_MultithermalMultibaric() {}
  double GaussianSwitchingFunc(const double, const double, const double) const;
};


PLUMED_REGISTER_ACTION(TD_MultithermalMultibaric,"TD_MULTITHERMAL_MULTIBARIC")


void TD_MultithermalMultibaric::registerKeywords(Keywords& keys) {
  TargetDistribution::registerKeywords(keys);
  keys.add("compulsory","THRESHOLD","1","Maximum exploration free energy in kT.");
  keys.add("compulsory","MIN_TEMP","Minimum energy.");
  keys.add("compulsory","MAX_TEMP","Maximum energy.");
  keys.add("compulsory","MIN_PRESSURE","Minimum pressure.");
  keys.add("compulsory","MAX_PRESSURE","Maximum pressure.");
  keys.add("compulsory","PRESSURE","Target pressure of the barostat used in the MD engine.");
  keys.add("compulsory","STEPS_TEMP","20","Number of temperature steps.");
  keys.add("compulsory","STEPS_PRESSURE","20","Number of pressure steps.");
  keys.add("optional","SIGMA","The standard deviation parameters of the Gaussian kernels used for smoothing the target distribution. One value must be specified for each argument, i.e. one value per CV. A value of 0.0 means that no smooting is performed, this is the default behavior.");
}


TD_MultithermalMultibaric::TD_MultithermalMultibaric(const ActionOptions& ao):
  PLUMED_VES_TARGETDISTRIBUTION_INIT(ao),
  threshold_(1.0),
  min_temp_(0.0),
  max_temp_(1000.0),
  min_press_(0.0),
  max_press_(1000.0),
  steps_temp_(20),
  steps_pressure_(20),
  sigma_(0.0)
{
  log.printf("  Multithermal-multibaric target distribution");
  log.printf("\n");

  log.printf("  Please read and cite ");
  log << plumed.cite("Piaggi and Parrinello, Phys. Rev. Lett. 122 (5), 050601 (2019)");
  log.printf(" and ");
  log << plumed.cite("Piaggi and Parrinello, arXiv preprint arXiv:1904.05624 (2019)");
  log.printf("\n");


  parse("THRESHOLD",threshold_);
  if(threshold_<=0.0) {
    plumed_merror("TD_MULTITHERMAL_MULTIBARIC target distribution: the value of the threshold should be positive.");
  }
  parse("MIN_TEMP",min_temp_);
  parse("MAX_TEMP",max_temp_);
  parse("MIN_PRESSURE",min_press_);
  parse("MAX_PRESSURE",max_press_);
  parse("PRESSURE",press_);
  parseVector("SIGMA",sigma_);
  if(sigma_.size()<2 || sigma_.size()>3) plumed_merror(getName()+": SIGMA takes 2 or 3 values as input.");
  parse("STEPS_TEMP",steps_temp_);
  parse("STEPS_PRESSURE",steps_pressure_);
  steps_temp_ += 1;
  steps_pressure_ += 1;

  setDynamic();
  setFesGridNeeded();
  checkRead();
}


double TD_MultithermalMultibaric::getValue(const std::vector<double>& argument) const {
  plumed_merror("getValue not implemented for TD_MultithermalMultibaric");
  return 0.0;
}


void TD_MultithermalMultibaric::updateGrid() {
  if (getStep() == 0) {
    if(targetDistGrid().getDimension()>3 && targetDistGrid().getDimension()<2) plumed_merror(getName()+" works only with 2 or 3 arguments, i.e. energy and volume, or energy, volume, and CV");
    if(sigma_.size()!=targetDistGrid().getDimension()) plumed_merror(getName()+": mismatch between SIGMA dimension and number of arguments");
    // Use uniform TD
    std::vector<double> integration_weights = GridIntegrationWeights::getIntegrationWeights(getTargetDistGridPntr());
    double norm = 0.0;
    for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++) {
      double value = 1.0;
      norm += integration_weights[l]*value;
      targetDistGrid().setValue(l,value);
    }
    targetDistGrid().scaleAllValuesAndDerivatives(1.0/norm);
    logTargetDistGrid().setMinToZero();
  } else {
    double beta = getBeta();
    double beta_prime_min = 1./(plumed.getAtoms().getKBoltzmann()*min_temp_);
    double beta_prime_max = 1./(plumed.getAtoms().getKBoltzmann()*max_temp_);
    plumed_massert(getFesGridPntr()!=NULL,"the FES grid has to be linked to use TD_MultithermalMultibaric!");
    // Set all to zero
    for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++) {
      double value = 0.0;
      targetDistGrid().setValue(l,value);
    }
    // Loop over pressures and temperatures
    for(unsigned i=0; i<steps_temp_; i++) {
      double beta_prime=beta_prime_min + (beta_prime_max-beta_prime_min)*i/(steps_temp_-1);
      for(unsigned j=0; j<steps_pressure_; j++) {
        double pressure_prime=min_press_ + (max_press_-min_press_)*j/(steps_pressure_-1);
        // Find minimum for this pressure and temperature
        double minval=DBL_MAX;
        for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++) {
          double energy = targetDistGrid().getPoint(l)[0];
          double volume = targetDistGrid().getPoint(l)[1];
          double value = getFesGridPntr()->getValue(l);
          value = beta*value + (beta_prime-beta)*energy + (beta_prime*pressure_prime-beta*press_)*volume;
          if(value<minval) {
            minval=value;
          }
        }
        // Now check which energies and volumes are below X kt
        for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++) {
          double energy = targetDistGrid().getPoint(l)[0];
          double volume = targetDistGrid().getPoint(l)[1];
          double value = getFesGridPntr()->getValue(l);
          value = beta*value + (beta_prime-beta)*energy + (beta_prime*pressure_prime-beta*press_)*volume - minval;
          if (value<threshold_) {
            double value = 1.0;
            targetDistGrid().setValue(l,value);
          }
        }
      }
    }
    std::vector<unsigned> nbin=targetDistGrid().getNbin();
    std::vector<double> dx=targetDistGrid().getDx();
    unsigned dim=targetDistGrid().getDimension();
    // Smoothening
    for(Grid::index_t index=0; index<targetDistGrid().getSize(); index++) {
      std::vector<unsigned> indices = targetDistGrid().getIndices(index);
      std::vector<double> point = targetDistGrid().getPoint(index);
      double value = targetDistGrid().getValue(index);
      if (value>(1-1.e-5)) { // Apply only if this grid point was 1.
        // Apply gaussians around
        std::vector<int> minBin(dim), maxBin(dim); // These cannot be unsigned
        // Only consider contributions less than n*sigma bins apart from the actual distance
        for(unsigned k=0; k<dim; k++) {
          int deltaBin=std::floor(5*sigma_[k]/dx[k]);
          minBin[k]=indices[k] - deltaBin;
          if (minBin[k] < 0) minBin[k]=0;
          if (minBin[k] > (nbin[k]-1)) minBin[k]=nbin[k]-1;
          maxBin[k]=indices[k] + deltaBin;
          if (maxBin[k] > (nbin[k]-1)) maxBin[k]=nbin[k]-1;
        }
        if (dim==2) {
          for(unsigned l=minBin[0]; l<maxBin[0]+1; l++) {
            for(unsigned m=minBin[1]; m<maxBin[1]+1; m++) {
              std::vector<unsigned> indices_prime(dim);
              indices_prime[0]=l;
              indices_prime[1]=m;
              Grid::index_t index_prime = targetDistGrid().getIndex(indices_prime);
              std::vector<double> point_prime = targetDistGrid().getPoint(index_prime);
              double value_prime = targetDistGrid().getValue(index_prime);
              // Apply gaussian
              double gaussian_value = 1;
              for(unsigned k=0; k<dim; k++) {
                gaussian_value *= GaussianSwitchingFunc(point_prime[k],point[k],sigma_[k]);
              }
              if (value_prime<gaussian_value) {
                targetDistGrid().setValue(index_prime,gaussian_value);
              }
            }
          }
        } else if (dim==3) {
          for(unsigned l=minBin[0]; l<maxBin[0]+1; l++) {
            for(unsigned m=minBin[1]; m<maxBin[1]+1; m++) {
              for(unsigned n=minBin[2]; n<maxBin[2]+1; n++) {
                std::vector<unsigned> indices_prime(dim);
                indices_prime[0]=l;
                indices_prime[1]=m;
                indices_prime[2]=n;
                Grid::index_t index_prime = targetDistGrid().getIndex(indices_prime);
                std::vector<double> point_prime = targetDistGrid().getPoint(index_prime);
                double value_prime = targetDistGrid().getValue(index_prime);
                // Apply gaussian
                double gaussian_value = 1;
                for(unsigned k=0; k<dim; k++) {
                  gaussian_value *= GaussianSwitchingFunc(point_prime[k],point[k],sigma_[k]);
                }
                if (value_prime<gaussian_value) {
                  targetDistGrid().setValue(index_prime,gaussian_value);
                }
              }
            }
          }
        }
      }
    }
    // Normalize
    std::vector<double> integration_weights = GridIntegrationWeights::getIntegrationWeights(getTargetDistGridPntr());
    double norm = 0.0;
    for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++) {
      double value = targetDistGrid().getValue(l);
      norm += integration_weights[l]*value;
    }
    targetDistGrid().scaleAllValuesAndDerivatives(1.0/norm);
    logTargetDistGrid().setMinToZero();
  }
}

inline
double TD_MultithermalMultibaric::GaussianSwitchingFunc(const double argument, const double center, const double sigma) const {
  if(sigma>0.0) {
    double arg=(argument-center)/sigma;
    return exp(-0.5*arg*arg);
  }
  else {
    return 0.0;
  }
}


}
}
