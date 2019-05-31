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

//+PLUMEDOC VES_TARGETDIST TD_MULTICANONICAL
/*
Multicanonical target distribution (dynamic).

Use the target distribution to sample the multicanonical ensemble \cite Berg-PRL-1992 \cite Piaggi-PRL-2019.
In this way, in a single molecular dynamics simulation one can obtain information about the system in a range of temperatures.
This range is determined through the keywords MIN_TEMP and MAX_TEMP.

The collective variables (CVs) used to construct the bias potential must be:
 1. the energy or,
 2. the energy and an order parameter.

Other choices of CVs or a different order of the above mentioned CVs are nonsensical.
The second CV, the order parameter, must be used when one aims at studying a first order phase transition in the chosen temperature interval \cite Piaggi-arXiv-2019.

The algorithm will explore the free energy at each temperature up to a predefined free
 energy threshold \f$\epsilon\f$ specified through the keyword THRESHOLD (in kT units).
If only the energy is biased, i.e. no phase transition is considered, then TRESHOLD can be set to 1.
If also an order parameter is used then the THRESHOLD should be greater than the barrier for the transformation in kT.
For small systems undergoing a freezing transition THRESHOLD is typically between 20 and 50.

When only the potential energy is used as CV the method is equivalent to the Wang-Landau algorithm \cite wanglandau.
The advantage with respect to Wang-Landau is that instead of sampling the potential energy indiscriminately, an interval is chosen on the fly based on the minimum and maximum targeted temperatures.

The algorithm works as follows.
The target distribution for the potential energy is chosen to be:
\f[
p(E)= \begin{cases}
         \frac{1}{E_2-E_1} & \mathrm{if} \quad E_1<E<E_2 \\
         0 & \mathrm{otherwise}
      \end{cases}
\f]
where the energy limits \f$E_1\f$ and \f$E_2\f$ are yet to be determined.
Clearly the interval \f$E_1–E_2\f$ chosen is related to the interval of temperatures \f$T_1-T_2\f$.
To link these two intervals we make use of the following relation:
\f[
\beta' F_{\beta'}(E) = \beta F_{\beta}(E) + (\beta' - \beta) E + C,
\f]
where \f$F_{\beta}(E)\f$ is determined during the optimization and we shall choose \f$C\f$ such that \f$F_{\beta'}(E_{m})=0\f$ with \f$E_{m}\f$ the position of the free energy minimum.
Using this relation we employ an iterative procedure to find the energy interval.
At iteration \f$k\f$ we have the estimates \f$E_1^k\f$ and \f$E_2^k\f$ for \f$E_1\f$ and \f$E_2\f$, and the target distribution is:
\f[
p^k(E)=\frac{1}{E_2^k-E_1^k} \quad \mathrm{for} \quad E_1^k<E<E_2^k.
\f]
\f$E_1^k\f$ and \f$E_2^k\f$ are obtained from the leftmost solution of \f$\beta_2 F_{\beta_2}^{k-1}(E_1^k)=\epsilon\f$ and the rightmost solution of \f$\beta_1 F_{\beta_1}^{k-1}(E_2^k)=\epsilon\f$.
The procedure is repeated until convergence.
This iterative approach is similar to that in \ref TD_WELLTEMPERED.

The version of this algorithm in which the energy and an order parameter are biased is similar to the one described in \ref TD_MULTITHERMAL_MULTIBARIC.

The output of these simulations can be reweighted in order to obtain information at all temperatures in the targeted temperature interval.
The reweighting can be performed using the action \ref REWEIGHT_TEMP_PRESS.

\par Examples

The following input can be used to run a simulation in the multicanonical ensemble.
The temperature interval to be explored is 400-600 K.
The energy is used as collective variable.
Legendre polynomials are used to construct the bias potential.
The averaged stochastic gradient descent algorithm is chosen to optimize the VES functional.
The target distribution is updated every 100 optimization steps (200 ps here) using the last estimation of the free energy.

\plumedfile
# Use energy and volume as CVs
energy: ENERGY

# Basis functions
bf1: BF_LEGENDRE ORDER=20 MINIMUM=-25000 MAXIMUM=-23500

# Target distributions
TD_MULTICANONICAL ...
 LABEL=td_multi
 SIGMA=50.0
 MIN_TEMP=400
 MAX_TEMP=600
 THRESHOLD=1
... TD_MULTICANONICAL

# Expansion
VES_LINEAR_EXPANSION ...
 ARG=energy
 BASIS_FUNCTIONS=bf1
 TEMP=500.0
 GRID_BINS=1000
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
  COEFFS_OUTPUT=10
  TARGETDIST_STRIDE=100
... OPT_AVERAGED_SGD

\endplumedfile

The multicanonical target distribution can also be used to explore a temperature interval in which a first order phase transitions is observed.

*/
//+ENDPLUMEDOC

class TD_Multicanonical: public TargetDistribution {
private:
  double threshold_, min_temp_, max_temp_;
  std::vector<double> sigma_;
  unsigned steps_temp_;
public:
  static void registerKeywords(Keywords&);
  explicit TD_Multicanonical(const ActionOptions& ao);
  void updateGrid();
  double getValue(const std::vector<double>&) const;
  ~TD_Multicanonical() {}
  double GaussianSwitchingFunc(const double, const double, const double) const;
};


PLUMED_REGISTER_ACTION(TD_Multicanonical,"TD_MULTICANONICAL")


void TD_Multicanonical::registerKeywords(Keywords& keys) {
  TargetDistribution::registerKeywords(keys);
  keys.add("compulsory","THRESHOLD","1","Maximum exploration free energy in kT.");
  keys.add("compulsory","MIN_TEMP","Minimum temperature.");
  keys.add("compulsory","MAX_TEMP","Maximum temperature.");
  keys.add("optional","STEPS_TEMP","Number of temperature steps. Only for the 2D version, i.e. energy and order parameter.");
  keys.add("optional","SIGMA","The standard deviation parameters of the Gaussian kernels used for smoothing the target distribution. One value must be specified for each argument, i.e. one value per CV. A value of 0.0 means that no smooting is performed, this is the default behavior.");
}


TD_Multicanonical::TD_Multicanonical(const ActionOptions& ao):
  PLUMED_VES_TARGETDISTRIBUTION_INIT(ao),
  threshold_(1.0),
  min_temp_(0.0),
  max_temp_(1000.0),
  steps_temp_(20),
  sigma_(0.0)
{
  log.printf("  Multicanonical target distribution");
  log.printf("\n");
  log.printf("  Please read and cite ");
  log << plumed.cite("Piaggi and Parrinello, Phys. Rev. Lett. 122 (5), 050601 (2019)");
  log.printf(" and ");
  log << plumed.cite("Piaggi and Parrinello, arXiv preprint arXiv:1904.05624 (2019)");
  log.printf("\n");
  parse("THRESHOLD",threshold_);
  if(threshold_<=0.0) {
    plumed_merror("TD_MULTICANONICAL target distribution: the value of the threshold should be positive.");
  }

  parse("MIN_TEMP",min_temp_);
  parse("MAX_TEMP",max_temp_);
  parseVector("SIGMA",sigma_);
  if(sigma_.size()<1 || sigma_.size()>2) plumed_merror(getName()+": SIGMA takes 1 or 2 values as input.");
  parse("STEPS_TEMP",steps_temp_); // Only used in the 2D version
  steps_temp_ += 1;

  setDynamic();
  setFesGridNeeded();
  checkRead();
}


double TD_Multicanonical::getValue(const std::vector<double>& argument) const {
  plumed_merror("getValue not implemented for TD_Multicanonical");
  return 0.0;
}


void TD_Multicanonical::updateGrid() {
  if (getStep() == 0) {
    if(targetDistGrid().getDimension()>2 && targetDistGrid().getDimension()<1) plumed_merror(getName()+" works only with 1 or 2 arguments, i.e. energy, or energy and CV");
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
    // Two variants: 1D and 2D
    // This could be done with one variant but the 1D variant is useful for pedagogical purposes.
    if(targetDistGrid().getDimension()==1) {
      // 1D variant: Multicanonical without order parameter
      // In this variant we find the minimum and maximum relevant potential energies.
      // Using this information we construct a uniform target distribution inbetween these two.
      double beta = getBeta();
      double beta_prime_min = 1./(plumed.getAtoms().getKBoltzmann()*min_temp_);
      double beta_prime_max = 1./(plumed.getAtoms().getKBoltzmann()*max_temp_);
      plumed_massert(getFesGridPntr()!=NULL,"the FES grid has to be linked to use TD_Multicanonical!");
      // Find minimum of F(U) at temperature min
      double minval=DBL_MAX;
      Grid::index_t minindex = (targetDistGrid().getSize())/2;
      double minpos = targetDistGrid().getPoint(minindex)[0];
      for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++) {
        double value = getFesGridPntr()->getValue(l);
        double argument = targetDistGrid().getPoint(l)[0];
        value = beta*value + (beta_prime_min-beta)*argument;
        if(value<minval) {
          minval=value;
          minpos=argument;
          minindex=l;
        }
      }
      // Find minimum energy at low temperature
      double minimum_low = minpos;
      for(Grid::index_t l=minindex; l>1; l-=1) {
        double argument = targetDistGrid().getPoint(l)[0];
        double argument_next = targetDistGrid().getPoint(l-1)[0];
        double value = getFesGridPntr()->getValue(l);
        double value_next = getFesGridPntr()->getValue(l-1);
        value = beta*value + (beta_prime_min-beta)*argument - minval;
        value_next = beta*value_next + (beta_prime_min-beta)*argument_next - minval;
        if (value<threshold_ && value_next>threshold_) {
          minimum_low = argument_next;
          break;
        }
      }
      // Find maximum energy at low temperature
      double maximum_low = minpos;
      for(Grid::index_t l=minindex; l<(targetDistGrid().getSize()-1); l++) {
        double argument = targetDistGrid().getPoint(l)[0];
        double argument_next = targetDistGrid().getPoint(l+1)[0];
        double value = getFesGridPntr()->getValue(l);
        double value_next = getFesGridPntr()->getValue(l+1);
        value = beta*value + (beta_prime_min-beta)*argument - minval;
        value_next = beta*value_next + (beta_prime_min-beta)*argument_next - minval;
        if (value<threshold_ && value_next>threshold_) {
          maximum_low = argument_next;
          break;
        }
      }
      // Find minimum of F(U) at temperature max
      minval=DBL_MAX;
      for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++) {
        double value = getFesGridPntr()->getValue(l);
        double argument = targetDistGrid().getPoint(l)[0];
        value = beta*value + (beta_prime_max-beta)*argument;
        if(value<minval) {
          minval=value;
          minpos=argument;
          minindex=l;
        }
      }
      // Find minimum energy at high temperature
      double minimum_high = minpos;
      for(Grid::index_t l=minindex; l>1; l-=1) {
        double argument = targetDistGrid().getPoint(l)[0];
        double argument_next = targetDistGrid().getPoint(l-1)[0];
        double value = getFesGridPntr()->getValue(l);
        double value_next = getFesGridPntr()->getValue(l-1);
        value = beta*value + (beta_prime_max-beta)*argument - minval;
        value_next = beta*value_next + (beta_prime_max-beta)*argument_next - minval;
        if (value<threshold_ && value_next>threshold_) {
          minimum_high = argument_next;
          break;
        }
      }
      // Find maximum energy at high temperature
      double maximum_high = minpos;
      for(Grid::index_t l=minindex; l<(targetDistGrid().getSize()-1); l++) {
        double argument = targetDistGrid().getPoint(l)[0];
        double argument_next = targetDistGrid().getPoint(l+1)[0];
        double value = getFesGridPntr()->getValue(l);
        double value_next = getFesGridPntr()->getValue(l+1);
        value = beta*value + (beta_prime_max-beta)*argument - minval;
        value_next = beta*value_next + (beta_prime_max-beta)*argument_next - minval;
        if (value<threshold_ && value_next>threshold_) {
          maximum_high = argument_next;
          break;
        }
      }
      double minimum = minimum_low;
      if (minimum_high<minimum_low) minimum=minimum_high;
      double maximum = maximum_low;
      if (maximum_high>maximum_low) maximum=maximum_high;
      // Construct uniform TD in the interval between minimum and maximum
      std::vector<double> integration_weights = GridIntegrationWeights::getIntegrationWeights(getTargetDistGridPntr());
      double norm = 0.0;
      for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++) {
        double argument = targetDistGrid().getPoint(l)[0];
        double value = 1.0;
        double tmp;
        if(argument < minimum) {
          tmp = GaussianSwitchingFunc(argument,minimum,sigma_[0]);
        }
        else if(argument > maximum) {
          tmp = GaussianSwitchingFunc(argument,maximum,sigma_[0]);
        }
        else {
          tmp = 1.0;
        }
        value *= tmp;
        norm += integration_weights[l]*value;
        targetDistGrid().setValue(l,value);
      }
      targetDistGrid().scaleAllValuesAndDerivatives(1.0/norm);
      logTargetDistGrid().setMinToZero();
    } else if(targetDistGrid().getDimension()==2) {
      // 2D variant: Multicanonical with order parameter
      // In this variant we find for each temperature the relevant region of potential energy and order parameter.
      // The target distribution will be the union of the relevant regions at all temperatures in the temperature interval.
      double beta = getBeta();
      double beta_prime_min = 1./(plumed.getAtoms().getKBoltzmann()*min_temp_);
      double beta_prime_max = 1./(plumed.getAtoms().getKBoltzmann()*max_temp_);
      plumed_massert(getFesGridPntr()!=NULL,"the FES grid has to be linked to use TD_MulticanonicalWithCV!");
      // Set all to zero
      for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++) {
        double value = 0.0;
        targetDistGrid().setValue(l,value);
      }
      // Loop over temperatures
      for(unsigned i=0; i<steps_temp_; i++) {
        double beta_prime=beta_prime_min + (beta_prime_max-beta_prime_min)*i/(steps_temp_-1);
        // Find minimum for this temperature
        double minval=DBL_MAX;
        for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++) {
          double energy = targetDistGrid().getPoint(l)[0];
          double value = getFesGridPntr()->getValue(l);
          value = beta*value + (beta_prime-beta)*energy;
          if(value<minval) {
            minval=value;
          }
        }
        // Now check which energies and volumes are below X kt
        for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++) {
          double energy = targetDistGrid().getPoint(l)[0];
          double value = getFesGridPntr()->getValue(l);
          value = beta*value + (beta_prime-beta)*energy - minval;
          if (value<threshold_) {
            double value = 1.0;
            targetDistGrid().setValue(l,value);
          }
        }
      }
      std::vector<unsigned> nbin=targetDistGrid().getNbin();
      std::vector<double> dx=targetDistGrid().getDx();
      // Smoothening
      for(unsigned i=0; i<nbin[0]; i++) {
        for(unsigned j=0; j<nbin[1]; j++) {
          std::vector<unsigned> indices(2);
          indices[0]=i;
          indices[1]=j;
          Grid::index_t index = targetDistGrid().getIndex(indices);
          double energy = targetDistGrid().getPoint(index)[0];
          double volume = targetDistGrid().getPoint(index)[1];
          double value = targetDistGrid().getValue(index);
          if (value>(1-1.e-5)) { // Apply only if this grid point was 1.
            // Apply gaussians around
            std::vector<int> minBin(2), maxBin(2), deltaBin(2); // These cannot be unsigned
            // Only consider contributions less than n*sigma bins apart from the actual distance
            deltaBin[0]=std::floor(5*sigma_[0]/dx[0]);;
            deltaBin[1]=std::floor(5*sigma_[1]/dx[1]);;
            // For energy
            minBin[0]=i - deltaBin[0];
            if (minBin[0] < 0) minBin[0]=0;
            if (minBin[0] > (nbin[0]-1)) minBin[0]=nbin[0]-1;
            maxBin[0]=i +  deltaBin[0];
            if (maxBin[0] > (nbin[0]-1)) maxBin[0]=nbin[0]-1;
            // For volume
            minBin[1]=j - deltaBin[1];
            if (minBin[1] < 0) minBin[1]=0;
            if (minBin[1] > (nbin[1]-1)) minBin[1]=nbin[1]-1;
            maxBin[1]=j +  deltaBin[1];
            if (maxBin[1] > (nbin[1]-1)) maxBin[1]=nbin[1]-1;
            for(unsigned l=minBin[0]; l<maxBin[0]+1; l++) {
              for(unsigned m=minBin[1]; m<maxBin[1]+1; m++) {
                std::vector<unsigned> indices_prime(2);
                indices_prime[0]=l;
                indices_prime[1]=m;
                Grid::index_t index_prime = targetDistGrid().getIndex(indices_prime);
                double energy_prime = targetDistGrid().getPoint(index_prime)[0];
                double volume_prime = targetDistGrid().getPoint(index_prime)[1];
                double value_prime = targetDistGrid().getValue(index_prime);
                // Apply gaussian
                double gaussian_value = GaussianSwitchingFunc(energy_prime,energy,sigma_[0])*GaussianSwitchingFunc(volume_prime,volume,sigma_[1]);
                if (value_prime<gaussian_value) {
                  targetDistGrid().setValue(index_prime,gaussian_value);
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
    } else plumed_merror(getName()+": Number of arguments for this target distribution must be 1 or 2");
  }
}

inline
double TD_Multicanonical::GaussianSwitchingFunc(const double argument, const double center, const double sigma) const {
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
