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

#include "Optimizer.h"
#include "CoeffsVector.h"

#include "core/ActionRegister.h"
#include "core/ActionSet.h"


namespace PLMD {
namespace ves {

//+PLUMEDOC VES_OPTIMIZER OPT_ADAM
/*
Adaptive moment estimation (ADAM) optimizer.

!!! attention ""

    __This optimizer is still experimental. The syntax might change. Restarting does not work. We recommend to use the averaged stochastic gradient decent optimizer ([OPT_AVERAGED_SGD](OPT_AVERAGED_SGD.md)) for now__.

This optimizer updates the coefficients of a VES bias using the ADAM algorithm discussed in the paper cited below.
Running estimates of the first moment $m$ and of the second moment $v$ of the gradient $g$ are updated at every iteration $n$ as

$$
m_n = \beta_1 m_{n-1} + (1-\beta_1) g_n \qquad v_n = \beta_2 v_{n-1} + (1-\beta_2) g_n^2
$$

and each coefficient $c$ is then updated with

$$
c_n = c_{n-1} - \mu_n \frac{m_n}{\sqrt{v_n}+\epsilon} \qquad \mu_n = \mu \frac{\sqrt{1-\beta_2^n}}{1-\beta_1^n}
$$

where $\mu$ is given by STEPSIZE and $\beta_1$, $\beta_2$ and $\epsilon$ are given by BETA_1, BETA_2 and EPSILON.
The AMSGRAD flag uses the maximum of all the past values of $v$ in place of $v_n$, while ADAMW_WEIGHT_DECAY enables the AdamW variant in which the coefficients are also decayed at each iteration.

## Examples

In the following example the coefficients of a VES bias acting on a torsion are optimized with ADAM every 1000 steps.

```plumed
phi:   TORSION ATOMS=5,7,9,15

bf1: BF_FOURIER ORDER=5 MINIMUM=-pi MAXIMUM=pi

VES_LINEAR_EXPANSION ...
 ARG=phi
 BASIS_FUNCTIONS=bf1
 LABEL=ves1
 TEMP=300.0
 GRID_BINS=100
... VES_LINEAR_EXPANSION

OPT_ADAM ...
  BIAS=ves1
  STRIDE=1000
  LABEL=o1
  STEPSIZE=0.001
  COEFFS_FILE=coefficients.data
  COEFFS_OUTPUT=10
  FES_OUTPUT=100
  BIAS_OUTPUT=100
... OPT_ADAM
```

*/
//+ENDPLUMEDOC

class Opt_Adam: public Optimizer {
private:
  // number of updates of each coefficient set (for the bias correction of the moments)
  std::vector<unsigned int> time_;
  double beta_1_;
  double beta_2_;
  double epsilon_;
  double weight_decay_;
  bool amsgrad_;
  bool adamw_;
  // 1st gradient moment uses the "AuxCoeffs", so only 2nd moment needs new CoeffVectors
  std::vector<std::unique_ptr<CoeffsVector>> var_coeffs_pntrs_;
  // used only for AMSGrad variant
  std::vector<std::unique_ptr<CoeffsVector>> varmax_coeffs_pntrs_;
protected:
  CoeffsVector& VarCoeffs(const unsigned int coeffs_id = 0) const;
  CoeffsVector& VarmaxCoeffs(const unsigned int coeffs_id = 0) const;
public:
  static void registerKeywords(Keywords&);
  explicit Opt_Adam(const ActionOptions&);
  void coeffsUpdate(const unsigned int c_id = 0) override;
};

inline
CoeffsVector& Opt_Adam::VarCoeffs(const unsigned int coeffs_id) const {
  return *var_coeffs_pntrs_[coeffs_id];
}

inline
CoeffsVector& Opt_Adam::VarmaxCoeffs(const unsigned int coeffs_id) const {
  return *varmax_coeffs_pntrs_[coeffs_id];
}


PLUMED_REGISTER_ACTION(Opt_Adam,"OPT_ADAM")


void Opt_Adam::registerKeywords(Keywords& keys) {
  Optimizer::registerKeywords(keys);
  Optimizer::useFixedStepSizeKeywords(keys);
  Optimizer::useMultipleWalkersKeywords(keys);
  Optimizer::useMaskKeywords(keys);
  Optimizer::useDynamicTargetDistributionKeywords(keys);
  keys.add("optional","BETA_1","Parameter for the first moment estimate. Defaults to 0.9");
  keys.add("optional","BETA_2","Parameter for the second moment estimate. Defaults to 0.999");
  keys.add("optional","EPSILON","Small parameter to avoid division by zero. Defaults to 1e-8");
  keys.add("optional","ADAMW_WEIGHT_DECAY","Weight decay parameter for the AdamW variant. Defaults to 0");
  keys.addFlag("AMSGRAD", false, "Use the AMSGrad variant");
  keys.addDOI("10.48550/arXiv.1412.6980");
}


Opt_Adam::Opt_Adam(const ActionOptions&ao):
  PLUMED_VES_OPTIMIZER_INIT(ao),
  time_(numberOfCoeffsSets(),0),
  beta_1_(0.9),
  beta_2_(0.999),
  epsilon_(0.00000001),
  weight_decay_(0.0),
  amsgrad_(false),
  adamw_(false),
  var_coeffs_pntrs_(0) {
  // add citation and print it to log
  log << "  Adam type stochastic gradient decent\n";
  parseFlag("AMSGRAD",amsgrad_);
  if (amsgrad_) {
    log << "  Using the AMSGrad variant of the Adam algorithm, see and cite\n";
  }

  double tmp_weight_decay = 0.0;
  parse("ADAMW_WEIGHT_DECAY",tmp_weight_decay);
  if (tmp_weight_decay != 0.0) {
    adamw_ = true;
    log << "  Using the AdamW variant (Adam with weight decay), see and cite\n";
    weight_decay_ = tmp_weight_decay;
    log << "    weight decay parameter: " << tmp_weight_decay << "\n";
  }

  log << "  Adam parameters:\n";
  parse("BETA_1",beta_1_);
  plumed_massert(beta_1_ > 0 && beta_1_ <= 1, "BETA_1 must be between 0 and 1");
  log << "    beta_1: " << beta_1_ << "\n";

  parse("BETA_2",beta_2_);
  plumed_massert(beta_2_ > 0 && beta_2_ <= 1, "BETA_2 must be between 0 and 1");
  log << "    beta_2: " << beta_2_ << "\n";

  parse("EPSILON",epsilon_);
  plumed_massert(epsilon_ > 0 && epsilon_ <= 1, "EPSILON must be between 0 and 1");
  log << "    epsilon: " << epsilon_ << "\n";


  // set up the coeff vector for the 2nd moment of the gradient (variance)
  for (unsigned i = 0; i < numberOfCoeffsSets(); ++i) {
    var_coeffs_pntrs_.emplace_back(std::unique_ptr<CoeffsVector>(new CoeffsVector(Coeffs(i))));
    VarCoeffs(i).replaceLabelString("coeffs","grad_var");
    VarCoeffs(i).setAllValuesToZero(); // can Coeffs(i) even be non-zero at this point?

    // add second set of coefficients to store the maximum values of the 2nd moment
    if (amsgrad_) {
      varmax_coeffs_pntrs_.emplace_back(std::unique_ptr<CoeffsVector>(new CoeffsVector(VarCoeffs(i))));
      VarmaxCoeffs(i).replaceLabelString("coeffs","grad_varmax");
    }

    // also rename the Coeffs used for the mean of the gradient
    AuxCoeffs(i).replaceLabelString("coeffs","grad_mean");
  }

  checkRead();
}


void Opt_Adam::coeffsUpdate(const unsigned int c_id) {
  time_[c_id]++;
  // AuxCoeffs is used for first moment (mean)
  AuxCoeffs(c_id) *= beta_1_;
  AuxCoeffs(c_id) += (1 - beta_1_ ) * Gradient(c_id) * CoeffsMask(c_id);
  VarCoeffs(c_id) *= beta_2_;
  VarCoeffs(c_id) += (1 - beta_2_ ) * Gradient(c_id) * Gradient(c_id) * CoeffsMask(c_id);

  if (amsgrad_) {
    for (size_t i = 0; i< VarCoeffs(c_id).getSize(); ++i) {
      if (VarCoeffs(c_id).getValue(i) > VarmaxCoeffs(c_id).getValue(i)) {
        VarmaxCoeffs(c_id)[i] = VarCoeffs(c_id).getValue(i);
      }
    }
  }

  // store sqrt of VarCoeffs in vector, easier than writing a CoeffsVector::sqrt() function
  // also directly add epsilon and invert to multiply with the Coeffs in last step
  std::vector<double> var_coeffs_sqrt;
  if (!amsgrad_) {
    for (size_t i = 0; i< VarCoeffs(c_id).getSize(); ++i) {
      var_coeffs_sqrt.push_back(1 / (sqrt(VarCoeffs(c_id).getValue(i)) + epsilon_));
    }
  } else { // use VarmaxCoffs instead of VarCoeffs
    for (size_t i = 0; i< VarmaxCoeffs(c_id).getSize(); ++i) {
      var_coeffs_sqrt.push_back(1 / (sqrt(VarmaxCoeffs(c_id).getValue(i)) + epsilon_));
    }
  }

  // bias correction
  double scalefactor = StepSize(c_id) * sqrt(1 - pow(beta_2_, time_[c_id])) / (1 - pow(beta_1_, time_[c_id]));

  if (adamw_) { // check is not necessary but probably faster than always multiplying by 1
    // only the coefficients that are optimized decay (multiplying by the mask set the others to zero)
    for (size_t i = 0; i< Coeffs(c_id).getSize(); ++i) {
      Coeffs(c_id)[i] *= 1 - weight_decay_ * CoeffsMask(c_id).getValue(i);
    }
  }

  // coeff update
  Coeffs(c_id) -= scalefactor * AuxCoeffs(c_id) * var_coeffs_sqrt * CoeffsMask(c_id);
}


}
}
