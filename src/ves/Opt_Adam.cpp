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

#include "Optimizer.h"
#include "CoeffsVector.h"

#include "core/ActionRegister.h"
#include "core/ActionSet.h"
#include <memory>


namespace PLMD {
namespace ves {

//+PLUMEDOC VES_OPTIMIZER OPT_ADAM
/*
Adaptive moment estimation (adam) optimizer.


\par Examples

*/
//+ENDPLUMEDOC

class Opt_Adam: public Optimizer {
private:
  unsigned int time_;
  double beta_1_;
  double beta_2_;
  double epsilon_;
  // 1st moment uses the "AuxCoeffs", so only 2nd moment needs new coeff vectors
  std::vector<std::unique_ptr<CoeffsVector>> var_coeffs_pntrs_;
protected:
  CoeffsVector& VarCoeffs(const unsigned int coeffs_id = 0) const;
public:
  static void registerKeywords(Keywords&);
  explicit Opt_Adam(const ActionOptions&);
  ~Opt_Adam() override =default;
  void coeffsUpdate(const unsigned int c_id = 0) override;
};

inline
CoeffsVector& Opt_Adam::VarCoeffs(const unsigned int coeffs_id) const {return *var_coeffs_pntrs_[coeffs_id];}


PLUMED_REGISTER_ACTION(Opt_Adam,"OPT_ADAM")


void Opt_Adam::registerKeywords(Keywords& keys) {
  Optimizer::registerKeywords(keys);
  Optimizer::useFixedStepSizeKeywords(keys);
  Optimizer::useMultipleWalkersKeywords(keys);
  Optimizer::useMaskKeywords(keys);
  Optimizer::useDynamicTargetDistributionKeywords(keys);
  keys.add("optional","BETA_1","parameter for the first moment estimate. Defaults to 0.9");
  keys.add("optional","BETA_2","parameter for the second moment estimate. Defaults to 0.999");
  keys.add("optional","EPSILON","-parameter for the second moment estimate. Defaults to 1e-8");
}


Opt_Adam::Opt_Adam(const ActionOptions&ao):
  PLUMED_VES_OPTIMIZER_INIT(ao),
  time_(0),
  beta_1_(0.9),
  beta_2_(0.999),
  epsilon_(0.00000001),
  var_coeffs_pntrs_(0)
{
  // add citation and print it to log
  log.printf("  Adam type stochastic gradient decent\n");

  parse("BETA_1",beta_1_);
  parse("BETA_2",beta_1_);
  parse("EPSILON",epsilon_);

  // set up the coeff vector for the 2nd moment (variance)
  for (unsigned i = 0; i < numberOfCoeffsSets(); ++i)
  {
      auto var_coeffs_tmp = new CoeffsVector(Coeffs(i));
      std::string var_label = Coeffs(i).getLabel();
      if(var_label.find("coeffs")!=std::string::npos) {
        var_label.replace(var_label.find("coeffs"), std::string("coeffs").length(), "var_coeffs");
      }
      else {
        var_label += "_var";
      }
      var_coeffs_tmp->setLabels(var_label);
      var_coeffs_pntrs_.push_back(std::unique_ptr<CoeffsVector>(var_coeffs_tmp));
      VarCoeffs(i).setValues(0.0);
  }

  checkRead();
}


void Opt_Adam::coeffsUpdate(const unsigned int c_id) {
  time_++;
  // AuxCoeffs is used for first moment (mean)
  AuxCoeffs(c_id) *= beta_1_;
  AuxCoeffs(c_id) += (1 - beta_1_ ) * Gradient(c_id) * CoeffsMask(c_id);
  VarCoeffs(c_id) *= beta_2_;
  VarCoeffs(c_id) += (1 - beta_2_ ) * Gradient(c_id) * Gradient(c_id) * CoeffsMask(c_id);

  // store sqrt of VarCoeffs in vector, easier than writing a CoeffsVector::sqrt() function
  // also directly add epsilon and invert to multiply with the Coeffs in last step
  std::vector<double> var_coeffs_sqrt;
  for (size_t i = 0; i< VarCoeffs(c_id).getSize(); ++i) {
    var_coeffs_sqrt.push_back(1 / (sqrt(VarCoeffs(c_id).getValue(i)) + epsilon));
  }

  // bias correction
  double scalefactor = StepSize(c_id) * sqrt(1 - pow(beta_2_, time_)) / (1 - pow(beta_1_, time_));

  // coeff update
  Coeffs(c_id) -= scalefactor * AuxCoeffs(c_id) * var_coeffs_sqrt;
}


}
}
