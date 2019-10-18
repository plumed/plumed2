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
public:
  static void registerKeywords(Keywords&);
  explicit Opt_Adam(const ActionOptions&);
  void coeffsUpdate(const unsigned int c_id = 0);
};


PLUMED_REGISTER_ACTION(Opt_Adam,"OPT_ADAM")


void Opt_Adam::registerKeywords(Keywords& keys) {
  Optimizer::registerKeywords(keys);
  Optimizer::useFixedStepSizeKeywords(keys);
  Optimizer::useMultipleWalkersKeywords(keys);
  Optimizer::useMaskKeywords(keys);
  keys.add("optional","BETA_1","parameter for the first moment estimate. Defaults to 0.9");
  keys.add("optional","BETA_2","parameter for the second moment estimate. Defaults to 0.999");
  keys.add("optional","EPSILON","-parameter for the second moment estimate. Defaults to 1e-8");
}


Opt_Adam::Opt_Adam(const ActionOptions&ao):
  PLUMED_VES_OPTIMIZER_INIT(ao),
  time_(0),
  beta_1_(0.9),
  beta_2_(0.999),
  epsilon_(0.00000001)
{
  // add citation and print it to log
  log.printf("  Adam type stochastic gradient decent\n");

  parse("BETA_1",beta_1_);
  //if (beta_1_ != 0.9) {addKeywordToList("BETA_1",beta_1_);}
  parse("BETA_2",beta_1_);
  parse("EPSILON",epsilon_);

  checkRead();
}


void Opt_Adam::coeffsUpdate(const unsigned int c_id) {
  Coeffs(c_id) += - StepSize(c_id) * CoeffsMask(c_id) * Gradient(c_id);
  //
  double aver_decay = 1.0 / ( getIterationCounterDbl() + 1.0 );
  AuxCoeffs(c_id) += aver_decay * ( Coeffs(c_id)-AuxCoeffs(c_id) );

}


}
}
