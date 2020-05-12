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

#include "core/ActionRegister.h"


namespace PLMD {
namespace ves {

//+PLUMEDOC VES_OPTIMIZER OPT_DUMMY
/*
Dummy optimizer for debugging.

This is dummy optimizer that can be used for debugging. It will not update the
coefficients but can be used to monitor the gradient and Hessian for a given
VES bias.

\par Examples

In the following input we use the OPT_DUMMY to monitor the gradient and
Hessian for a given VES bias every 1 iteration.
\plumedfile
phi:   TORSION ATOMS=5,7,9,15

bf1: BF_FOURIER ORDER=5 MINIMUM=-pi MAXIMUM=pi

VES_LINEAR_EXPANSION ...
 ARG=phi
 BASIS_FUNCTIONS=bf1
 LABEL=ves1
 TEMP=300.0
 GRID_BINS=100
... VES_LINEAR_EXPANSION

OPT_DUMMY ...
  BIAS=ves1
  STRIDE=1000
  LABEL=o1
  MONITOR_HESSIAN
  GRADIENT_FILE=gradient.data
  GRADIENT_OUTPUT=1
  GRADIENT_FMT=%12.6f
  HESSIAN_FILE=hessian.data
  HESSIAN_OUTPUT=1
  HESSIAN_FMT=%12.6f
... OPT_DUMMY
\endplumedfile

*/
//+ENDPLUMEDOC



class Opt_Dummy : public Optimizer {

public:
  static void registerKeywords(Keywords&);
  explicit Opt_Dummy(const ActionOptions&);
  void coeffsUpdate(const unsigned int c_id = 0) override;
};


PLUMED_REGISTER_ACTION(Opt_Dummy,"OPT_DUMMY")


void Opt_Dummy::registerKeywords(Keywords& keys) {
  Optimizer::registerKeywords(keys);
  //
  Optimizer::useMultipleWalkersKeywords(keys);
  Optimizer::useHessianKeywords(keys);
  Optimizer::useMonitorAverageGradientKeywords(keys);
  keys.addFlag("MONITOR_HESSIAN",false,"also monitor the Hessian");
}


Opt_Dummy::Opt_Dummy(const ActionOptions&ao):
  PLUMED_VES_OPTIMIZER_INIT(ao)
{
  log.printf("  fake optimizer that does not update coefficients\n");
  log.printf("  can be used to monitor gradient and Hessian for debugging purposes\n");
  bool monitor_hessian = false;
  parseFlag("MONITOR_HESSIAN",monitor_hessian);
  if(monitor_hessian) {
    turnOnHessian();
    log.printf("  the Hessian will also be monitored\n");
  }
  else {
    turnOffHessian();
  }
  turnOffCoeffsOutputFiles();
  checkRead();
}


void Opt_Dummy::coeffsUpdate(const unsigned int c_id) {}


}
}
