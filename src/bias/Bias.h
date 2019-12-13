/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2019 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_bias_Bias_h
#define __PLUMED_bias_Bias_h

#include "core/ActionPilot.h"
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"

#define PLUMED_BIAS_INIT(ao) Action(ao),Bias(ao)

namespace PLMD {
namespace bias {

/**
\ingroup INHERIT
This is the abstract base class to use for implementing new simulation biases, within it there is
information as to how to go about implementing a new bias.
*/

class Bias :
  public ActionPilot,
  public ActionWithValue,
  public ActionWithArguments
{
/// the vector of the forces
  std::vector<double> outputForces;
/// the pointer to the bias component
  Value *valueBias;
protected:
/// set the force from the bias on argument i, this automatically set the partial derivative of the bias with respect to i to -f
  void setOutputForce(int i,double f);
/// set the value of the bias
  void setBias(double bias);
public:
  static void registerKeywords(Keywords&);
  explicit Bias(const ActionOptions&ao);
  void apply() override;
  unsigned getNumberOfDerivatives() override;
};

inline
void Bias::setOutputForce(int i,double f) {
  outputForces[i]=f;
  valueBias->addDerivative(i,-f);
}

inline
void Bias::setBias(double bias) {
  valueBias->set(bias);
}

inline
unsigned Bias::getNumberOfDerivatives() {
  return getNumberOfArguments();
}

}
}

#endif

