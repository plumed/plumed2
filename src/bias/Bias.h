/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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

namespace PLMD{
namespace bias{

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
  std::vector<double> outputForces;
protected:
  void resetOutputForces();
  void setOutputForce(int i,double g);
public:
  static void registerKeywords(Keywords&);
  explicit Bias(const ActionOptions&ao);
  void apply();
  unsigned getNumberOfDerivatives();
  void turnOnDerivatives();
};

inline
void Bias::setOutputForce(int i,double f){
  outputForces[i]=f;
}

inline
void Bias::resetOutputForces(){
  for(unsigned i=0;i<outputForces.size();++i) outputForces[i]=0.0;
}

inline
unsigned Bias::getNumberOfDerivatives(){
  return getNumberOfArguments();
}

}
}

#endif

