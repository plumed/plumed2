/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2023 The plumed team
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
#include "ReweightBase.h"
#include "core/ActionRegister.h"

//+PLUMEDOC REWEIGHTING REWEIGHT_BIAS
/*
Calculate weights for ensemble averages that negate the effect the bias has on the region of phase space explored

If a static or pseudo-static bias $V(x,t')$ is acting on
the system we can remove the bias and get the unbiased probability distribution using:

$$
\langle P(s',t) \rangle = \frac{ \sum_{t'}^t \delta( s(x) - s' ) \exp\left( +\frac{V(x,t')}{k_B T} \right) }{ \sum_t'^t \exp\left( +\frac{V(x,t')}{k_B T} \right) }
$$

The weights calculated by this action are equal to $\exp\left( +\frac{V(x,t')}{k_B T} \right)$ these weights can then be used in any action
that computes ensemble averages.  For example this action can be used in tandem with [HISTOGRAM](HISTOGRAM.md) or [AVERAGE](AVERAGE.md).

## Examples

In the following example there is a fixed restraint on the distance between atoms 1 and 2.  Clearly, this
restraint will have an effect on the region of phase space that will be sampled when an MD simulation is
run using this variable.  Consequently, when the histogram as a function of the distance, $x$, is accumulated,
we use reweighting into order to discount the effect of the bias from our final histogram.

```plumed
x: DISTANCE ATOMS=1,2
RESTRAINT ARG=x SLOPE=1.0 AT=0.0
bias: REWEIGHT_BIAS TEMP=300

hB: HISTOGRAM ...
  ARG=x
  GRID_MIN=0.0
  GRID_MAX=3.0
  GRID_BIN=100
  BANDWIDTH=0.1
  LOGWEIGHTS=bias
...

DUMPGRID ARG=hB FILE=histoB STRIDE=1
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace bias {

class ReweightBias : public ReweightBase {
public:
  static void registerKeywords(Keywords&);
  explicit ReweightBias(const ActionOptions&ao);
  double getLogWeight() override;
};

PLUMED_REGISTER_ACTION(ReweightBias,"REWEIGHT_BIAS")

void ReweightBias::registerKeywords(Keywords& keys ) {
  ReweightBase::registerKeywords( keys );
  keys.addInputKeyword("compulsory","ARG","scalar","*.bias","the biases that must be taken into account when reweighting");
  keys.setValueDescription("scalar","the weight to use for this frame to negate the effect the bias");
}

ReweightBias::ReweightBias(const ActionOptions&ao):
  Action(ao),
  ReweightBase(ao) {
}

double ReweightBias::getLogWeight() {
  // Retrieve the bias
  double bias=0.0;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    bias+=getArgument(i);
  }
  return bias / simtemp;
}

}
}
