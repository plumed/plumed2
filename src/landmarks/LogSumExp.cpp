/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2023 The plumed team
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
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"

//+PLUMEDOC ANALYSIS LOGSUMEXP
/*
This action takes the exponential of a vector of logarithms and divides each element of the vector by the sum of the exponentials.

If you has performed a simulation with a time-independent bias, $V(x)$, you can recover the unbiased distribution from the trajectory
by constructing a histogram in which the $i$th frame of the trajectory is given weight of:

$$
w_i = \frac{e^{\beta V(x_i)}}{\sum_j e^{\beta V(x_j)}}
$$

where $x_i$ is used to indicate the coordinates for frame $i$ and $\beta$ is the inverse temperature in units of energy.  If the bias
is large then it is easy for $e^{\beta V(x_i)}$ to overflow.  This action thus allows you to use the following trick is to compute the $w_i$
values in the above expression:

$$
w_i = e^{\beta(V(x_i) - c} \qquad \textrm{where} \qquad c = \beta V_\textrm{max} + \log\left[ \sum_j e^{\beta V(x_j) - \beta V_\textrm{max}} \right]
$$

In this expression $V_\textrm{max}$ is the maximum value the bias took during the simulation.

The following example shows how you can write a PLUMED input that exploits this trick.

```plumed
# Calculate some CVs for later analysis
t1: TORSION ATOMS=1,2,3,4
t2: TORSION ATOMS=5,6,7,8
t3: TORSION ATOMS=9,10,11,12

# This computes the bias potential
r: RESTRAINT ARG=t1 AT=pi/2 KAPPA=10
# This calculates the instantaneous reweighting weight
# given the bias potential that is acting upon the system
bw: REWEIGHT_BIAS TEMP=300

# This collects our data for later analysis and the reweighting weights
cc: COLLECT_FRAMES ARG=t1,t2,t3 LOGWEIGHTS=bw

# And this determines the final vector of weights that should be used for histograms and so on.
weights: LOGSUMEXP ARG=cc_logweights

# This outputs the time series of reweighting weights to a file
DUMPVECTOR ARG=weights FILE=weights
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace landmarks {

class LogSumExp : public ActionShortcut {
private:
  std::string fixArgumentName( const std::string& argin );
public:
  static void registerKeywords( Keywords& keys );
  explicit LogSumExp( const ActionOptions& ao );
};

PLUMED_REGISTER_ACTION(LogSumExp,"LOGSUMEXP")

void LogSumExp::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","ARG","the vector of logweights that you would like to normalise using the logsumexp trick");
  keys.setValueDescription("vector","the logarithms of the input weights logweights that are computed with the log-sum weights formula");
  keys.needsAction("HIGHEST");
  keys.needsAction("CUSTOM");
  keys.needsAction("SUM");
}


LogSumExp::LogSumExp( const ActionOptions& ao ):
  Action(ao),
  ActionShortcut(ao) {
  // Find the argument name
  std::string argn;
  parse("ARG",argn);
  // Find the maximum weight
  readInputLine( getShortcutLabel() + "_maxlogweight: HIGHEST ARG=" + argn );
  // Calculate the maximum
  readInputLine( getShortcutLabel() + "_shiftw: CUSTOM ARG=" + argn + "," + getShortcutLabel() + "_maxlogweight FUNC=exp(x-y) PERIODIC=NO");
  // compute the sum of all the exponentials
  readInputLine( getShortcutLabel() + "_sumw: SUM ARG=" + getShortcutLabel() + "_shiftw PERIODIC=NO");
  // and the logsum
  readInputLine( getShortcutLabel() + "_logsum: CUSTOM ARG=" + getShortcutLabel() + "_sumw," + getShortcutLabel() + "_maxlogweight FUNC=y+log(x) PERIODIC=NO");
  // And the final weights
  readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + argn + "," +  getShortcutLabel() + "_logsum FUNC=exp(x-y) PERIODIC=NO");
}

}
}
