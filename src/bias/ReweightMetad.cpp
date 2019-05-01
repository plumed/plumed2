/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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
#include "core/ActionRegister.h"
#include "ReweightBase.h"

//+PLUMEDOC REWEIGHTING REWEIGHT_METAD
/*
Calculate the weights configurations should contribute to the histogram in a simulation in which a metadynamics bias acts upon the system.

This command allows you to use the reweighting algorithm discussed in \cite PratyushReweighting when constructing a histogram
of the configurations visited during a metadynamics simulation.

\par Examples

In the following example there is a metadynamics bias acting on the distance between atoms 1 and 2.  Clearly, this
bias will have an effect on the region of phase space that will be sampled when an MD simulation is
run using this variable.  Consequently, when the histogram as a function of the angle, \f$a\f$, is accumulated,
we use reweighting into order to discount the effect of the bias from our final histogram.  We do not use
\ref REWEIGHT_BIAS here, however, as the bias changes with time.  We thus use the reweighting algorithm for
metadynamics instead.  Notice also that we have to specify how often we would like to calculate the c(t) reweighting
factor and the grid over which we calculate c(t) in the input to the METAD command.

\plumedfile
a: ANGLE ATOMS=1,2,3
x: DISTANCE ATOMS=1,2
METAD ARG=x PACE=100 SIGMA=0.1 HEIGHT=1.5 BIASFACTOR=5 GRID_MIN=0 GRID_MAX=10 GRID_BIN=100 REWEIGHTING_NGRID=100 REWEIGHTING_NHILLS=50

bias: REWEIGHT_METAD TEMP=300

HISTOGRAM ...
  ARG=a
  GRID_MIN=0.0
  GRID_MAX=pi
  GRID_BIN=100
  BANDWIDTH=0.1
  LOGWEIGHTS=bias
  LABEL=hB
... HISTOGRAM

DUMPGRID GRID=hB FILE=histoB STRIDE=1 FMT=%8.4f
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace bias {

class ReweightMetad : public ReweightBase {
public:
  static void registerKeywords(Keywords&);
  explicit ReweightMetad(const ActionOptions&ao);
  double getLogWeight();
};

PLUMED_REGISTER_ACTION(ReweightMetad,"REWEIGHT_METAD")

void ReweightMetad::registerKeywords(Keywords& keys ) {
  ReweightBase::registerKeywords( keys ); keys.remove("ARG");
  keys.add("compulsory","ARG","*.rbias","the biases that must be taken into account when reweighting");
}

ReweightMetad::ReweightMetad(const ActionOptions&ao):
  Action(ao),
  ReweightBase(ao)
{
}

double ReweightMetad::getLogWeight() {
  // Retrieve the bias
  double bias=0.0; for(unsigned i=0; i<getNumberOfArguments(); ++i) bias+=getArgument(i);
  return bias / simtemp;
}

}
}
