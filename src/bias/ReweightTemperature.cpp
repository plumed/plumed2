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
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/Atoms.h"
#include "ReweightBase.h"

//+PLUMEDOC REWEIGHTING REWEIGHT_TEMP
/*
Calculate weights for ensemble averages allow for the computing of ensemble averages at temperatures lower/higher than that used in your original simulation.

We can use our knowledge of the Boltzmann distribution in the canonical ensemble to reweight the data
contained in trajectories.  Using this procedure we can take trajectory at temperature \f$T_1\f$ and use it to
extract probabilities at a different temperature, \f$T_2\f$, using:

\f[
P(s',t) = \frac{ \sum_{t'}^t \delta( s(x) - s' ) \exp\left( +( \left[\frac{1}{T_1} - \frac{1}{T_2}\right] \frac{U(x,t')}{k_B} \right) }{ \sum_{t'}^t \exp\left( +\left[\frac{1}{T_1} - \frac{1}{T_2}\right] \frac{U(x,t')}{k_B} \right) }
\f]

The weights calculated by this action are equal to \f$\exp\left( +( \left[\frac{1}{T_1} - \frac{1}{T_2}\right] \frac{U(x,t')}{k_B} \right)\f$ and take
the effect the bias has on the system into account.  These weights can be used in any action
that computes ensemble averages.  For example this action can be used in tandem with \ref HISTOGRAM or \ref AVERAGE.

\par Examples

The following input can be used to post process a molecular dynamics trajectory calculated at a temperature of 500 K.
The \ref HISTOGRAM as a function of the distance between atoms 1 and 2 that would have been obtained if the simulation
had been run at the lower temperature of 300 K is estimated using the data from the higher temperature trajectory and output
to a file.

\plumedfile
x: DISTANCE ATOMS=1,2
aa: REWEIGHT_TEMP TEMP=500 REWEIGHT_TEMP=300
hB: HISTOGRAM ARG=x GRID_MIN=0.0 GRID_MAX=3.0 GRID_BIN=100 BANDWIDTH=0.1 LOGWEIGHTS=aa
DUMPGRID GRID=hB FILE=histoB
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace bias {

class ReweightTemperature : public ReweightBase {
private:
///
  double rtemp;
/// The biases we are using in reweighting and the args we store them separately
  std::vector<Value*> biases;
public:
  static void registerKeywords(Keywords&);
  explicit ReweightTemperature(const ActionOptions&ao);
  void prepare();
  double getLogWeight();
};

PLUMED_REGISTER_ACTION(ReweightTemperature,"REWEIGHT_TEMP")

void ReweightTemperature::registerKeywords(Keywords& keys ) {
  ReweightBase::registerKeywords( keys );
  keys.add("compulsory","REWEIGHT_TEMP","reweight data from a trajectory at one temperature and output the probability "
           "distribution at a second temperature. This is not possible during post processing.");
}

ReweightTemperature::ReweightTemperature(const ActionOptions&ao):
  Action(ao),
  ReweightBase(ao)
{
  parse("REWEIGHT_TEMP",rtemp);
  log.printf("  reweighting simulation to probabilities at temperature %f\n",rtemp);
  rtemp*=plumed.getAtoms().getKBoltzmann();

  std::vector<ActionWithValue*> all=plumed.getActionSet().select<ActionWithValue*>();
  if( all.empty() ) error("your input file is not telling plumed to calculate anything");
  log.printf("  using the following biases in reweighting ");
  for(unsigned j=0; j<all.size(); j++) {
    std::string flab; flab=all[j]->getLabel() + ".bias";
    if( all[j]->exists(flab) ) {
      biases.push_back( all[j]->copyOutput(flab) );
      log.printf(" %s", flab.c_str());
    }
  }
  log.printf("\n");
}

void ReweightTemperature::prepare() {
  plumed.getAtoms().setCollectEnergy(true);
}

double ReweightTemperature::getLogWeight() {
  // Retrieve the bias
  double bias=0.0; for(unsigned i=0; i<biases.size(); ++i) bias+=biases[i]->get();
  double energy=plumed.getAtoms().getEnergy()+bias;
  return -( (1.0/rtemp) - (1.0/simtemp) )*(energy+bias);
}

}
}
