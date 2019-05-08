/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The plumed team
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
#include "core/Atoms.h"
#include "ReweightBase.h"

//+PLUMEDOC REWEIGHTING REWEIGHT_TEMP
/*
Calculate weights for ensemble averages allow for the computing of ensemble averages at temperatures lower/higher than that used in your original simulation.

We can use our knowledge of the Boltzmann distribution in the cannonical ensemble to reweight the data
contained in trajectories.  Using this procedure we can take trajectory at temperature \f$T_1\f$ and use it to
extract probabilities at a different temperature, \f$T_2\f$, using:

\f[
P(s',t) = \frac{ \sum_{t'}^t \delta( s(x) - s' ) \exp\left( +( \left[\frac{1}{T_1} - \frac{1}{T_2}\right] \frac{U(x,t')}{k_B} \right) }{ \sum_{t'}^t \exp\left( +\left[\frac{1}{T_1} - \frac{1}{T_2}\right] \frac{U(x,t')}{k_B} \right) }
\f]

The weights calculated by this action are equal to \f$\exp\left( +( \left[\frac{1}{T_1} - \frac{1}{T_2}\right] \frac{U(x,t')}{k_B} \right)\f$ and take
the effect the bias has on the system into account.  These weights can be used in any action
that computes ensemble averages.  For example this action can be used in tandem with \ref HISTOGRAM or \ref AVERAGE.

\par Examples

The following input can be used to postprocess a molecular dynamics trajectory calculated at a temperature of 500 K.
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

class ReweightTemperaturePressure : public ReweightBase {
private:
///
  double rpress_, press_, rtemp_;
/// The biases we are using in reweighting and the args we store them separately
  std::vector<Value*> biases;
public:
  static void registerKeywords(Keywords&);
  explicit ReweightTemperaturePressure(const ActionOptions&ao);
  double getLogWeight();
};

PLUMED_REGISTER_ACTION(ReweightTemperaturePressure,"REWEIGHT_TEMPERATURE_PRESSURE")

void ReweightTemperaturePressure::registerKeywords(Keywords& keys ) {
  ReweightBase::registerKeywords( keys );
  keys.add("compulsory","REWEIGHT_PRESSURE","reweight data from a trajectory at one pressure and output the probability "
           "distribution at a second temperature. This is not possible during postprocessing.");
  keys.add("compulsory","PRESSURE","Pressure of the trajectory");
  //keys.add("compulsory","DOF","Degrees of freedom of the kinetic energy");
  keys.add("compulsory","REWEIGHT_TEMP","reweight data from a trajectory at one temperature and output the probability ");
  keys.remove("ARG"); keys.add("compulsory","ARG","energy,volume");
}

ReweightTemperaturePressure::ReweightTemperaturePressure(const ActionOptions&ao):
  Action(ao),
  ReweightBase(ao)
{
  parse("REWEIGHT_PRESSURE",rpress_);
  parse("PRESSURE",press_);
  parse("REWEIGHT_TEMP",rtemp_);
  //parse("DOF",dof_);
  log.printf("  reweighting simulation to probabilities at temperature %f\n",rtemp_);
  rtemp_*=plumed.getAtoms().getKBoltzmann();
  //rtemp*=plumed.getAtoms().getKBoltzmann();
  if(getNumberOfArguments()!=2) error(" Number of arguments should be 2. The energy and the volume, in that order.");
}

double ReweightTemperaturePressure::getLogWeight() {
  double energy=getArgument(0);
  double volume=getArgument(1);
  //double deltaPressure=dof_*(rtemp_-simtemp)/(3*volume);
  //double correctedPressure=press_+deltaPressure;
  //log.printf("Press %f deltaPressure %f correctedPressure %f \n",press_/0.0602,deltaPressure/0.0602,correctedPressure/0.0602);
  return ((1.0/simtemp)- (1.0/rtemp_) )*energy + ((1.0/simtemp)*press_ - (1.0/rtemp_)*rpress_ )*volume;
}

}
}
