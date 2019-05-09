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

//+PLUMEDOC REWEIGHTING REWEIGHT_TEMP_PRESS
/*
Calculate weights for ensemble averages at temperatures and pressures different than that used in your original simulation.

We can use our knowledge of the probability distribution in the isothermal-isobaric ensemble to reweight the data
contained in trajectories and obtain ensemble averages at different temperatures and pressures.

Consider the ensemble average of an observable \f$O(\mathbf{R},\mathcal{V})\f$ that depends on the atomic coordinates \f$\mathbf{R}\f$ and the volume \f$\mathcal{V}\f$.
This observable is in practice any collective variable (CV) calculated by Plumed.
The ensemble average of this observable in the isothermal-isobaric ensemble at (inverse) temperature \f$\beta'\f$ and pressure \f$P'\f$ can be calculated from a simulation performed at temperature \f$\beta\f$ and pressure \f$P\f$ using:
\f[
\langle O(\mathbf{R},\mathcal{V}) \rangle_{\beta',P'} = \frac{\langle O(\mathbf{R},\mathcal{V}) e^{(\beta-\beta')E(\mathbf{R}) + (\beta P - \beta' P') \mathcal{V}}  \rangle_{\beta,P}}
                                                     {\langle e^{(\beta-\beta')E(\mathbf{R}) + (\beta P - \beta' P') \mathcal{V}}  \rangle_{\beta,P}}
\f]
where \f$\langle \cdot \rangle_{\beta,P}\f$ is a mean value in the isothermal-isobaric ensemble at temperature \f$\beta\f$ and pressure \f$P\f$, and \f$ E(\mathbf{R}) \f$ is the potential energy of the system.

The weights calculated by this action are equal to \f$ e^{(\beta-\beta')E(\mathbf{R}) + (\beta P - \beta' P') \mathcal{V}} \f$.
These weights can be used in any action that computes ensemble averages. 
For example this action can be used in tandem with \ref HISTOGRAM or \ref AVERAGE.

The above equation is often impractical since the overlap between the distributions of energy and volume at different temperatures and pressures is only significant for neighboring temperatures and pressures.
For this reason an unbiased simulation is of little use to reweight at different temperatures and pressures.
A successful approach has been altering the probability in the isothermal-isobaric simulation in order to increase this overlap \cite wanglandau.
This is done through a bias potential \f$ V(\mathbf{s}) \f$ where \f$ \mathbf{s} \f$ is a set of CVs, that often is the energy (and possibly the volume).
In order to calculate ensemble averages, also the effect of this bias must be taken into account.
The ensemble average of the observable in the isothermal-isobaric ensemble at (inverse) temperature \f$\beta'\f$ and pressure \f$P'\f$ can be calculated from a biased simulation performed at temperature \f$\beta\f$ and pressure \f$P\f$ with bias \f$ V(\mathbf{s}) \f$ using:
\f[
\langle O(\mathbf{R},\mathcal{V}) \rangle_{\beta',P'} = \frac{\langle O(\mathbf{R},\mathcal{V}) e^{(\beta-\beta')E(\mathbf{R}) + (\beta P - \beta' P') \mathcal{V}} e^{\beta V(\mathbf{s})}  \rangle_{\beta,P,V}}
                                                     {\langle e^{(\beta-\beta')E(\mathbf{R}) + (\beta P - \beta' P') \mathcal{V}} e^{\beta V(\mathbf{s})}  \rangle_{\beta,P,V}}
\f]
where \f$\langle \cdot \rangle_{\beta,P,V}\f$ is a mean value in the biased ensemble at temperature \f$\beta\f$ and pressure \f$P\f$ with static bias \f$ V(\mathbf{s}) \f$.
Therefore in order to reweight the trajectory at different temperatures and pressures one must use the weights calculated by this action \f$ e^{(\beta-\beta')E(\mathbf{R}) + (\beta P - \beta' P') \mathcal{V}} \f$ together with the weights of \ref REWEIGHT_BIAS (see example below).

The bias potential \f$ V(\mathbf{s}) \f$ can be constructed with \ref METAD using \ref ENERGY as a CV \cite mich+04prl.
More specialized tools are available, for instance using bespoke target distributions such as \ref TD_MULTICANONICAL and \ref TD_MULTITHERMAL_MULTIBARIC \cite Piaggi-PRL-2019 \cite Piaggi-arXiv-2019 within \ref VES.
In the latter algorithms the interval of temperatures and pressures in which the trajectory can be reweighted is chosen explicitly.

\par Examples

The following input can be used to postprocess a molecular dynamics trajectory of a system of 1000 particles run at 500 K and 1 bar.
We read from a file COLVAR the potential energy, the volume, and the value of the bias potential and calculate the ensemble average of the (particle) density at 300 K and 300 MPa.

\plumedfile
energy: READ FILE=COLVAR VALUES=energy  IGNORE_TIME
volume: READ FILE=COLVAR VALUES=volume  IGNORE_TIME
mybias: READ FILE=COLVAR VALUES=mybias.bias  IGNORE_TIME

# Shift energy and volume
rvol: COMBINE ARG=vol PARAMETERS=7.8 PERIODIC=NO
renergy: COMBINE ARG=energy PARAMETERS=-13250 PERIODIC=NO

# Weights
bias_weights: REWEIGHT_BIAS TEMP=500 ARG=mybias.bias
temp_press_weights: REWEIGHT_TEMP_PRESS TEMP=500 REWEIGHT_TEMP=300 PRESSURE=0.06022140857 REWEIGHT_PRESSURE=180.66422571 ARG=renergy,rvolume

# Ensemble average of the volume at 300 K and 300 MPa
avg_vol: AVERAGE ARG=volume LOGWEIGHTS=bias_weights,temp_press_weights
# Ensemble average of the density at 300 K and 300 MPa
avg_density: CUSTOM ARG=avg_vol FUNC=1000/x PERIODIC=NO

PRINT ARG=avg_density FILE=COLVAR_REWEIGHT STRIDE=1
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace bias {

class ReweightTemperaturePressure : public ReweightBase {
private:
///
  double rpress_, press_, rtemp_;
public:
  static void registerKeywords(Keywords&);
  explicit ReweightTemperaturePressure(const ActionOptions&ao);
  double getLogWeight();
};

PLUMED_REGISTER_ACTION(ReweightTemperaturePressure,"REWEIGHT_TEMP_PRESS")

void ReweightTemperaturePressure::registerKeywords(Keywords& keys ) {
  ReweightBase::registerKeywords( keys );
  keys.add("compulsory","REWEIGHT_PRESSURE","Reweighting pressure");
  keys.add("compulsory","PRESSURE","The system pressure");
  keys.add("compulsory","REWEIGHT_TEMP","Reweighting temperature");
  keys.remove("ARG"); keys.add("compulsory","ARG","This action takes two arguments: the energy and the volume, in that order");
}

ReweightTemperaturePressure::ReweightTemperaturePressure(const ActionOptions&ao):
  Action(ao),
  ReweightBase(ao)
{
  parse("REWEIGHT_PRESSURE",rpress_);
  parse("PRESSURE",press_);
  parse("REWEIGHT_TEMP",rtemp_);
  log.printf("  reweighting simulation to probabilities at temperature %f and pressure %f \n",rtemp_,rpress_);
  rtemp_*=plumed.getAtoms().getKBoltzmann();
  if(getNumberOfArguments()!=2) error(" Number of arguments should be 2. The energy and the volume, in that order.");
}

double ReweightTemperaturePressure::getLogWeight() {
  double energy=getArgument(0);
  double volume=getArgument(1);
  return ((1.0/simtemp)- (1.0/rtemp_) )*energy + ((1.0/simtemp)*press_ - (1.0/rtemp_)*rpress_ )*volume;
}

}
}
