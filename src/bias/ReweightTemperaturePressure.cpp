/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2019 The plumed team
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
Calculate weights for ensemble averages at temperatures and/or pressures different than those used in your original simulation.

We can use our knowledge of the probability distribution in the canonical (N\f$\mathcal{V}\f$T) or the isothermal-isobaric ensemble (NPT) to reweight the data
contained in trajectories and obtain ensemble averages at different temperatures and/or pressures.

Consider the ensemble average of an observable \f$O(\mathbf{R},\mathcal{V})\f$ that depends on the atomic coordinates \f$\mathbf{R}\f$ and the volume \f$\mathcal{V}\f$.
This observable is in practice any collective variable (CV) calculated by Plumed.
The ensemble average of the observable in an ensemble \f$ \xi' \f$  can be calculated from a simulation performed in an ensemble \f$ \xi \f$ using:
\f[
\langle O(\mathbf{R},\mathcal{V}) \rangle_{\xi'} = \frac{\langle O(\mathbf{R},\mathcal{V}) w(\mathbf{R},\mathcal{V}) \rangle_{\xi}}
                                                     {\langle w(\mathbf{R},\mathcal{V}) \rangle_{\xi}}
\f]
where \f$\langle \cdot \rangle_{\xi}\f$ and  \f$\langle \cdot \rangle_{\xi'}\f$ are mean values in the simulated and targeted ensemble, respectively, \f$ E(\mathbf{R}) \f$ is the potential energy of the system, and \f$ w (\mathbf{R},\mathcal{V}) \f$ are the appropriate weights to take from \f$ \xi \f$ to \f$ \xi' \f$.
This action calculates the weights  \f$ w (\mathbf{R},\mathcal{V}) \f$ and handles 4 different cases:
  1. Change of temperature from T to T' at constant volume. That is to say, from a simulation performed in the N\f$\mathcal{V}\f$T (canonical) ensemble, obtain an ensemble average in the N\f$\mathcal{V}\f$T' ensemble. The weights in this case are  \f$ w(\mathbf{R},\mathcal{V}) = e^{(\beta-\beta')E(\mathbf{R})} \f$ with \f$ \beta \f$ and \f$ \beta' \f$ the inverse temperatures.
  2. Change of temperature from T to T' at constant pressure. That is to say, from a simulation performed in the NPT (isothermal-isobaric) ensemble, obtain an ensemble average in the NPT' ensemble. The weights in this case are \f$ w(\mathbf{R},\mathcal{V}) = e^{(\beta-\beta')(E(\mathbf{R}) + P\mathcal{V}) } \f$.
  3. Change of pressure from P to P' at constant temperature. That is to say, from a simulation performed in the NPT (isothermal-isobaric) ensemble, obtain an ensemble average in the NP'T ensemble. The weights in this case are \f$ w(\mathbf{R},\mathcal{V}) = e^{\beta (P - P') \mathcal{V}} \f$.
  4. Change of temperature and pressure from T,P to T',P'. That is to say, from a simulation performed in the NPT (isothermal-isobaric) ensemble, obtain an ensemble average in the NP'T' ensemble. The weights in this case are \f$ w(\mathbf{R},\mathcal{V}) = e^{(\beta-\beta')E(\mathbf{R}) + (\beta P - \beta' P') \mathcal{V}} \f$.

These weights can be used in any action that computes ensemble averages.
For example this action can be used in tandem with \ref HISTOGRAM or \ref AVERAGE.


The above equation is often impractical since the overlap between the distributions of energy and volume at different temperatures and pressures is only significant for neighboring temperatures and pressures.
For this reason an unbiased simulation is of little use to reweight at different temperatures and/or pressures.
A successful approach has been altering the probability of observing a configuration in order to increase this overlap \cite wanglandau.
This is done through a bias potential \f$ V(\mathbf{s}) \f$ where \f$ \mathbf{s} \f$ is a set of CVs, that often is the energy (and possibly the volume).
In order to calculate ensemble averages, also the effect of this bias must be taken into account.
The ensemble average of the observable in the ensemble \f$ \xi' \f$ can be calculated from a biased simulation performed in the ensemble \f$\xi\f$ with bias \f$ V(\mathbf{s}) \f$ using:
\f[
\langle O(\mathbf{R},\mathcal{V}) \rangle_{\xi'} = \frac{\langle O(\mathbf{R},\mathcal{V})  w (\mathbf{R},\mathcal{V}) e^{\beta V(\mathbf{s})}  \rangle_{\xi,V}}
                                                     {\langle w (\mathbf{R},\mathcal{V})  e^{\beta V(\mathbf{s})}  \rangle_{\xi,V}}
\f]
where \f$\langle \cdot \rangle_{\xi,V}\f$ is a mean value in the biased ensemble with static bias \f$ V(\mathbf{s}) \f$.
Therefore in order to reweight the trajectory at different temperatures and/or pressures one must use the weights calculated by this action \f$ w (\mathbf{R},\mathcal{V}) \f$ together with the weights of \ref REWEIGHT_BIAS (see the examples below).

The bias potential \f$ V(\mathbf{s}) \f$ can be constructed with \ref METAD using \ref ENERGY as a CV \cite mich+04prl.
More specialized tools are available, for instance using bespoke target distributions such as \ref TD_MULTICANONICAL and \ref TD_MULTITHERMAL_MULTIBARIC \cite Piaggi-PRL-2019 \cite Piaggi-arXiv-2019 within \ref VES.
In the latter algorithms the interval of temperatures and pressures in which the trajectory can be reweighted is chosen explicitly.

\par Examples

We consider the 4 cases described above.

The following input can be used to postprocess a molecular dynamics trajectory of a system of 1000 particles run at 500 K and constant volume using a static bias potential.
We read from a file COLVAR the potential energy, a distance, and the value of the bias potential and calculate the ensemble average of the distance at 300 K.

\plumedfile
energy: READ FILE=COLVAR VALUES=energy  IGNORE_TIME
distance: READ FILE=COLVAR VALUES=distance  IGNORE_TIME
mybias: READ FILE=COLVAR VALUES=mybias.bias  IGNORE_TIME

# Shift energy (to avoid numerical issues)
renergy: COMBINE ARG=energy PARAMETERS=-13250 PERIODIC=NO

# Weights
bias_weights: REWEIGHT_BIAS TEMP=500 ARG=mybias.bias
temp_press_weights: REWEIGHT_TEMP_PRESS TEMP=500 REWEIGHT_TEMP=300 ENERGY=renergy

# Ensemble average of the distance at 300 K
avg_dist: AVERAGE ARG=distance LOGWEIGHTS=bias_weights,temp_press_weights

PRINT ARG=avg_dist FILE=COLVAR_REWEIGHT STRIDE=1
\endplumedfile

The next three inputs can be used to postprocess a molecular dynamics trajectory of a system of 1000 particles run at 500 K and 1 bar using a static bias potential.

We read from a file COLVAR the potential energy, the volume, and the value of the bias potential and calculate the ensemble average of the (particle) density at 300 K and 1 bar (the simulation temperature was 500 K).

\plumedfile
energy: READ FILE=COLVAR VALUES=energy  IGNORE_TIME
volume: READ FILE=COLVAR VALUES=volume  IGNORE_TIME
mybias: READ FILE=COLVAR VALUES=mybias.bias  IGNORE_TIME

# Shift energy and volume (to avoid numerical issues)
rvol: COMBINE ARG=vol PARAMETERS=7.8 PERIODIC=NO
renergy: COMBINE ARG=energy PARAMETERS=-13250 PERIODIC=NO

# Weights
bias_weights: REWEIGHT_BIAS TEMP=500 ARG=mybias.bias
temp_press_weights: REWEIGHT_TEMP_PRESS TEMP=500 REWEIGHT_TEMP=300 PRESSURE=0.06022140857 ENERGY=renergy VOLUME=rvolume

# Ensemble average of the volume at 300 K
avg_vol: AVERAGE ARG=volume LOGWEIGHTS=bias_weights,temp_press_weights
# Ensemble average of the density at 300 K
avg_density: CUSTOM ARG=avg_vol FUNC=1000/x PERIODIC=NO

PRINT ARG=avg_density FILE=COLVAR_REWEIGHT STRIDE=1
\endplumedfile

In the next example we calculate the ensemble average of the (particle) density at 500 K and 300 MPa (the simulation pressure was 1 bar).

\plumedfile
volume: READ FILE=COLVAR VALUES=volume  IGNORE_TIME
mybias: READ FILE=COLVAR VALUES=mybias.bias  IGNORE_TIME

# Shift volume (to avoid numerical issues)
rvol: COMBINE ARG=vol PARAMETERS=7.8 PERIODIC=NO

# Weights
bias_weights: REWEIGHT_BIAS TEMP=500 ARG=mybias.bias
temp_press_weights: REWEIGHT_TEMP_PRESS TEMP=500 PRESSURE=0.06022140857 REWEIGHT_PRESSURE=180.66422571 VOLUME=volume

# Ensemble average of the volume at 300 K and 300 MPa
avg_vol: AVERAGE ARG=volume LOGWEIGHTS=bias_weights,temp_press_weights
# Ensemble average of the density at 300 K and 300 MPa
avg_density: CUSTOM ARG=avg_vol FUNC=1000/x PERIODIC=NO

PRINT ARG=avg_density FILE=COLVAR_REWEIGHT STRIDE=1
\endplumedfile


In this final example we calculate the ensemble average of the (particle) density at 300 K and 300 MPa (the simulation temperature and pressure were 500 K and 1 bar).

\plumedfile
energy: READ FILE=COLVAR VALUES=energy  IGNORE_TIME
volume: READ FILE=COLVAR VALUES=volume  IGNORE_TIME
mybias: READ FILE=COLVAR VALUES=mybias.bias  IGNORE_TIME

# Shift energy and volume (to avoid numerical issues)
rvol: COMBINE ARG=vol PARAMETERS=7.8 PERIODIC=NO
renergy: COMBINE ARG=energy PARAMETERS=-13250 PERIODIC=NO

# Weights
bias_weights: REWEIGHT_BIAS TEMP=500 ARG=mybias.bias
temp_press_weights: REWEIGHT_TEMP_PRESS TEMP=500 REWEIGHT_TEMP=300 PRESSURE=0.06022140857 REWEIGHT_PRESSURE=180.66422571 ENERGY=renergy VOLUME=rvolume

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
  std::vector<Value*> myenergy, myvol;
public:
  static void registerKeywords(Keywords&);
  explicit ReweightTemperaturePressure(const ActionOptions&ao);
  double getLogWeight();
};

PLUMED_REGISTER_ACTION(ReweightTemperaturePressure,"REWEIGHT_TEMP_PRESS")

void ReweightTemperaturePressure::registerKeywords(Keywords& keys ) {
  ReweightBase::registerKeywords( keys );
  keys.remove("ARG");
  keys.add("optional","ENERGY","Energy");
  keys.add("optional","VOLUME","Volume");
  keys.add("optional","REWEIGHT_PRESSURE","Reweighting pressure");
  keys.add("optional","PRESSURE","The system pressure");
  keys.add("optional","REWEIGHT_TEMP","Reweighting temperature");
}

ReweightTemperaturePressure::ReweightTemperaturePressure(const ActionOptions&ao):
  Action(ao),
  ReweightBase(ao)
{
  // Initialize to not defined (negative)
  rpress_=-1;
  press_=-1;
  rtemp_=-1;
  parse("REWEIGHT_PRESSURE",rpress_);
  parse("PRESSURE",press_);
  parse("REWEIGHT_TEMP",rtemp_);
  rtemp_*=plumed.getAtoms().getKBoltzmann();

  parseArgumentList("ENERGY",myenergy);
  if(!myenergy.empty()) {
    log.printf("  with energies: ");
    for(unsigned i=0; i<myenergy.size(); i++) log.printf(" %s",myenergy[i]->getName().c_str());
    log.printf("\n");
  }
  //requestArguments(myenergy);

  parseArgumentList("VOLUME",myvol);
  if(!myvol.empty()) {
    log.printf("  with volumes: ");
    for(unsigned i=0; i<myvol.size(); i++) log.printf(" %s",myvol[i]->getName().c_str());
    log.printf("\n");
  }

  std::vector<Value*> conc;
  conc.insert(conc.begin(), myenergy.begin(), myenergy.end());
  conc.insert(conc.end(), myvol.begin(), myvol.end());
  requestArguments(conc);

  // 4 possible cases
  // Case 1) Reweight from T to T' with V=const (canonical)
  if (rtemp_>=0 && press_<0 && rpress_<0 && !myenergy.empty() && myvol.empty() ) {
    log.printf("  reweighting simulation from temperature %f to temperature %f at constant volume \n",simtemp/plumed.getAtoms().getKBoltzmann(),rtemp_/plumed.getAtoms().getKBoltzmann() );
    log.printf("  WARNING: If the simulation is performed at constant pressure add the keywords PRESSURE and VOLUME \n" );
  }
  // Case 2) Reweight from T to T' with P=const (isothermal-isobaric)
  else if (rtemp_>=0 && press_>=0 && rpress_<0 && !myenergy.empty() && !myvol.empty() ) log.printf("  reweighting simulation from temperature %f to temperature %f at constant pressure %f \n",simtemp/plumed.getAtoms().getKBoltzmann(),rtemp_/plumed.getAtoms().getKBoltzmann(), press_ );
  // Case 3) Reweight from P to P' with T=const (isothermal-isobaric)
  else if (rtemp_<0 && press_>=0 && rpress_>=0 && myenergy.empty() && !myvol.empty() ) log.printf("  reweighting simulation from pressure %f to pressure %f at constant temperature %f\n",press_,rpress_,simtemp/plumed.getAtoms().getKBoltzmann() );
  // Case 4) Reweight from T,P to T',P' (isothermal-isobaric)
  else if (rtemp_>0 && press_>=0 && rpress_>=0 && !myenergy.empty() && !myvol.empty() ) log.printf("  reweighting simulation from temperature %f and pressure %f to temperature %f and pressure %f \n",simtemp/plumed.getAtoms().getKBoltzmann(), press_, rtemp_/plumed.getAtoms().getKBoltzmann(), rpress_);
  else error("Combination of ENERGY, VOLUME, REWEIGHT_PRESSURE, PRESSURE and REWEIGHT_TEMP not supported. Please refer to the manual for supported combinations.");
}

double ReweightTemperaturePressure::getLogWeight() {
  double energy=0.0; for(unsigned i=0; i<myenergy.size(); ++i) energy+=getArgument(i);
  double volume=0.0; for(unsigned i=0; i<myvol.size(); ++i) volume+=getArgument(myenergy.size()+i);
  // 4 possible cases
  // Case 1) Reweight from T to T' with V=const (canonical)
  if (rtemp_>=0 && press_<0 && rpress_<0) return ((1.0/simtemp)- (1.0/rtemp_) )*energy;
  // Case 2) Reweight from T to T' with P=const (isothermal-isobaric)
  else if (rtemp_>=0 && press_>=0 && rpress_<0)  return ((1.0/simtemp)- (1.0/rtemp_) )*energy + ((1.0/simtemp) - (1.0/rtemp_))*press_*volume;
  // Case 3) Reweight from P to P' with T=const (isothermal-isobaric)
  else if (rtemp_<0 && press_>=0 && rpress_>=0)  return (1.0/simtemp)*(press_ - rpress_)*volume;
  // Case 4) Reweight from T,P to T',P' (isothermal-isobaric)
  else if (rtemp_>0 && press_>=0 && rpress_>=0) return ((1.0/simtemp)- (1.0/rtemp_) )*energy + ((1.0/simtemp)*press_ - (1.0/rtemp_)*rpress_ )*volume;
  else return 0;
}

}
}
