/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#include "Bias.h"
#include "ActionRegister.h"
#include "tools/Random.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

#include <iostream>


using namespace std;


namespace PLMD {
namespace bias {

//+PLUMEDOC BIAS EXTENDED_LAGRANGIAN
/*
Add extended Lagrangian.

This action can be used to create fictitious collective variables coupled to the real ones.
Given \f$x_i\f$ the \f$i\f$th argument of this bias potential, potential
and kinetic contributions are added to the energy of the system as
\f[
  V=\sum_i \frac{k_i}{2} (x_i-s_i)^2 + \sum_i \frac{\dot{s}_i^2}{2m_i}
\f].

The resulting potential is thus similar to a \ref RESTRAINT,
but the restraint center moved with time following Hamiltonian
dynamics with mass \f$m_i\f$.

This bias potential accepts thus vectorial keywords (one element per argument)
to define the coupling constant (KAPPA) and a relaxation time \f$tau\f$ (TAU).
The mass is them computed as \f$m=k(\frac{\tau}{2\pi})^2\f$.

Notice that this action creates several components.
The ones named XX_fict are the fictitious coordinates. It is possible
to add further forces on them by means of other bias potential,
e.g. to obtain an indirect \ref METAD as in \cite continua .
Also notice that the velocities of the fictitious coordinates
are reported (XX_vfict). However, printed velocities are the ones
at the previous step.

It is also possible to provide a non-zero friction (one value per component).
This is then used to implement a Langevin thermostat, so as to implement
TAMD/dAFED method \cite Maragliano2006 \cite AbramsJ2008 . Notice that
here a massive Langevin thermostat is used, whereas usually
TAMD employs an overamped Langevin dynamics and dAFED
a Gaussian thermostat.

\warning
The bias potential is reported in the component bias.
Notice that this bias potential, although formally compatible with
replica exchange framework, probably does not work as expected in that case.
Indeed, since fictitious coordinates are not swapped upon exchange,
acceptace can be expected to be extremely low unless (by chance) two neighboring
replicas have the fictitious variables located properly in space.

\warning
\ref RESTART is not properly supported by this action. Indeed,
at every start the postion of the fictitious variable is reset to the value
of the real variable, and its velocity is set to zero.
This is not expected to introduce big errors, but certainly is
introducing a small inconsistency between a single long run
and many shorter runs.

\par Examples

The following input tells plumed to perform a metadynamics
with an extended Lagrangian on two torsional angles.
\plumedfile
phi: TORSION ATOMS=5,7,9,15
psi: TORSION ATOMS=7,9,15,17
ex: EXTENDED_LAGRANGIAN ARG=phi,psi KAPPA=20,20.0 TAU=0.1,0.1
METAD ARG=ex.phi_fict,ex.psi_fict PACE=100 SIGMA=0.35,0.35 HEIGHT=0.1
# monitor the two variables
PRINT STRIDE=10 ARG=phi,psi,ex.phi_fict,ex.psi_fict FILE=COLVAR
\endplumedfile

The following input tells plumed to perform a TAMD (or dAFED)
calculation on two torsional angles, keeping the two variables
at a fictitious temperature of 3000K with a Langevin thermostat
with friction 10
\plumedfile
phi: TORSION ATOMS=5,7,9,15
psi: TORSION ATOMS=7,9,15,17
ex: EXTENDED_LAGRANGIAN ARG=phi,psi KAPPA=20,20.0 TAU=0.1,0.1 FRICTION=10,10 TEMP=3000
# monitor the two variables
PRINT STRIDE=10 ARG=phi,psi,ex.phi_fict,ex.psi_fict FILE=COLVAR
\endplumedfile

*/
//+ENDPLUMEDOC

class ExtendedLagrangian : public Bias {
  bool firsttime;
  std::vector<double> fict;
  std::vector<double> vfict;
  std::vector<double> vfict_laststep;
  std::vector<double> ffict;
  std::vector<double> kappa;
  std::vector<double> tau;
  std::vector<double> friction;
  std::vector<Value*> fictValue;
  std::vector<Value*> vfictValue;
  double kbt;
  Random rand;
public:
  explicit ExtendedLagrangian(const ActionOptions&);
  void calculate() override;
  void update() override;
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(ExtendedLagrangian,"EXTENDED_LAGRANGIAN")

void ExtendedLagrangian::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","KAPPA","specifies that the restraint is harmonic and what the values of the force constants on each of the variables are");
  keys.add("compulsory","TAU","specifies that the restraint is harmonic and what the values of the force constants on each of the variables are");
  keys.add("compulsory","FRICTION","0.0","add a friction to the variable");
  keys.add("optional","TEMP","the system temperature - needed when FRICTION is present. If not provided will be taken from MD code (if available)");
  componentsAreNotOptional(keys);
  keys.addOutputComponent("_fict","default","one or multiple instances of this quantity can be referenced elsewhere in the input file. "
                          "These quantities will named with the arguments of the bias followed by "
                          "the character string _tilde. It is possible to add forces on these variable.");
  keys.addOutputComponent("_vfict","default","one or multiple instances of this quantity can be referenced elsewhere in the input file. "
                          "These quantities will named with the arguments of the bias followed by "
                          "the character string _tilde. It is NOT possible to add forces on these variable.");
}

ExtendedLagrangian::ExtendedLagrangian(const ActionOptions&ao):
  PLUMED_BIAS_INIT(ao),
  firsttime(true),
  fict(getNumberOfArguments(),0.0),
  vfict(getNumberOfArguments(),0.0),
  vfict_laststep(getNumberOfArguments(),0.0),
  ffict(getNumberOfArguments(),0.0),
  kappa(getNumberOfArguments(),0.0),
  tau(getNumberOfArguments(),0.0),
  friction(getNumberOfArguments(),0.0),
  fictValue(getNumberOfArguments(),NULL),
  vfictValue(getNumberOfArguments(),NULL),
  kbt(0.0)
{
  parseVector("TAU",tau);
  parseVector("FRICTION",friction);
  parseVector("KAPPA",kappa);
  double temp=-1.0;
  parse("TEMP",temp);
  if(temp>=0.0) kbt=plumed.getAtoms().getKBoltzmann()*temp;
  else kbt=plumed.getAtoms().getKbT();
  checkRead();

  log.printf("  with harmonic force constant");
  for(unsigned i=0; i<kappa.size(); i++) log.printf(" %f",kappa[i]);
  log.printf("\n");

  log.printf("  with relaxation time");
  for(unsigned i=0; i<tau.size(); i++) log.printf(" %f",tau[i]);
  log.printf("\n");

  bool hasFriction=false;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) if(friction[i]>0.0) hasFriction=true;

  if(hasFriction) {
    log.printf("  with friction");
    for(unsigned i=0; i<friction.size(); i++) log.printf(" %f",friction[i]);
    log.printf("\n");
  }

  log.printf("  and kbt");
  log.printf(" %f",kbt);
  log.printf("\n");

  for(unsigned i=0; i<getNumberOfArguments(); i++) {
    std::string comp=getPntrToArgument(i)->getName()+"_fict";
    addComponentWithDerivatives(comp);
    if(getPntrToArgument(i)->isPeriodic()) {
      std::string a,b;
      getPntrToArgument(i)->getDomain(a,b);
      componentIsPeriodic(comp,a,b);
    } else componentIsNotPeriodic(comp);
    fictValue[i]=getPntrToComponent(comp);
    comp=getPntrToArgument(i)->getName()+"_vfict";
    addComponent(comp);
    componentIsNotPeriodic(comp);
    vfictValue[i]=getPntrToComponent(comp);
  }

  log<<"  Bibliography "<<plumed.cite("Iannuzzi, Laio, and Parrinello, Phys. Rev. Lett. 90, 238302 (2003)");
  if(hasFriction) {
    log<<plumed.cite("Maragliano and Vanden-Eijnden, Chem. Phys. Lett. 426, 168 (2006)");
    log<<plumed.cite("Abrams and Tuckerman, J. Phys. Chem. B 112, 15742 (2008)");
  }
  log<<"\n";
}


void ExtendedLagrangian::calculate() {

  if(firsttime) {
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      fict[i]=getArgument(i);
    }
    firsttime=false;
  }
  double ene=0.0;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    const double cv=difference(i,fict[i],getArgument(i));
    const double k=kappa[i];
    const double f=-k*cv;
    ene+=0.5*k*cv*cv;
    setOutputForce(i,f);
    ffict[i]=-f;
  };
  setBias(ene);
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    fict[i]=fictValue[i]->bringBackInPbc(fict[i]);
    fictValue[i]->set(fict[i]);
    vfictValue[i]->set(vfict_laststep[i]);
  }
}

void ExtendedLagrangian::update() {
  double dt=getTimeStep()*getStride();
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    double mass=kappa[i]*tau[i]*tau[i]/(4*pi*pi); // should be k/omega**2
    double c1=exp(-0.5*friction[i]*dt);
    double c2=sqrt(kbt*(1.0-c1*c1)/mass);
// consider additional forces on the fictitious particle
// (e.g. MetaD stuff)
    ffict[i]+=fictValue[i]->getForce();

// update velocity (half step)
    vfict[i]+=ffict[i]*0.5*dt/mass;
// thermostat (half step)
    vfict[i]=c1*vfict[i]+c2*rand.Gaussian();
// save full step velocity to be dumped at next step
    vfict_laststep[i]=vfict[i];
// thermostat (half step)
    vfict[i]=c1*vfict[i]+c2*rand.Gaussian();
// update velocity (half step)
    vfict[i]+=ffict[i]*0.5*dt/mass;
// update position (full step)
    fict[i]+=vfict[i]*dt;
  }
}

}

}
