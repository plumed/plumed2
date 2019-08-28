/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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
#include "Colvar.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR DIPOLE
/*
Calculate the dipole moment for a group of atoms.

When running with periodic boundary conditions, the atoms should be
in the proper periodic image. This is done automatically since PLUMED 2.5,
by considering the ordered list of atoms and rebuilding the molecule with a procedure
that is equivalent to that done in \ref WHOLEMOLECULES . Notice that
rebuilding is local to this action. This is different from \ref WHOLEMOLECULES
which actually modifies the coordinates stored in PLUMED.

In case you want to recover the old behavior you should use the NOPBC flag.
In that case you need to take care that atoms are in the correct
periodic image.

\par Examples

The following tells plumed to calculate the dipole of the group of atoms containing
the atoms from 1-10 and print it every 5 steps
\plumedfile
d: DIPOLE GROUP=1-10
PRINT FILE=output STRIDE=5 ARG=d
\endplumedfile

\attention
If the total charge Q of the group in non zero, then a charge Q/N will be subtracted to every atom,
where N is the number of atoms. This implies that the dipole (which for a charged system depends
on the position) is computed on the geometric center of the group.


*/
//+ENDPLUMEDOC

class Dipole : public Colvar {
  vector<AtomNumber> ga_lista;
  bool components;
  bool nopbc;
public:
  explicit Dipole(const ActionOptions&);
  void calculate() override;
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(Dipole,"DIPOLE")

void Dipole::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.add("atoms","GROUP","the group of atoms we are calculating the dipole moment for");
  keys.addFlag("COMPONENTS",false,"calculate the x, y and z components of the dipole separately and store them as label.x, label.y and label.z");
  keys.addOutputComponent("x","COMPONENTS","the x-component of the dipole");
  keys.addOutputComponent("y","COMPONENTS","the y-component of the dipole");
  keys.addOutputComponent("z","COMPONENTS","the z-component of the dipole");
}

Dipole::Dipole(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  components(false)
{
  parseAtomList("GROUP",ga_lista);
  parseFlag("COMPONENTS",components);
  parseFlag("NOPBC",nopbc);
  checkRead();
  if(components) {
    addComponentWithDerivatives("x"); componentIsNotPeriodic("x");
    addComponentWithDerivatives("y"); componentIsNotPeriodic("y");
    addComponentWithDerivatives("z"); componentIsNotPeriodic("z");
  } else {
    addValueWithDerivatives(); setNotPeriodic();
  }

  log.printf("  of %u atoms\n",static_cast<unsigned>(ga_lista.size()));
  for(unsigned int i=0; i<ga_lista.size(); ++i) {
    log.printf("  %d", ga_lista[i].serial());
  }
  log.printf("  \n");
  if(nopbc) log.printf("  without periodic boundary conditions\n");
  else      log.printf("  using periodic boundary conditions\n");

  requestAtoms(ga_lista);
}

// calculator
void Dipole::calculate()
{
  if(!nopbc) makeWhole();
  double ctot=0.;
  unsigned N=getNumberOfAtoms();
  vector<double> charges(N);
  Vector dipje;

  for(unsigned i=0; i<N; ++i) {
    charges[i]=getCharge(i);
    ctot+=charges[i];
  }
  ctot/=(double)N;

  for(unsigned i=0; i<N; ++i) {
    charges[i]-=ctot;
    dipje += charges[i]*getPosition(i);
  }

  if(!components) {
    double dipole = dipje.modulo();
    double idip = 1./dipole;

    for(unsigned i=0; i<N; i++) {
      double dfunc=charges[i]*idip;
      setAtomsDerivatives(i,dfunc*dipje);
    }
    setBoxDerivativesNoPbc();
    setValue(dipole);
  } else {
    Value* valuex=getPntrToComponent("x");
    Value* valuey=getPntrToComponent("y");
    Value* valuez=getPntrToComponent("z");
    for(unsigned i=0; i<N; i++) {
      setAtomsDerivatives(valuex,i,charges[i]*Vector(1.0,0.0,0.0));
      setAtomsDerivatives(valuey,i,charges[i]*Vector(0.0,1.0,0.0));
      setAtomsDerivatives(valuez,i,charges[i]*Vector(0.0,0.0,1.0));
    }
    setBoxDerivativesNoPbc(valuex);
    setBoxDerivativesNoPbc(valuey);
    setBoxDerivativesNoPbc(valuez);
    valuex->set(dipje[0]);
    valuey->set(dipje[1]);
    valuez->set(dipje[2]);
  }
}

}
}
