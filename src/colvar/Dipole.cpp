/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include "ColvarShortcut.h"
#include "MultiColvarTemplate.h"
#include "core/ActionRegister.h"

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

//+PLUMEDOC COLVAR DIPOLE_SCALAR
/*
Calculate the dipole moment for a group of atoms.

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR DIPOLE_VECTOR
/*
Calculate a vector of dipole moments for a set of groups of atoms.

\par Examples

*/
//+ENDPLUMEDOC

class Dipole : public Colvar {
  std::vector<AtomNumber> ga_lista;
  bool components;
  bool nopbc;
  std::vector<double> value, masses, charges;
  std::vector<std::vector<Vector> > derivs;
  std::vector<Tensor> virial;
  Value* valuex=nullptr;
  Value* valuey=nullptr;
  Value* valuez=nullptr;
public:
  explicit Dipole(const ActionOptions&);
  static void parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa );
  static unsigned getModeAndSetupValues( ActionWithValue* av );
  void calculate() override;
  static void registerKeywords(Keywords& keys);
  static void calculateCV( const unsigned& mode, const std::vector<double>& masses, std::vector<double>& charges,
                           const std::vector<Vector>& pos, std::vector<double>& vals, std::vector<std::vector<Vector> >& derivs,
                           std::vector<Tensor>& virial, const ActionAtomistic* aa );
};

typedef ColvarShortcut<Dipole> DipoleShortcut;
PLUMED_REGISTER_ACTION(DipoleShortcut,"DIPOLE")
PLUMED_REGISTER_ACTION(Dipole,"DIPOLE_SCALAR")
typedef MultiColvarTemplate<Dipole> DipoleMulti;
PLUMED_REGISTER_ACTION(DipoleMulti,"DIPOLE_VECTOR")

void Dipole::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.setDisplayName("DIPOLE");
  keys.add("atoms","GROUP","the group of atoms we are calculating the dipole moment for");
  keys.addFlag("COMPONENTS",false,"calculate the x, y and z components of the dipole separately and store them as label.x, label.y and label.z");
  keys.addOutputComponent("x","COMPONENTS","the x-component of the dipole");
  keys.addOutputComponent("y","COMPONENTS","the y-component of the dipole");
  keys.addOutputComponent("z","COMPONENTS","the z-component of the dipole");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  keys.setValueDescription("the DIPOLE for these atoms");
}

Dipole::Dipole(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  components(false),
  value(1),
  derivs(1),
  virial(1) {
  parseAtomList(-1,ga_lista,this);
  charges.resize(ga_lista.size());
  components=(getModeAndSetupValues(this)==1);
  if( components ) {
    value.resize(3);
    derivs.resize(3);
    virial.resize(3);
    valuex=getPntrToComponent("x");
    valuey=getPntrToComponent("y");
    valuez=getPntrToComponent("z");
  }
  for(unsigned i=0; i<derivs.size(); ++i) {
    derivs[i].resize( ga_lista.size() );
  }
  parseFlag("NOPBC",nopbc);
  checkRead();

  if(nopbc) {
    log.printf("  without periodic boundary conditions\n");
  } else {
    log.printf("  using periodic boundary conditions\n");
  }

  requestAtoms(ga_lista);
}

void Dipole::parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa ) {
  aa->parseAtomList("GROUP",num,t);
  if( t.size()>0 ) {
    aa->log.printf("  of %u atoms\n",static_cast<unsigned>(t.size()));
    for(unsigned int i=0; i<t.size(); ++i) {
      aa->log.printf("  %d", t[i].serial());
    }
    aa->log.printf("  \n");
  }
}

unsigned Dipole::getModeAndSetupValues( ActionWithValue* av ) {
  bool c;
  av->parseFlag("COMPONENTS",c);
  if( c ) {
    av->addComponentWithDerivatives("x");
    av->componentIsNotPeriodic("x");
    av->addComponentWithDerivatives("y");
    av->componentIsNotPeriodic("y");
    av->addComponentWithDerivatives("z");
    av->componentIsNotPeriodic("z");
    return 1;
  }
  av->addValueWithDerivatives();
  av->setNotPeriodic();
  return 0;
}

// calculator
void Dipole::calculate() {
  if(!nopbc) {
    makeWhole();
  }
  unsigned N=getNumberOfAtoms();
  for(unsigned i=0; i<N; ++i) {
    charges[i]=getCharge(i);
  }

  if(!components) {
    calculateCV( 0, masses, charges, getPositions(), value, derivs, virial, this );
    for(unsigned i=0; i<N; i++) {
      setAtomsDerivatives(i,derivs[0][i]);
    }
    setBoxDerivatives(virial[0]);
    setValue(value[0]);
  } else {
    calculateCV( 1, masses, charges, getPositions(), value, derivs, virial, this );
    for(unsigned i=0; i<N; i++) {
      setAtomsDerivatives(valuex,i,derivs[0][i]);
      setAtomsDerivatives(valuey,i,derivs[1][i]);
      setAtomsDerivatives(valuez,i,derivs[2][i]);
    }
    setBoxDerivatives(valuex,virial[0]);
    setBoxDerivatives(valuey,virial[1]);
    setBoxDerivatives(valuez,virial[2]);
    valuex->set(value[0]);
    valuey->set(value[1]);
    valuez->set(value[2]);
  }
}

void Dipole::calculateCV( const unsigned& mode, const std::vector<double>& masses, std::vector<double>& charges,
                          const std::vector<Vector>& pos, std::vector<double>& vals, std::vector<std::vector<Vector> >& derivs,
                          std::vector<Tensor>& virial, const ActionAtomistic* aa ) {
  unsigned N=pos.size();
  double ctot=0.;
  for(unsigned i=0; i<N; ++i) {
    ctot += charges[i];
  }
  ctot/=(double)N;

  Vector dipje;
  for(unsigned i=0; i<N; ++i) {
    charges[i]-=ctot;
    dipje += charges[i]*pos[i];
  }

  if( mode==1 ) {
    for(unsigned i=0; i<N; i++) {
      derivs[0][i]=charges[i]*Vector(1.0,0.0,0.0);
      derivs[1][i]=charges[i]*Vector(0.0,1.0,0.0);
      derivs[2][i]=charges[i]*Vector(0.0,0.0,1.0);
    }
    for(unsigned i=0; i<3; ++i ) {
      vals[i] = dipje[i];
    }
  } else {
    vals[0] = dipje.modulo();
    double idip = 1./vals[0];
    for(unsigned i=0; i<N; i++) {
      double dfunc=charges[i]*idip;
      derivs[0][i] = dfunc*dipje;
    }
  }
  setBoxDerivativesNoPbc( pos, derivs, virial );
}

}
}
