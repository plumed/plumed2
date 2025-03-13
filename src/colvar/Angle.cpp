/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#ifdef __PLUMED_HAS_OPENACC
#define __PLUMED_USE_OPENACC 1
#endif //__PLUMED_HAS_OPENACC
#include "Colvar.h"
#include "ColvarShortcut.h"
#include "MultiColvarTemplate.h"
#include "core/ActionRegister.h"
#include "tools/Angle.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR ANGLE
/*
Calculate an angle.

This command can be used to compute the angle between three atoms. Alternatively
if four atoms appear in the atom
specification it calculates the angle between
two vectors identified by two pairs of atoms.

If _three_ atoms are given, the angle is defined as:
\f[
\theta=\arccos\left(\frac{ {\bf r}_{21}\cdot {\bf r}_{23}}{
|{\bf r}_{21}| |{\bf r}_{23}|}\right)
\f]
Here \f$ {\bf r}_{ij}\f$ is the distance vector among the
\f$i\f$th and the \f$j\f$th listed atom.

If _four_ atoms are given, the angle is defined as:
\f[
\theta=\arccos\left(\frac{ {\bf r}_{21}\cdot {\bf r}_{34}}{
|{\bf r}_{21}| |{\bf r}_{34}|}\right)
\f]

Notice that angles defined in this way are non-periodic variables and
their value is limited by definition between 0 and \f$\pi\f$.

The vectors \f$ {\bf r}_{ij}\f$ are by default evaluated taking
periodic boundary conditions into account.
This behavior can be changed with the NOPBC flag.

\par Examples

This command tells plumed to calculate the angle between the vector connecting atom 1 to atom 2 and
the vector connecting atom 2 to atom 3 and to print it on file COLVAR1. At the same time,
the angle between vector connecting atom 1 to atom 2 and the vector connecting atom 3 to atom 4 is printed
on file COLVAR2.
\plumedfile

a: ANGLE ATOMS=1,2,3
# equivalently one could state:
# a: ANGLE ATOMS=1,2,2,3

b: ANGLE ATOMS=1,2,3,4

PRINT ARG=a FILE=COLVAR1
PRINT ARG=b FILE=COLVAR2
\endplumedfile


*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR ANGLE_SCALAR
/*
Calculate an angle.

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR ANGLE_VECTOR
/*
Calculate multiple angles.

\par Examples

*/
//+ENDPLUMEDOC

class Angle : public Colvar {
  bool pbc;
  std::vector<double> value;
  std::vector<double> derivs;
public:
  explicit Angle(const ActionOptions&);
// active methods:
  void calculate() override;
  static void registerKeywords( Keywords& keys );
  static void parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa );
  static unsigned getModeAndSetupValues( ActionWithValue* av );
  static void calculateCV( const ColvarInput& cvin, ColvarOutput& cvout );
};

typedef ColvarShortcut<Angle> AngleShortcut;
PLUMED_REGISTER_ACTION(AngleShortcut,"ANGLE")
PLUMED_REGISTER_ACTION(Angle,"ANGLE_SCALAR")
typedef MultiColvarTemplate<Angle> AngleMulti;
PLUMED_REGISTER_ACTION(AngleMulti,"ANGLE_VECTOR")

void Angle::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords(keys);
  keys.setDisplayName("ANGLE");
  keys.add("atoms","ATOMS","the list of atoms involved in this collective variable (either 3 or 4 atoms)");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  keys.setValueDescription("scalar/vector","the ANGLE involving these atoms");
}

void Angle::parseAtomList( const int& num, std::vector<AtomNumber>& atoms, ActionAtomistic* aa ) {
  aa->parseAtomList("ATOMS",num,atoms);
  if(atoms.size()==3) {
    aa->log.printf("  between atoms %d %d %d\n",atoms[0].serial(),atoms[1].serial(),atoms[2].serial());
    atoms.resize(4);
    atoms[3]=atoms[2];
    atoms[2]=atoms[1];
  } else if(atoms.size()==4) {
    aa->log.printf("  between lines %d-%d and %d-%d\n",atoms[0].serial(),atoms[1].serial(),atoms[2].serial(),atoms[3].serial());
  } else if( num<0 || atoms.size()>0 ) {
    aa->error("Number of specified atoms should be either 3 or 4");
  }
}

unsigned Angle::getModeAndSetupValues( ActionWithValue* av ) {
  av->addValueWithDerivatives();
  av->setNotPeriodic();
  return 0;
}

Angle::Angle(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  value(1) {
  std::vector<AtomNumber> atoms;
  parseAtomList( -1, atoms, this );
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  if(pbc) {
    log.printf("  using periodic boundary conditions\n");
  } else {
    log.printf("  without periodic boundary conditions\n");
  }

  addValueWithDerivatives();
  setNotPeriodic();
  requestAtoms(atoms);
  checkRead();
}

// calculator
void Angle::calculate() {

  if(pbc) {
    makeWhole();
  }
  ColvarOutput cvout = ColvarOutput::createColvarOutput(value,derivs,this);
  calculateCV( ColvarInput::createColvarInput( 0, getPositions(), this ), cvout );
  setValue( value[0] );
  for(unsigned i=0; i<getPositions().size(); ++i) {
    setAtomsDerivatives( i, cvout.getAtomDerivatives(0,i) );
  }
  setBoxDerivatives( cvout.virial[0] );
}

void Angle::calculateCV( const ColvarInput& cvin, ColvarOutput& cvout ) {
  Vector dij,dik;
  dij=delta(cvin.pos[2],cvin.pos[3]);
  dik=delta(cvin.pos[1],cvin.pos[0]);
  Vector ddij,ddik;
  PLMD::Angle a;
  cvout.values[0]=a.compute(dij,dik,ddij,ddik);
  cvout.derivs[0][0]=ddik;
  cvout.derivs[0][1]=-ddik;
  cvout.derivs[0][2]=-ddij;
  cvout.derivs[0][3]=ddij;
  ColvarInput::setBoxDerivativesNoPbc( cvin, cvout );
}

}
}



