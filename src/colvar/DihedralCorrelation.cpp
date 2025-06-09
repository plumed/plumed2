/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#include "tools/Torsion.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR DIHEDRAL_CORRELATION
/*
Measure the correlation between a pair of dihedral angles

This CV measures the correlation between two dihedral angles, $\phi$ and $\psi$, as follows:

$$
s = \frac{1}{2} \left[ 1 + \cos( \phi - \psi ) \right]
$$

An example input that computes and calculates this quantity with $\phi$ being the dihedral
angle involving atoms 1, 2, 3 and 4 and $\psi$ being the angle involving atoms 7, 8, 9 and 10
is shown below:

```plumed
d: DIHEDRAL_CORRELATION ATOMS=1,2,3,4,5,6,7,8
PRINT ARG=d FILE=colvar
```

If you want to calculate the dihedral correlations between multiple pairs of dihedral angles using
this action you would use an input like this one shown below:

```plumed
d: DIHEDRAL_CORRELATION ...
  ATOMS1=1,2,3,4,5,6,7,8
  ATOMS2=9,10,11,12,13,14,15,16
...
PRINT ARG=d FILE=colvar
```

This input calculates and outputs a two dimensional vector that contains two of these dihedral correlation
values. Commands similar to these are used within the [DIHCOR](DIHCOR.md) shortcut.

The last thing to note is that by default a procedure akin to that used in [WHOLEMOLECULES](WHOLEMOLECULES.md)
is used to ensure that the sets of atoms that are specified to each ATOMS keyword are not broken by the periodic
boundary conditions.  If you would like to turn this off for any reason you add the NOPBC in your input file as shown
below:

```plumed
d: DIHEDRAL_CORRELATION ATOMS=1,2,3,4,5,6,7,8 NOPBC
PRINT ARG=d FILE=colvar
```

*/
//+ENDPLUMEDOC

class DihedralCorrelation : public Colvar {
private:
  bool pbc;
  std::vector<double> value;
  std::vector<double> derivs;
public:
  static void registerKeywords( Keywords& keys );
  explicit DihedralCorrelation(const ActionOptions&);
  static void parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa );
  static unsigned getModeAndSetupValues( ActionWithValue* av );
  void calculate() override;
  static void calculateCV( const ColvarInput& cvin, ColvarOutput& cvout );
};

typedef ColvarShortcut<DihedralCorrelation> DihedralCorrelationShortcut;
PLUMED_REGISTER_ACTION(DihedralCorrelationShortcut,"DIHEDRAL_CORRELATION")
PLUMED_REGISTER_ACTION(DihedralCorrelation,"DIHEDRAL_CORRELATION_SCALAR")
typedef MultiColvarTemplate<DihedralCorrelation> DihedralCorrelationMulti;
PLUMED_REGISTER_ACTION(DihedralCorrelationMulti,"DIHEDRAL_CORRELATION_VECTOR")

void DihedralCorrelation::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.setDisplayName("DIHEDRAL_CORRELATION");
  keys.add("atoms","ATOMS","the set of 8 atoms that are being used to calculate this quantity");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  keys.setValueDescription("scalar/vector","the DIHEDRAL_CORRELATION for these atoms");
  keys.reset_style("NUMERICAL_DERIVATIVES","hidden");
}

DihedralCorrelation::DihedralCorrelation(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  value(1) {
  std::vector<AtomNumber> atoms;
  parseAtomList(-1,atoms,this);
  if( atoms.size()!=8 ) {
    error("Number of specified atoms should be 8");
  }

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
}

void DihedralCorrelation::parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa ) {
  aa->parseAtomList("ATOMS",num,t);
  if( num<0 && t.size()!=8 ) {
    aa->error("Number of specified atoms should be 8");
  }
  if( t.size()==8 ) {
    aa->log.printf("  correlation between dihedral angle for atoms %d %d %d %d and atoms %d %d %d %d\n",
                   t[0].serial(),t[1].serial(),t[2].serial(),t[3].serial(),t[4].serial(),t[5].serial(),t[6].serial(),t[7].serial());
  }
}

unsigned DihedralCorrelation::getModeAndSetupValues( ActionWithValue* av ) {
  av->addValueWithDerivatives();
  av->setNotPeriodic();
  return 0;
}

void DihedralCorrelation::calculate() {

  if(pbc) {
    makeWhole();
  }
  ColvarOutput cvout = ColvarOutput::createColvarOutput(value,derivs,this);
  calculateCV( ColvarInput::createColvarInput( 0, getPositions(), this ), cvout );
  setValue( value[0] );
  for(unsigned i=0; i<getPositions().size(); ++i) {
    setAtomsDerivatives( i, cvout.getAtomDerivatives(0, i) );
  }
  setBoxDerivatives( cvout.virial[0] );
}

void DihedralCorrelation::calculateCV( const ColvarInput& cvin, ColvarOutput& cvout ) {
  const Vector d10=delta(cvin.pos[1],cvin.pos[0]);
  const Vector d11=delta(cvin.pos[2],cvin.pos[1]);
  const Vector d12=delta(cvin.pos[3],cvin.pos[2]);

  Vector dd10,dd11,dd12;
  PLMD::Torsion t1;
  const double phi1  = t1.compute(d10,d11,d12,dd10,dd11,dd12);

  const Vector d20=delta(cvin.pos[5],cvin.pos[4]);
  const Vector d21=delta(cvin.pos[6],cvin.pos[5]);
  const Vector d22=delta(cvin.pos[7],cvin.pos[6]);

  Vector dd20,dd21,dd22;
  PLMD::Torsion t2;
  const double phi2 = t2.compute( d20, d21, d22, dd20, dd21, dd22 );

  // Calculate value
  const double diff = phi2 - phi1;
  cvout.values[0] = 0.5*(1.+std::cos(diff));
  // Derivatives wrt phi1
  const double dval = 0.5*std::sin(diff);
  dd10 *= dval;
  dd11 *= dval;
  dd12 *= dval;
  // And add
  cvout.derivs[0][0]=dd10;
  cvout.derivs[0][1]=dd11-dd10;
  cvout.derivs[0][2]=dd12-dd11;
  cvout.derivs[0][3]=-dd12;
  // Derivative wrt phi2
  dd20 *= -dval;
  dd21 *= -dval;
  dd22 *= -dval;
  // And add
  cvout.derivs[0][4]=dd20;
  cvout.derivs[0][5]=dd21-dd20;
  cvout.derivs[0][6]=dd22-dd21;
  cvout.derivs[0][7]=-dd22;
  cvout.virial.set( 0,
                    -(extProduct(d10,dd10)+extProduct(d11,dd11)+extProduct(d12,dd12))
                    - (extProduct(d20,dd20)+extProduct(d21,dd21)+extProduct(d22,dd22)) );
}

}
}
