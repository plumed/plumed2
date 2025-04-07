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
#include "Colvar.h"
#include "ColvarShortcut.h"
#include "MultiColvarTemplate.h"
#include "core/ActionRegister.h"
#include "tools/Torsion.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR TORSION
/*
Calculate one or multiple torsional angles.

This command can be used to compute the torsion between four atoms as shown by the input below:

```plumed
t: TORSION ATOMS=1,2,3,4
PRINT ARG=t FILE=COLVAR
```

Alternatively you can use this action to calculate the angle between two vectors projected on the plane
orthogonal to an axis.  The example below uses this syntax and computes the same torsion as the first example
input above.

```plumed
t: TORSION VECTOR1=2,1 AXIS=2,3 VECTOR2=3,4
PRINT ARG=t FILE=COLVAR
```

If you combine this sytax with the functionality in [FIXEDATOM](FIXEDATOM.md) you can see how we can calculate the
torsional angle between two bond vectors around the z-axis as shown below:

```plumed
a0: FIXEDATOM AT=0,0,0
az: FIXEDATOM AT=0,0,1
t1: TORSION VECTOR1=1,2 AXIS=a0,az VECTOR2=5,6
PRINT ARG=t1 FILE=colvar STRIDE=20
```

If you are working with a protein you can specify the special named torsion angles $\phi$, $\psi$, $\omega$ and $\chi_1$
by using TORSION in combination with the [MOLINFO](MOLINFO.md) command.  This can be done by using the following
syntax.

```plumed
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO MOLTYPE=protein STRUCTURE=regtest/basic/rt32/helix.pdb
t1: TORSION ATOMS=@phi-3
t2: TORSION ATOMS=@psi-4
PRINT ARG=t1,t2 FILE=colvar STRIDE=10
```

Here, `@phi-3` tells plumed that you would like to calculate the $\phi$ angle in the third residue of the protein.
Similarly `@psi-4` tells plumed that you want to calculate the $\psi$ angle of the fourth residue of the protein.

If you would like to calculate multiple torsion angles at the same time you can use a command like the one shown below:

```plumed
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO MOLTYPE=protein STRUCTURE=regtest/basic/rt32/helix.pdb
t1: TORSION ATOMS1=@phi-3 ATOMS2=@phi-4 ATOMS3=@phi-5 ATOMS4=@phi-6 ATOMS5=@phi-7
PRINT ARG=t1 FILE=colvar STRIDE=10
```

This input tells PLUMED to calculate the $\phi$ angles in residues 3-7 of the protein.  The output, `t1`, is a 5 dimensional vector.


*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR TORSION_SCALAR
/*
Calculate a torsional angle.

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR TORSION_VECTOR
/*
Calculate multiple torsional angles.

\par Examples

*/
//+ENDPLUMEDOC

class Torsion : public Colvar {
  bool pbc;
  bool do_cosine;

  std::vector<double> value, masses, charges;
  std::vector<std::vector<Vector> > derivs;
  std::vector<Tensor> virial;
public:
  explicit Torsion(const ActionOptions&);
  static void parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa );
  static unsigned getModeAndSetupValues( ActionWithValue* av );
// active methods:
  void calculate() override;
  static void calculateCV( const unsigned& mode, const std::vector<double>& masses, const std::vector<double>& charges,
                           const std::vector<Vector>& pos, std::vector<double>& vals, std::vector<std::vector<Vector> >& derivs,
                           std::vector<Tensor>& virial, const ActionAtomistic* aa );
  static void registerKeywords(Keywords& keys);
};

typedef ColvarShortcut<Torsion> TorsionShortcut;
PLUMED_REGISTER_ACTION(TorsionShortcut,"TORSION")
PLUMED_REGISTER_ACTION(Torsion,"TORSION_SCALAR")
typedef MultiColvarTemplate<Torsion> TorsionMulti;
PLUMED_REGISTER_ACTION(TorsionMulti,"TORSION_VECTOR")

void Torsion::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords( keys );
  keys.setDisplayName("TORSION");
  keys.add("atoms-1","ATOMS","the four atoms involved in the torsional angle");
  keys.add("atoms-2","AXIS","two atoms that define an axis.  You can use this to find the angle in the plane perpendicular to the axis between the vectors specified using the VECTORA and VECTORB keywords.");
  keys.add("atoms-2","VECTORA","two atoms that define a vector.  You can use this in combination with VECTOR2 and AXIS");
  keys.add("atoms-2","VECTORB","two atoms that define a vector.  You can use this in combination with VECTOR1 and AXIS");
  keys.add("atoms-3","VECTOR1","two atoms that define a vector.  You can use this in combination with VECTOR2 and AXIS");
  keys.add("atoms-3","VECTOR2","two atoms that define a vector.  You can use this in combination with VECTOR1 and AXIS");
  keys.addFlag("COSINE",false,"calculate cosine instead of dihedral");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  keys.setValueDescription("scalar/vector","the TORSION involving these atoms");
}

Torsion::Torsion(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  do_cosine(false),
  value(1),
  derivs(1),
  virial(1) {
  derivs[0].resize(6);
  std::vector<AtomNumber> atoms;
  std::vector<AtomNumber> v1;
  ActionAtomistic::parseAtomList("VECTOR1",v1);
  if( v1.size()>0 ) {
    std::vector<AtomNumber> v2;
    ActionAtomistic::parseAtomList("VECTOR2",v2);
    std::vector<AtomNumber> axis;
    ActionAtomistic::parseAtomList("AXIS",axis);
    if( !(v1.size()==2 && v2.size()==2 && axis.size()==2)) {
      error("VECTOR1, VECTOR2 and AXIS should specify 2 atoms each");
    }
    atoms.resize(6);
    atoms[0]=v1[1];
    atoms[1]=v1[0];
    atoms[2]=axis[0];
    atoms[3]=axis[1];
    atoms[4]=v2[0];
    atoms[5]=v2[1];
    log.printf("  between lines %d-%d and %d-%d, projected on the plane orthogonal to line %d-%d\n",
               v1[0].serial(),v1[1].serial(),v2[0].serial(),v2[1].serial(),axis[0].serial(),axis[1].serial());
  } else {
    parseAtomList(-1,atoms,this);
  }
  unsigned mode=getModeAndSetupValues(this);
  if( mode==1 ) {
    do_cosine=true;
  }

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  checkRead();

  if(pbc) {
    log.printf("  using periodic boundary conditions\n");
  } else {
    log.printf("  without periodic boundary conditions\n");
  }
  requestAtoms(atoms);
}

void Torsion::parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa ) {
  std::vector<AtomNumber> v1,v2,axis;
  aa->parseAtomList("ATOMS",num,t);
  aa->parseAtomList("VECTORA",num,v1);
  aa->parseAtomList("VECTORB",num,v2);
  aa->parseAtomList("AXIS",num,axis);

  if(t.size()==4) {
    if(!(v1.empty() && v2.empty() && axis.empty())) {
      aa->error("ATOMS keyword is not compatible with VECTORA, VECTORB and AXIS keywords");
    }
    aa->log.printf("  between atoms %d %d %d %d\n",t[0].serial(),t[1].serial(),t[2].serial(),t[3].serial());
    t.resize(6);
    t[5]=t[3];
    t[4]=t[2];
    t[3]=t[2];
    t[2]=t[1];
  } else if(t.empty()) {
    if( num>0 && v1.empty() && v2.empty() && axis.empty() ) {
      return;
    }
    if(!(v1.size()==2 && v2.size()==2 && axis.size()==2)) {
      aa->error("VECTORA, VECTORB and AXIS should specify 2 atoms each");
    }
    aa->log.printf("  between lines %d-%d and %d-%d, projected on the plane orthogonal to line %d-%d\n",
                   v1[0].serial(),v1[1].serial(),v2[0].serial(),v2[1].serial(),axis[0].serial(),axis[1].serial());
    t.resize(6);
    t[0]=v1[1];
    t[1]=v1[0];
    t[2]=axis[0];
    t[3]=axis[1];
    t[4]=v2[0];
    t[5]=v2[1];
  } else if( t.size()!=4 ) {
    aa->error("ATOMS should specify 4 atoms");
  }
}

unsigned Torsion::getModeAndSetupValues( ActionWithValue* av ) {
  bool do_cos;
  av->parseFlag("COSINE",do_cos);
  if(do_cos) {
    av->log.printf("  calculating cosine instead of torsion\n");
  }

  av->addValueWithDerivatives();
  if(!do_cos) {
    av->setPeriodic("-pi","pi");
    return 0;
  }
  av->setNotPeriodic();
  return 1;
}

// calculator
void Torsion::calculate() {
  if(pbc) {
    makeWhole();
  }
  if(do_cosine) {
    calculateCV( 1, masses, charges, getPositions(), value, derivs, virial, this );
  } else {
    calculateCV( 0, masses, charges, getPositions(), value, derivs, virial, this );
  }
  for(unsigned i=0; i<6; ++i) {
    setAtomsDerivatives(i,derivs[0][i] );
  }
  setValue(value[0]);
  setBoxDerivatives( virial[0] );
}

void Torsion::calculateCV( const unsigned& mode, const std::vector<double>& masses, const std::vector<double>& charges,
                           const std::vector<Vector>& pos, std::vector<double>& vals, std::vector<std::vector<Vector> >& derivs,
                           std::vector<Tensor>& virial, const ActionAtomistic* aa ) {
  Vector d0=delta(pos[1],pos[0]);
  Vector d1=delta(pos[3],pos[2]);
  Vector d2=delta(pos[5],pos[4]);
  Vector dd0,dd1,dd2;
  PLMD::Torsion t;
  vals[0] = t.compute(d0,d1,d2,dd0,dd1,dd2);
  if(mode==1) {
    dd0 *= -std::sin(vals[0]);
    dd1 *= -std::sin(vals[0]);
    dd2 *= -std::sin(vals[0]);
    vals[0] = std::cos(vals[0]);
  }
  derivs[0][0] = dd0;
  derivs[0][1] = -dd0;
  derivs[0][2] = dd1;
  derivs[0][3] = -dd1;
  derivs[0][4] = dd2;
  derivs[0][5] = -dd2;
  setBoxDerivativesNoPbc( pos, derivs, virial );
}

}
}



