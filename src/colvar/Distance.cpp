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
#include "tools/Pbc.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR DISTANCE
/*
Calculate the distance/s between pairs of atoms.

The following example illustrates how this action can be used to calculate and print the distance between atom 1
and atom 2.

```plumed
d: DISTANCE ATOMS=1,2
PRINT ARG=d FILE=colvar
```

By default the distance is computed in a way that takes periodic
boundary conditions in account. This behavior can be changed by using the NOPBC flag.
Furthermore, if you wish to calculate the vector connecting a pair of atoms you can use the
`COMPONENTS` flag as shown below:

```plumed
d: DISTANCE ATOMS=1,2 COMPONENTS
PRINT ARG=d.x,d.y,d.z FILE=colvar
```

Alternatively, you can calculate the components projected on the lattice vector by using the `SCALED_COMPONENTS`
flag as shown below;

```plumed
d: DISTANCE ATOMS=1,2 SCALED_COMPONENTS
PRINT ARG=d.a,d.b,d.c FILE=colvar
```

The advantage of using `SCALED_COMPONENTS` over `COMPONENTS` is that the a, b and c variables
that are calculated when `SCALED_COMPONENTS` is employed have the proper periodicity. This feature is useful
if you wish to study the motion of a molecule across a membrane.

You can also use this command to calculate multiple indistinguishable distances or vectors with a single
line of PLUMED input.  For example, the following input calculates and outputs the distances between four
pairs of atoms:

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8
PRINT ARG=d FILE=colvar
```

By a similar token, the following input outputs three four dimensional vectors that contain the x, y and z
components of the vectors connecting the four atoms:

```plumed
d: DISTANCE COMPONENTS ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8
PRINT ARG=d.x,d.y,d.z FILE=colvar
```

You can also replace COMPONENTS with SCALED_COMPONENTS in the above input and obtain the projects of these vectors
on the lattice vectors.

## Managing periodic boundary conditions

When using the DISTANCE command to calculate the end-to-end distance for a large polymer you need to ensure that you
are managing PBCs correctly.  This problems that can occur with these calculations are explained at length in the
early parts of the document that is referenced in the bibliography. Notice, however, that the input
provides an example of an input that could be used to compute the end-to-end distance for a polymer
of 100 atoms and keeps it at a value around 5.

```plumed
WHOLEMOLECULES ENTITY0=1-100
e2e: DISTANCE ATOMS=1,100 NOPBC
RESTRAINT ARG=e2e KAPPA=1 AT=5
```

Notice that NOPBC is used here so as to ensure that the distance is calculated correctely even if the end-to-end distance is larger than half the simulation
Also notice that, since many MD codes break molecules across cell boundary, it might be necessary to
use the [WHOLEMOLECULES](WHOLEMOLECULES.md) keyword (also notice that it should be _before_ distance). The list of atoms provided to [WHOLEMOLECULES](WHOLEMOLECULES.md)
here contains all the atoms between 1 and 100. Strictly speaking, this
is not necessary. If you know for sure that atoms with difference in
the index say equal to 10 are _not_ going to be farther than half cell
you can e.g. use

```plumed
WHOLEMOLECULES ENTITY0=1,10,20,30,40,50,60,70,80,90,100
e2e: DISTANCE ATOMS=1,100 NOPBC
RESTRAINT ARG=e2e KAPPA=1 AT=5
```

Just be sure that the ordered list provide to [WHOLEMOLECULES](WHOLEMOLECULES.md) has the following
properties:

- Consecutive atoms should be closer than half-cell throughout the entire simulation.
- Atoms required later for the distance (e.g. 1 and 100) should be included in the list

The following example shows how to take periodicity into account when computing the z-component of a distance

```plumed
# this is a center of mass of a large group
c: COM ATOMS=1-100
# this is the distance between atom 101 and the group
d: DISTANCE ATOMS=c,101 COMPONENTS
# this makes a new variable, dd, equal to d and periodic, with domain -10,10
# this is the right choise if e.g. the cell is orthorombic and its size in
# z direction is 20.
dz: COMBINE ARG=d.z PERIODIC=-10,10
# metadynamics on dd
METAD ARG=dz SIGMA=0.1 HEIGHT=0.1 PACE=200
```

You can use the same input even if the DISTANCE command is calculating the vectors connecting multiple pairs of atoms.
However, using SCALED_COMPONENTS ensures this problem does not arise because these variables are always periodic
with domain (-0.5,+0.5).

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR DISTANCE_SCALAR
/*
Calculate the distance between a pair of atoms

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR DISTANCE_VECTOR
/*
Calculate a vector containing the distances between various pairs of atoms

\par Examples

*/
//+ENDPLUMEDOC

class Distance : public Colvar {
  bool components;
  bool scaled_components;
  bool pbc;

  std::vector<double> value, masses, charges;
  std::vector<std::vector<Vector> > derivs;
  std::vector<Tensor> virial;
public:
  static void registerKeywords( Keywords& keys );
  explicit Distance(const ActionOptions&);
  static void parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa );
  static unsigned getModeAndSetupValues( ActionWithValue* av );
// active methods:
  void calculate() override;
  static void calculateCV( const unsigned& mode, const std::vector<double>& masses, const std::vector<double>& charges,
                           const std::vector<Vector>& pos, std::vector<double>& vals, std::vector<std::vector<Vector> >& derivs,
                           std::vector<Tensor>& virial, const ActionAtomistic* aa );
};

typedef ColvarShortcut<Distance> DistanceShortcut;
PLUMED_REGISTER_ACTION(DistanceShortcut,"DISTANCE")
PLUMED_REGISTER_ACTION(Distance,"DISTANCE_SCALAR")
typedef MultiColvarTemplate<Distance> DistanceMulti;
PLUMED_REGISTER_ACTION(DistanceMulti,"DISTANCE_VECTOR")

void Distance::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.setDisplayName("DISTANCE");
  constexpr auto scalarOrVector = Keywords::componentType::scalar | Keywords::componentType::vector;
  keys.add("atoms","ATOMS","the pair of atom that we are calculating the distance between");
  keys.addFlag("COMPONENTS",false,"calculate the x, y and z components of the distance separately and store them as label.x, label.y and label.z");
  keys.addFlag("SCALED_COMPONENTS",false,"calculate the a, b and c scaled components of the distance separately and store them as label.a, label.b and label.c");
  keys.addOutputComponent("x","COMPONENTS",scalarOrVector,"the x-component of the vector connecting the two atoms");
  keys.addOutputComponent("y","COMPONENTS",scalarOrVector,"the y-component of the vector connecting the two atoms");
  keys.addOutputComponent("z","COMPONENTS",scalarOrVector,"the z-component of the vector connecting the two atoms");
  keys.addOutputComponent("a","SCALED_COMPONENTS",scalarOrVector,"the normalized projection on the first lattice vector of the vector connecting the two atoms");
  keys.addOutputComponent("b","SCALED_COMPONENTS",scalarOrVector,"the normalized projection on the second lattice vector of the vector connecting the two atoms");
  keys.addOutputComponent("c","SCALED_COMPONENTS",scalarOrVector,"the normalized projection on the third lattice vector of the vector connecting the two atoms");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  keys.setValueDescription(scalarOrVector,"the DISTANCE between this pair of atoms");
  keys.addDOI("10.1007/978-1-4939-9608-7_21");
}

Distance::Distance(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  components(false),
  scaled_components(false),
  pbc(true),
  value(1),
  derivs(1),
  virial(1) {
  derivs[0].resize(2);
  std::vector<AtomNumber> atoms;
  parseAtomList(-1,atoms,this);
  if(atoms.size()!=2) {
    error("Number of specified atoms should be 2");
  }

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  if(pbc) {
    log.printf("  using periodic boundary conditions\n");
  } else {
    log.printf("  without periodic boundary conditions\n");
  }

  unsigned mode = getModeAndSetupValues( this );
  if(mode==1) {
    components=true;
  } else if(mode==2) {
    scaled_components=true;
  }
  if( components || scaled_components ) {
    value.resize(3);
    derivs.resize(3);
    virial.resize(3);
    for(unsigned i=0; i<3; ++i) {
      derivs[i].resize(2);
    }
  }
  requestAtoms(atoms);
}

void Distance::parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa ) {
  aa->parseAtomList("ATOMS",num,t);
  if( t.size()==2 ) {
    aa->log.printf("  between atoms %d %d\n",t[0].serial(),t[1].serial());
  }
}

unsigned Distance::getModeAndSetupValues( ActionWithValue* av ) {
  bool c;
  av->parseFlag("COMPONENTS",c);
  bool sc;
  av->parseFlag("SCALED_COMPONENTS",sc);
  if( c && sc ) {
    av->error("COMPONENTS and SCALED_COMPONENTS are not compatible");
  }

  if(c) {
    av->addComponentWithDerivatives("x");
    av->componentIsNotPeriodic("x");
    av->addComponentWithDerivatives("y");
    av->componentIsNotPeriodic("y");
    av->addComponentWithDerivatives("z");
    av->componentIsNotPeriodic("z");
    av->log<<"  WARNING: components will not have the proper periodicity - see manual\n";
    return 1;
  } else if(sc) {
    av->addComponentWithDerivatives("a");
    av->componentIsPeriodic("a","-0.5","+0.5");
    av->addComponentWithDerivatives("b");
    av->componentIsPeriodic("b","-0.5","+0.5");
    av->addComponentWithDerivatives("c");
    av->componentIsPeriodic("c","-0.5","+0.5");
    return 2;
  }
  av->addValueWithDerivatives();
  av->setNotPeriodic();
  return 0;
}

// calculator
void Distance::calculate() {

  if(pbc) {
    makeWhole();
  }

  if( components ) {
    calculateCV( 1, masses, charges, getPositions(), value, derivs, virial, this );
    Value* valuex=getPntrToComponent("x");
    Value* valuey=getPntrToComponent("y");
    Value* valuez=getPntrToComponent("z");

    for(unsigned i=0; i<2; ++i) {
      setAtomsDerivatives(valuex,i,derivs[0][i] );
    }
    setBoxDerivatives(valuex,virial[0]);
    valuex->set(value[0]);

    for(unsigned i=0; i<2; ++i) {
      setAtomsDerivatives(valuey,i,derivs[1][i] );
    }
    setBoxDerivatives(valuey,virial[1]);
    valuey->set(value[1]);

    for(unsigned i=0; i<2; ++i) {
      setAtomsDerivatives(valuez,i,derivs[2][i] );
    }
    setBoxDerivatives(valuez,virial[2]);
    valuez->set(value[2]);
  } else if( scaled_components ) {
    calculateCV( 2, masses, charges, getPositions(), value, derivs, virial, this );

    Value* valuea=getPntrToComponent("a");
    Value* valueb=getPntrToComponent("b");
    Value* valuec=getPntrToComponent("c");
    for(unsigned i=0; i<2; ++i) {
      setAtomsDerivatives(valuea,i,derivs[0][i] );
    }
    valuea->set(value[0]);
    for(unsigned i=0; i<2; ++i) {
      setAtomsDerivatives(valueb,i,derivs[1][i] );
    }
    valueb->set(value[1]);
    for(unsigned i=0; i<2; ++i) {
      setAtomsDerivatives(valuec,i,derivs[2][i] );
    }
    valuec->set(value[2]);
  } else  {
    calculateCV( 0, masses, charges, getPositions(), value, derivs, virial, this );
    for(unsigned i=0; i<2; ++i) {
      setAtomsDerivatives(i,derivs[0][i] );
    }
    setBoxDerivatives(virial[0]);
    setValue           (value[0]);
  }
}

void Distance::calculateCV( const unsigned& mode, const std::vector<double>& masses, const std::vector<double>& charges,
                            const std::vector<Vector>& pos, std::vector<double>& vals, std::vector<std::vector<Vector> >& derivs,
                            std::vector<Tensor>& virial, const ActionAtomistic* aa ) {
  Vector distance=delta(pos[0],pos[1]);
  const double value=distance.modulo();
  const double invvalue=1.0/value;

  if(mode==1) {
    derivs[0][0] = Vector(-1,0,0);
    derivs[0][1] = Vector(+1,0,0);
    vals[0] = distance[0];

    derivs[1][0] = Vector(0,-1,0);
    derivs[1][1] = Vector(0,+1,0);
    vals[1] = distance[1];

    derivs[2][0] = Vector(0,0,-1);
    derivs[2][1] = Vector(0,0,+1);
    vals[2] = distance[2];
    setBoxDerivativesNoPbc( pos, derivs, virial );
  } else if(mode==2) {
    Vector d=aa->getPbc().realToScaled(distance);
    derivs[0][0] = matmul(aa->getPbc().getInvBox(),Vector(-1,0,0));
    derivs[0][1] = matmul(aa->getPbc().getInvBox(),Vector(+1,0,0));
    vals[0] = Tools::pbc(d[0]);
    derivs[1][0] = matmul(aa->getPbc().getInvBox(),Vector(0,-1,0));
    derivs[1][1] = matmul(aa->getPbc().getInvBox(),Vector(0,+1,0));
    vals[1] = Tools::pbc(d[1]);
    derivs[2][0] = matmul(aa->getPbc().getInvBox(),Vector(0,0,-1));
    derivs[2][1] = matmul(aa->getPbc().getInvBox(),Vector(0,0,+1));
    vals[2] = Tools::pbc(d[2]);
  } else {
    derivs[0][0] = -invvalue*distance;
    derivs[0][1] = invvalue*distance;
    setBoxDerivativesNoPbc( pos, derivs, virial );
    vals[0] = value;
  }
}

}
}



