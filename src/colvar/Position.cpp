/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2023 The plumed team
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
#include "tools/Pbc.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR POSITION
/*
Calculate the components of the position of an atom or atoms.

To print the position of atom one to a file you can use an input like this:

```plumed
p: POSITION ATOM=1
PRINT ARG=p.x,p.y,p.z FILE=colvar
```

To print the position of four atoms you can use an input like this:

```plumed
p: POSITION ATOMS=1-4
PRINT ARG=p.x,p.y,p.z FILE=colvar
```

The three output values, p.x, p.y and p.z, here are all four dimensional vectors.  Furthermore, if you wish
to use a procedure akin to that described in the documentation for [WHOLEMOLECULES](WHOLEMOLECULES.md) to ensure
that molecules are reassembled if they are broken by periodic boundary conditions you can use the `WHOLEMOLECULES`
flag as shown below:

```plumed
p: POSITION ATOMS=1-4 WHOLEMOLECULES
PRINT ARG=p.x,p.y,p.z FILE=colvar
```

This is used in the [CENTER](CENTER.md) when we compute centers with arbitrary weights for shortcuts.

## Periodic boundary conditions

!!! warning ""

    Notice that single components will not have the proper periodicity!

If you need the values to be consistent through PBC you can use SCALED_COMPONENTS flag as shown below:

```plumed
p: POSITION ATOM=1 SCALED_COMPONENTS
PRINT ARG=p.a,p.b,p.c FILE=colvar
```

This flag ensures that values are output that are, by construction, in the -0.5,0.5 domain. The output is
similar to the equivalent flag for [DISTANCE](DISTANCE.md).

If it also important to note that by default the positions are output by calculating minimal image distance between
the instantaneous position of the atoms and the position of the origin, $(0.0,0.0,0.0)$.
This behaviour can be changed by using NOPBC as shown below, which will output the position of the atom directly.

```plumed
p: POSITION ATOM=1 NOPBC
PRINT ARG=p.x,p.y,p.z FILE=colvar
```

!!! warning ""

    This variable should be used with extreme care since it allows you to easily get in troubles.
    It can be only be used if the Hamiltonian is not invariant for translation (i.e. there are other absolute positions which are biased, e.g. by position restraints)
    and cell size and shapes are fixed through the simulation.

If you are not in this situation and still want to use the absolute position of an atom you should first fix the reference frame
by using [FIT_TO_TEMPLATE](FIT_TO_TEMPLATE.md) as shown in the example below

```plumed
#SETTINGS INPUTFILES=regtest/basic/rt63/align.pdb
# align to a template
FIT_TO_TEMPLATE REFERENCE=regtest/basic/rt63/align.pdb
p: POSITION ATOM=3
PRINT ARG=p.x,p.y,p.z
```

*/
//+ENDPLUMEDOC

class Position : public Colvar {
  bool scaled_components;
  bool pbc;
  std::vector<double> value;
  std::vector<double> derivs;
public:
  static void registerKeywords( Keywords& keys );
  explicit Position(const ActionOptions&);
  static void parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa );
  static unsigned getModeAndSetupValues( ActionWithValue* av );
// active methods:
  void calculate() override;
  static void calculateCV( const ColvarInput& cvin, ColvarOutput& cvout );
};

typedef ColvarShortcut<Position> PositionShortcut;
PLUMED_REGISTER_ACTION(PositionShortcut,"POSITION")
PLUMED_REGISTER_ACTION(Position,"POSITION_SCALAR")
typedef MultiColvarTemplate<Position> PositionMulti;
PLUMED_REGISTER_ACTION(PositionMulti,"POSITION_VECTOR")

void Position::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.setDisplayName("POSITION");
  keys.add("atoms","ATOM","the atom number");
  keys.add("atoms","ATOMS","the atom numbers that you would like to use the positions of");
  keys.addFlag("WHOLEMOLECULES",false,"if this is a vector of positions do you want to make the positions into a whole before");
  keys.addFlag("SCALED_COMPONENTS",false,"calculate the a, b and c scaled components of the position separately and store them as label.a, label.b and label.c");
  keys.addOutputComponent("x","default","scalar/vector","the x-component of the atom position");
  keys.addOutputComponent("y","default","scalar/vector","the y-component of the atom position");
  keys.addOutputComponent("z","default","scalar/vector","the z-component of the atom position");
  keys.addOutputComponent("a","SCALED_COMPONENTS","scalar/vector","the normalized projection on the first lattice vector of the atom position");
  keys.addOutputComponent("b","SCALED_COMPONENTS","scalar/vector","the normalized projection on the second lattice vector of the atom position");
  keys.addOutputComponent("c","SCALED_COMPONENTS","scalar/vector","the normalized projection on the third lattice vector of the atom position");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  keys.reset_style("NUMERICAL_DERIVATIVES","hidden");
}

Position::Position(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  scaled_components(false),
  pbc(true),
  value(3) {
  std::vector<AtomNumber> atoms;
  parseAtomList(-1,atoms,this);
  unsigned mode=getModeAndSetupValues(this);
  if( mode==1 ) {
    scaled_components=true;
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

void Position::parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa ) {
  aa->parseAtomList("ATOM",num,t);
  if( t.size()==1 ) {
    aa->log.printf("  for atom %d\n",t[0].serial());
  } else if( num<0 || t.size()!=0 ) {
    aa->error("Number of specified atoms should be 1");
  }
}

unsigned Position::getModeAndSetupValues( ActionWithValue* av ) {
  bool sc;
  av->parseFlag("SCALED_COMPONENTS",sc);
  if(sc) {
    av->addComponentWithDerivatives("a");
    av->componentIsPeriodic("a","-0.5","+0.5");
    av->addComponentWithDerivatives("b");
    av->componentIsPeriodic("b","-0.5","+0.5");
    av->addComponentWithDerivatives("c");
    av->componentIsPeriodic("c","-0.5","+0.5");
    return 1;
  }
  av->addComponentWithDerivatives("x");
  av->componentIsNotPeriodic("x");
  av->addComponentWithDerivatives("y");
  av->componentIsNotPeriodic("y");
  av->addComponentWithDerivatives("z");
  av->componentIsNotPeriodic("z");
  av->log<<"  WARNING: components will not have the proper periodicity - see manual\n";
  return 0;
}

// calculator
void Position::calculate() {

  std::vector<Vector> distance(1);
  if(pbc) {
    distance[0]=pbcDistance(Vector(0.0,0.0,0.0),getPosition(0));
  } else {
    distance[0]=delta(Vector(0.0,0.0,0.0),getPosition(0));
  }

  ColvarOutput cvout = ColvarOutput::createColvarOutput(value,derivs,this);
  if(scaled_components) {
    calculateCV( ColvarInput::createColvarInput( 1, distance, this ), cvout );
    Value* valuea=getPntrToComponent("a");
    Value* valueb=getPntrToComponent("b");
    Value* valuec=getPntrToComponent("c");
    setAtomsDerivatives (valuea,0,cvout.getAtomDerivatives(0,0));
    valuea->set(value[0]);
    setAtomsDerivatives (valueb,0,cvout.getAtomDerivatives(1,0));
    valueb->set(value[1]);
    setAtomsDerivatives (valuec,0,cvout.getAtomDerivatives(2,0));
    valuec->set(value[2]);
  } else {
    calculateCV( ColvarInput::createColvarInput( 0, distance, this ), cvout );
    Value* valuex=getPntrToComponent("x");
    Value* valuey=getPntrToComponent("y");
    Value* valuez=getPntrToComponent("z");

    setAtomsDerivatives (valuex,0,cvout.getAtomDerivatives(0,0));
    setBoxDerivatives   (valuex,cvout.virial[0]);
    valuex->set(value[0]);

    setAtomsDerivatives (valuey,0,cvout.getAtomDerivatives(1,0));
    setBoxDerivatives   (valuey,cvout.virial[1]);
    valuey->set(value[1]);

    setAtomsDerivatives (valuez,0,cvout.getAtomDerivatives(2,0));
    setBoxDerivatives   (valuez,cvout.virial[2]);
    valuez->set(value[2]);
  }
}

void Position::calculateCV( const ColvarInput& cvin, ColvarOutput& cvout ) {
  if( cvin.mode==1 ) {
    Vector d=cvin.pbc.realToScaled(Vector(cvin.pos[0][0],cvin.pos[0][1],cvin.pos[0][2]));
    cvout.values[0]=Tools::pbc(d[0]);
    cvout.values[1]=Tools::pbc(d[1]);
    cvout.values[2]=Tools::pbc(d[2]);
    cvout.derivs[0][0]=matmul(cvin.pbc.getInvBox(),Vector(+1,0,0));
    cvout.derivs[1][0]=matmul(cvin.pbc.getInvBox(),Vector(0,+1,0));
    cvout.derivs[2][0]=matmul(cvin.pbc.getInvBox(),Vector(0,0,+1));
  } else {
    for(unsigned i=0; i<3; ++i) {
      cvout.values[i]=cvin.pos[0][i];
    }
    cvout.derivs[0][0]=Vector(+1,0,0);
    cvout.derivs[1][0]=Vector(0,+1,0);
    cvout.derivs[2][0]=Vector(0,0,+1);
    cvout.virial.set(0, Tensor(Vector(cvin.pos[0][0],cvin.pos[0][1],cvin.pos[0][2]),Vector(-1,0,0)) );
    cvout.virial.set(1, Tensor(Vector(cvin.pos[0][0],cvin.pos[0][1],cvin.pos[0][2]),Vector(0,-1,0)) );
    cvout.virial.set(2, Tensor(Vector(cvin.pos[0][0],cvin.pos[0][1],cvin.pos[0][2]),Vector(0,0,-1)) );
  }
}

}
}



