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
#include "core/ActionRegister.h"
#include "MultiColvarTemplate.h"
#include "tools/Pbc.h"

#include <string>
#include <cmath>

//+PLUMEDOC COLVAR PLANE
/*
Calculate the plane perpendicular to two vectors in order to represent the orientation of a planar molecule.

To calculate the orientation of the plane connecting atoms 1, 2 and 3 you use an input like this:

```plumed
p: PLANE ATOMS=1,2,3
PRINT ARG=p.x,p.y,p.z FILE=colvar
```

The three components, p.x, p.y and p.z, output by the PLANE action here are the x, y and z components of the normal
vector to the plane that is obtained by taking the cross product between the vector connecting atoms 1 and 2 and
the vector connecting atoms 2 and 3.  Notice that the default here is to evaluate these two vectors in a way that takes
any periodic boundary conditions (PBC) into account. If you wish to disregard the PBC you can use the NOPBC flag as shown in the following input:

```plumed
p: PLANE ATOMS=1,2,3 NOPBC
PRINT ARG=p.x,p.y,p.z FILE=colvar
```

To calculate the cross product of the vector connecting atoms 1 and 2 and the vector connecting atoms 3 and 4 you use
an input like this:

```plumed
p: PLANE ATOMS=1,2,3,4
PRINT ARG=p.x,p.y,p.z FILE=colvar
```

If you have multiple molecules and wish to determine the orientations of the planes containing all them with one line of PLUMED
input you can use an input like this:

```plumed
p: PLANE ATOMS1=1,2,3 ATOMS2=4,5,6 ATOMS3=7,8,9 ATOMS4=10,11,12
PRINT ARG=p.x,p.y,p.z FILE=colvar
```

The output from this command consists of 3 vectors with 4 components. These vectors, p.x, p.y and p.z, contain the x, y and z components
of the normals to the planes of the molecules.  Commands similar to this are useful for variables that can be used to monitor
nucleation of molecular crystals such as [SMAC](SMAC.md).

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace colvar {

class Plane : public Colvar {
private:
  bool pbc;
  std::vector<double> value;
  std::vector<double> derivs;
public:
  static void registerKeywords( Keywords& keys );
  explicit Plane(const ActionOptions&);
  static void parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa );
  static unsigned getModeAndSetupValues( ActionWithValue* av );
// active methods:
  void calculate() override;
  static void calculateCV( const ColvarInput& cvin, ColvarOutput& cvout );
};

typedef ColvarShortcut<Plane> PlaneShortcut;
PLUMED_REGISTER_ACTION(PlaneShortcut,"PLANE")
PLUMED_REGISTER_ACTION(Plane,"PLANE_SCALAR")
typedef MultiColvarTemplate<Plane> PlaneMulti;
PLUMED_REGISTER_ACTION(PlaneMulti,"PLANE_VECTOR")

void Plane::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.setDisplayName("PLANE");
  keys.add("atoms","ATOMS","the three or four atoms whose plane we are computing");
  keys.addOutputComponent("x","default","scalar/vector","the x-component of the vector that is normal to the plane containing the atoms");
  keys.addOutputComponent("y","default","scalar/vector","the y-component of the vector that is normal to the plane containing the atoms");
  keys.addOutputComponent("z","default","scalar/vector","the z-component of the vector that is normal to the plane containing the atoms");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  keys.reset_style("NUMERICAL_DERIVATIVES","hidden");
}

void Plane::parseAtomList( const int& num, std::vector<AtomNumber>& atoms, ActionAtomistic* aa ) {
  aa->parseAtomList("ATOMS",num,atoms);
  if(atoms.size()==3) {
    aa->log.printf("  containing atoms %d %d %d\n",atoms[0].serial(),atoms[1].serial(),atoms[2].serial());
    atoms.resize(4);
    atoms[3]=atoms[2];
    atoms[2]=atoms[1];
  } else if(atoms.size()==4) {
    aa->log.printf("  containing lines %d-%d and %d-%d\n",atoms[0].serial(),atoms[1].serial(),atoms[2].serial(),atoms[3].serial());
  } else if( num<0 || atoms.size()>0 ) {
    aa->error("Number of specified atoms should be either 3 or 4");
  }
}

unsigned Plane::getModeAndSetupValues( ActionWithValue* av ) {
  av->addComponentWithDerivatives("x");
  av->componentIsNotPeriodic("x");
  av->addComponentWithDerivatives("y");
  av->componentIsNotPeriodic("y");
  av->addComponentWithDerivatives("z");
  av->componentIsNotPeriodic("z");
  return 0;
}

Plane::Plane(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  value(3) {
  std::vector<AtomNumber> atoms;
  parseAtomList(-1,atoms,this);
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  if(pbc) {
    log.printf("  using periodic boundary conditions\n");
  } else {
    log.printf("  without periodic boundary conditions\n");
  }

  getModeAndSetupValues( this );
  requestAtoms(atoms);
  checkRead();
}

void Plane::calculate() {

  if(pbc) {
    makeWhole();
  }
  ColvarOutput cvout = ColvarOutput::createColvarOutput(value,derivs,this);
  calculateCV( ColvarInput::createColvarInput( 0, getPositions(), this ), cvout );
  Value* valuex=getPntrToComponent("x");
  Value* valuey=getPntrToComponent("y");
  Value* valuez=getPntrToComponent("z");

  for(unsigned i=0; i<getPositions().size(); ++i) {
    setAtomsDerivatives( valuex, i, cvout.getAtomDerivatives(0,i) );
  }
  setBoxDerivatives( valuex, cvout.virial[0] );
  valuex->set( value[0] );

  for(unsigned i=0; i<getPositions().size(); ++i) {
    setAtomsDerivatives( valuey, i, cvout.getAtomDerivatives(1,i) );
  }
  setBoxDerivatives( valuey, cvout.virial[1] );
  valuey->set( value[1] );

  for(unsigned i=0; i<getPositions().size(); ++i) {
    setAtomsDerivatives( valuez, i, cvout.getAtomDerivatives(2,i) );
  }
  setBoxDerivatives( valuez, cvout.virial[2] );
  valuez->set( value[2] );
}

void Plane::calculateCV( const ColvarInput& cvin, ColvarOutput& cvout ) {
  Vector d1=delta( cvin.pos[1], cvin.pos[0] );
  Vector d2=delta( cvin.pos[2], cvin.pos[3] );
  Vector cp = crossProduct( d1, d2 );

  cvout.derivs[0][0] = crossProduct( Vector(-1.0,0,0), d2 );
  cvout.derivs[0][1] = crossProduct( Vector(+1.0,0,0), d2 );
  cvout.derivs[0][2] = crossProduct( Vector(-1.0,0,0), d1 );
  cvout.derivs[0][3] = crossProduct( Vector(+1.0,0,0), d1 );
  cvout.virial.set( 0, Tensor(d1,crossProduct(Vector(+1.0,0,0), d2)) + Tensor( d2, crossProduct(Vector(-1.0,0,0), d1)) );
  cvout.values[0] = cp[0];

  cvout.derivs[1][0] = crossProduct( Vector(0,-1.0,0), d2 );
  cvout.derivs[1][1] = crossProduct( Vector(0,+1.0,0), d2 );
  cvout.derivs[1][2] = crossProduct( Vector(0,-1.0,0), d1 );
  cvout.derivs[1][3] = crossProduct( Vector(0,+1.0,0), d1 );
  cvout.virial.set(1, Tensor(d1,crossProduct(Vector(0,+1.0,0), d2)) + Tensor( d2, crossProduct(Vector(0,-1.0,0), d1)) );
  cvout.values[1] = cp[1];

  cvout.derivs[2][0] = crossProduct( Vector(0,0,-1.0), d2 );
  cvout.derivs[2][1] = crossProduct( Vector(0,0,+1.0), d2 );
  cvout.derivs[2][2] = crossProduct( Vector(0,0,-1.0), d1 );
  cvout.derivs[2][3] = crossProduct( Vector(0,0,+1.0), d1 );
  cvout.virial.set(2, Tensor(d1,crossProduct(Vector(0,0,+1.0), d2)) + Tensor( d2, crossProduct(Vector(0,0,-1.0), d1)) );
  cvout.values[2] = cp[2];
}

}
}



