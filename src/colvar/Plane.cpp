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
#include "core/ActionRegister.h"
#include "MultiColvarTemplate.h"
#include "tools/Pbc.h"

#include <string>
#include <cmath>

//+PLUMEDOC COLVAR PLANE
/*
Calculate the plane perpendicular to two vectors in order to represent the orientation of a planar molecule.

\par Examples


*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR PLANE_SCALAR
/*
Calculate the plane perpendicular to two vectors in order to represent the orientation of a planar molecule.

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR PLANE_VECTOR
/*
Calculate the plane perpendicular to two vectors in order to represent the orientation of a planar molecule multiple times.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace colvar {

class Plane : public Colvar {
private:
  bool pbc;
  std::vector<double> value, masses, charges;
  std::vector<std::vector<Vector> > derivs;
  std::vector<Tensor> virial;
public:
  static void registerKeywords( Keywords& keys );
  explicit Plane(const ActionOptions&);
  static void parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa );
  static unsigned getModeAndSetupValues( ActionWithValue* av );
// active methods:
  void calculate() override;
  static void calculateCV( const unsigned& mode, const std::vector<double>& masses, const std::vector<double>& charges,
                           const std::vector<Vector>& pos, std::vector<double>& vals, std::vector<std::vector<Vector> >& derivs,
                           std::vector<Tensor>& virial, const ActionAtomistic* aa );
};

typedef ColvarShortcut<Plane> PlaneShortcut;
PLUMED_REGISTER_ACTION(PlaneShortcut,"PLANE")
PLUMED_REGISTER_ACTION(Plane,"PLANE_SCALAR")
typedef MultiColvarTemplate<Plane> PlaneMulti;
PLUMED_REGISTER_ACTION(PlaneMulti,"PLANE_VECTOR")

void Plane::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys ); keys.setDisplayName("PLANE");
  keys.add("atoms","ATOMS","the three or four atoms whose plane we are computing");
  keys.addOutputComponent("x","default","the x-component of the vector that is normal to the plane containing the atoms");
  keys.addOutputComponent("y","default","the y-component of the vector that is normal to the plane containing the atoms");
  keys.addOutputComponent("z","default","the z-component of the vector that is normal to the plane containing the atoms");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
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
  } else if( num<0 || atoms.size()>0 ) aa->error("Number of specified atoms should be either 3 or 4");
}

unsigned Plane::getModeAndSetupValues( ActionWithValue* av ) {
  av->addComponentWithDerivatives("x"); av->componentIsNotPeriodic("x");
  av->addComponentWithDerivatives("y"); av->componentIsNotPeriodic("y");
  av->addComponentWithDerivatives("z"); av->componentIsNotPeriodic("z");
  return 0;
}

Plane::Plane(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  value(3),
  derivs(3),
  virial(3)
{
  for(unsigned i=0; i<3; ++i) derivs[i].resize(4);
  std::vector<AtomNumber> atoms; parseAtomList(-1,atoms,this);
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");

  unsigned mode = getModeAndSetupValues( this );
  requestAtoms(atoms);
  checkRead();
}

void Plane::calculate() {

  if(pbc) makeWhole();
  calculateCV( 0, masses, charges, getPositions(), value, derivs, virial, this );
  setValue( value[0] );
  for(unsigned i=0; i<derivs[0].size(); ++i) setAtomsDerivatives( i, derivs[0][i] );
  setBoxDerivatives( virial[0] );
}

void Plane::calculateCV( const unsigned& mode, const std::vector<double>& masses, const std::vector<double>& charges,
                         const std::vector<Vector>& pos, std::vector<double>& vals, std::vector<std::vector<Vector> >& derivs,
                         std::vector<Tensor>& virial, const ActionAtomistic* aa ) {
  Vector d1=delta( pos[1], pos[0] );
  Vector d2=delta( pos[2], pos[3] );
  Vector cp = crossProduct( d1, d2 );

  derivs[0][0] = crossProduct( Vector(-1.0,0,0), d2 );
  derivs[0][1] = crossProduct( Vector(+1.0,0,0), d2 );
  derivs[0][2] = crossProduct( Vector(-1.0,0,0), d1 );
  derivs[0][3] = crossProduct( Vector(+1.0,0,0), d1 );
  virial[0] = Tensor(d1,crossProduct(Vector(+1.0,0,0), d2)) + Tensor( d2, crossProduct(Vector(-1.0,0,0), d1));
  vals[0] = cp[0];

  derivs[1][0] = crossProduct( Vector(0,-1.0,0), d2 );
  derivs[1][1] = crossProduct( Vector(0,+1.0,0), d2 );
  derivs[1][2] = crossProduct( Vector(0,-1.0,0), d1 );
  derivs[1][3] = crossProduct( Vector(0,+1.0,0), d1 );
  virial[1] = Tensor(d1,crossProduct(Vector(0,+1.0,0), d2)) + Tensor( d2, crossProduct(Vector(0,-1.0,0), d1));
  vals[1] = cp[1];

  derivs[2][0] = crossProduct( Vector(0,0,-1.0), d2 );
  derivs[2][1] = crossProduct( Vector(0,0,+1.0), d2 );
  derivs[2][2] = crossProduct( Vector(0,0,-1.0), d1 );
  derivs[2][3] = crossProduct( Vector(0,0,+1.0), d1 );
  virial[2] = Tensor(d1,crossProduct(Vector(0,0,+1.0), d2)) + Tensor( d2, crossProduct(Vector(0,0,-1.0), d1));
  vals[2] = cp[2];
}

}
}



