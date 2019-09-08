/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#include "MultiColvarBase.h"
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "tools/Pbc.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace multicolvar {

class Plane : public MultiColvarBase {
private:
public:
  static void registerKeywords( Keywords& keys );
  explicit Plane(const ActionOptions&);
// active methods:
  void compute( const std::vector<Vector>& pos, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(Plane,"PLANE")

void Plane::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys );
  keys.addOutputComponent("x","default","the x-component of the vector that is normal to the plane containing the atoms");
  keys.addOutputComponent("y","default","the y-component of the vector that is normal to the plane containing the atoms");
  keys.addOutputComponent("z","default","the z-component of the vector that is normal to the plane containing the atoms");
}

Plane::Plane(const ActionOptions&ao):
  Action(ao),
  MultiColvarBase(ao)
{
  if(getNumberOfAtomsInEachCV()==3 ) useFourAtomsForEachCV();
  if( getNumberOfAtomsInEachCV()!=4 ) error("Number of specified atoms should be 3 or 4");
  checkRead();
  addComponentWithDerivatives("x"); componentIsNotPeriodic("x");
  addComponentWithDerivatives("y"); componentIsNotPeriodic("y");
  addComponentWithDerivatives("z"); componentIsNotPeriodic("z");
}

// calculator
void Plane::compute( const std::vector<Vector>& pos, MultiValue& myvals  ) const {

  Vector d1=delta( pos[1], pos[0] );
  Vector d2=delta( pos[2], pos[3] );
  Vector cp = crossProduct( d1, d2 );

  addAtomsDerivatives( 0, 0, crossProduct( Vector(-1.0,0,0), d2 ), myvals);
  addAtomsDerivatives( 0, 1, crossProduct( Vector(+1.0,0,0), d2 ), myvals);
  addAtomsDerivatives( 0, 2, crossProduct( Vector(-1.0,0,0), d1 ), myvals);
  addAtomsDerivatives( 0, 3, crossProduct( Vector(+1.0,0,0), d1 ), myvals);
  addBoxDerivatives( 0, Tensor(d1,crossProduct(Vector(+1.0,0,0), d2)) + Tensor( d2, crossProduct(Vector(-1.0,0,0), d1)), myvals );
  setValue( 0, cp[0], myvals );

  addAtomsDerivatives( 1, 0, crossProduct( Vector(0,-1.0,0), d2 ), myvals);
  addAtomsDerivatives( 1, 1, crossProduct( Vector(0,+1.0,0), d2 ), myvals);
  addAtomsDerivatives( 1, 2, crossProduct( Vector(0,-1.0,0), d1 ), myvals);
  addAtomsDerivatives( 1, 3, crossProduct( Vector(0,+1.0,0), d1 ), myvals);
  addBoxDerivatives( 1, Tensor(d1,crossProduct(Vector(0,+1.0,0), d2)) + Tensor( d2, crossProduct(Vector(0,-1.0,0), d1)), myvals );
  setValue( 1, cp[1], myvals );

  addAtomsDerivatives( 2, 0, crossProduct( Vector(0,0,-1.0), d2 ), myvals);
  addAtomsDerivatives( 2, 1, crossProduct( Vector(0,0,+1.0), d2 ), myvals);
  addAtomsDerivatives( 2, 2, crossProduct( Vector(0,0,-1.0), d1 ), myvals);
  addAtomsDerivatives( 2, 3, crossProduct( Vector(0,0,+1.0), d1 ), myvals);
  addBoxDerivatives( 2, Tensor(d1,crossProduct(Vector(0,0,+1.0), d2)) + Tensor( d2, crossProduct(Vector(0,0,-1.0), d1)), myvals );
  setValue( 2, cp[2], myvals );
}

class PlaneShortcut : public ActionShortcut {
private:
  void createVectorNormInput( const std::string& ilab, const std::string& olab, const std::string& vlab );
public:
  static void registerKeywords( Keywords& keys );
  PlaneShortcut(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(PlaneShortcut,"PLANES")

void PlaneShortcut::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys );
  keys.addFlag("VMEAN",false,"calculate the norm of the mean vector.");
  keys.addOutputComponent("_vmean","VMEAN","the norm of the mean vector");
  keys.addFlag("VSUM",false,"calculate the norm of the sum of all the vectors");
  keys.addOutputComponent("_vsum","VSUM","the norm of the mean vector");
}

PlaneShortcut::PlaneShortcut(const ActionOptions&ao):
Action(ao),
ActionShortcut(ao)
{
  bool vmean, vsum; parseFlag("VMEAN",vmean); parseFlag("VSUM",vsum);
  readInputLine( getShortcutLabel() + ": PLANE " + convertInputLineToString() );
  if( vmean ) {
    readInputLine( getShortcutLabel() + "_xs: COMBINE ARG=" + getShortcutLabel() + ".x NORMALIZE PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_ys: COMBINE ARG=" + getShortcutLabel() + ".y NORMALIZE PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_zs: COMBINE ARG=" + getShortcutLabel() + ".z NORMALIZE PERIODIC=NO");
    // Now calculate the total length of the vector
    createVectorNormInput( getShortcutLabel(), getShortcutLabel() + "_vmean", "s" );
  }
  if( vsum ) {
    readInputLine( getShortcutLabel() + "_xz: COMBINE ARG=" + getShortcutLabel() + ".x PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_yz: COMBINE ARG=" + getShortcutLabel() + ".y PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_zz: COMBINE ARG=" + getShortcutLabel() + ".z PERIODIC=NO"); 
    // Now calculate the total length of the vector
    createVectorNormInput( getShortcutLabel(), getShortcutLabel() + "_vsum", "z" );
  }
}

void PlaneShortcut::createVectorNormInput( const std::string& ilab, const std::string& olab, const std::string& vlab ) {
  readInputLine( olab + "2: COMBINE ARG1=" + ilab + "_x" + vlab + " ARG2=" + ilab + "_y" + vlab + " ARG3=" + ilab + "_z" + vlab + " POWERS=2,2,2 PERIODIC=NO");
  readInputLine( olab + ": MATHEVAL ARG1=" + olab + "2 FUNC=sqrt(x) PERIODIC=NO"); 
}

}
}



