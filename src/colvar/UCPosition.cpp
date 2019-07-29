/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2019 The plumed team
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
#include "tools/Pbc.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR UCPOSITION
/*

*/
//+ENDPLUMEDOC

class UCPosition : public Colvar {
  vector<double> lattice_;

public:
  static void registerKeywords( Keywords& keys );
  explicit UCPosition(const ActionOptions&);
// active methods:
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(UCPosition,"UCPOSITION")

void UCPosition::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOM","the atom of which we are calculating the position in the unit cell");
  keys.add("compulsory","LATTICE","lattice parameters in x and y directions");
  keys.addOutputComponent("x","COMPONENTS","the x coordinate of lattice position");
  keys.addOutputComponent("y","COMPONENTS","the y coordinate of lattice position");
  keys.addOutputComponent("z","COMPONENTS","the z coordinate of lattice position");
}

UCPosition::UCPosition(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOM",atoms);
  if(atoms.size()!=1)
    error("Number of specified atoms should be 1");

  parseVector("LATTICE", lattice_);
  if(lattice_.size()!=2)
    error("LATTICE should be specified by 2 real numbers");

  checkRead();

  log.printf("  for atom %d \n",atoms[0].serial());
  log.printf("  lattice parameters in x and y directions: %lf %lf\n",lattice_[0],lattice_[1]);

  // convert to string
  std::string l0; Tools::convert(lattice_[0], l0);
  std::string l1; Tools::convert(lattice_[1], l1);
  // add components
  addComponentWithDerivatives("x"); componentIsPeriodic("x", "0.0", l0);
  addComponentWithDerivatives("y"); componentIsPeriodic("y", "0.0", l1);
  addComponentWithDerivatives("z"); componentIsNotPeriodic("z");

  requestAtoms(atoms);
}


// calculator
void UCPosition::calculate() {

  Vector pos = getPosition(0);

  Value* valuex=getPntrToComponent("x");
  Value* valuey=getPntrToComponent("y");
  Value* valuez=getPntrToComponent("z");

  setAtomsDerivatives (valuex,0,Vector(+1,0,0));
  setBoxDerivativesNoPbc(valuex);
  // calculate x-distance modulo lattice parameter
  double dx = pos[0] - lattice_[0] * static_cast<int>( pos[0] / lattice_[0] );
  if(dx<0) dx += lattice_[0];
  valuex->set(dx);

  setAtomsDerivatives (valuey,0,Vector(0,+1,0));
  setBoxDerivativesNoPbc(valuey);
  // calculate y-distance modulo lattice parameter
  double dy = pos[1] - lattice_[1] * static_cast<int>( pos[1] / lattice_[1] );
  if(dy<0) dy += lattice_[1];
  valuey->set(dy);

  setAtomsDerivatives (valuez,0,Vector(0,0,+1));
  setBoxDerivativesNoPbc(valuez);
  valuez->set(pos[2]);

}

}
}
