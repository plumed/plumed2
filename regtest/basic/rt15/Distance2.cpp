/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#include "colvar/Colvar.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{

class Distance : public colvar::Colvar {
  bool components;
  bool pbc;

public:
  static void registerKeywords( Keywords& keys );
  Distance(const ActionOptions&);
// active methods:
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(Distance,"DISTANCE2")

void Distance::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOMS","the pair of atom that we are calculating the distance between");
  keys.addFlag("COMPONENTS",false,"calculate the x, y and z components of the distance separately and store them as label.x, label.y and label.z");  
}

Distance::Distance(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
components(false),
pbc(true)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  if(atoms.size()!=2)
    error("Number of specified atoms should be 2");
  parseFlag("COMPONENTS",components);
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  checkRead();

  log.printf("  between atoms %d %d\n",atoms[0].serial(),atoms[1].serial());
  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");


  if(!components){

    addValueWithDerivatives(); setNotPeriodic();

  } else{
    addComponentWithDerivatives("x"); componentIsNotPeriodic("x");
    addComponentWithDerivatives("y"); componentIsNotPeriodic("y");
    addComponentWithDerivatives("z"); componentIsNotPeriodic("z");
  }

  requestAtoms(atoms);
}


// calculator
void Distance::calculate(){

  Vector distance;
  if(pbc){
    distance=pbcDistance(getPosition(0),getPosition(1));
  } else {
    distance=delta(getPosition(0),getPosition(1));
  }
  const double value=distance.modulo();
  const double invvalue=1.0/value;

  if(!components){

    setAtomsDerivatives(0,-invvalue*distance);
    setAtomsDerivatives(1,invvalue*distance);
    setBoxDerivatives  (-invvalue*Tensor(distance,distance));
    setValue           (value);

  }else{

    Value* valuex=getPntrToComponent("x");
    Value* valuey=getPntrToComponent("y");
    Value* valuez=getPntrToComponent("z");

    setAtomsDerivatives (valuex,0,Vector(-1,0,0));
    setAtomsDerivatives (valuex,1,Vector(+1,0,0));
    setBoxDerivatives   (valuex,Tensor(distance,Vector(-1,0,0)));
    valuex->set(distance[0]);

    setAtomsDerivatives (valuey,0,Vector(0,-1,0));
    setAtomsDerivatives (valuey,1,Vector(0,+1,0));
    setBoxDerivatives   (valuey,Tensor(distance,Vector(0,-1,0)));
    valuey->set(distance[1]);

    setAtomsDerivatives (valuez,0,Vector(0,0,-1));
    setAtomsDerivatives (valuez,1,Vector(0,0,+1));
    setBoxDerivatives   (valuez,Tensor(distance,Vector(0,0,-1)));
    valuez->set(distance[2]);
  };
}

}



