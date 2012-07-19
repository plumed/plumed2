/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The PLUMED team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of PLUMED, version 2.0.

   PLUMED is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   PLUMED is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with PLUMED.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Colvar.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR DISTANCE
/*
Calculate the distance between a pair of atoms.

\par Examples

The following input tells plumed to print the distance between atoms 3 and 5,
the distance between atoms 2 and 4 and the x component of the distance between atoms 2 and 4.
\verbatim
DISTANCE ATOMS=3,5             LABEL=d1
DISTANCE ATOMS=2,4 COMPONENTS  LABEL=d2
PRINT ARG=d1,d2,d2.x
\endverbatim
(See also \ref PRINT).

*/
//+ENDPLUMEDOC
   
class ColvarDistance : public Colvar {
  bool components;
  bool pbc;

public:
  static void registerKeywords( Keywords& keys );
  ColvarDistance(const ActionOptions&);
// active methods:
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(ColvarDistance,"DISTANCE")

void ColvarDistance::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOMS","the pair of atom that we are calculating the distance between");
  keys.addFlag("COMPONENTS",false,"calculate the x, y and z components of the distance separately and store them as label.x, label.y and label.z");  
}

ColvarDistance::ColvarDistance(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
components(false),
pbc(true)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  assert(atoms.size()==2);
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
void ColvarDistance::calculate(){

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



