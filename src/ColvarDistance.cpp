#include "Colvar.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR DISTANCE
/**
Calculate the distance between two atoms.

\par Syntax
\verbatim
DISTANCE ATOMS=x,y [COMPONENTS] [PBC]
\endverbatim
If the COMPONENTS flag is present, the three components of the distance
can be accessed respectively as label.x label.y and label.z .
If the PBC flag is present, distance is computed using periodic boundary conditions.

\par Example
The following input is printing the distance between atoms 3 and 5,
the distance between atoms 2 and 4 and its x component.
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

public:
  ColvarDistance(const ActionOptions&);
// active methods:
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(ColvarDistance,"DISTANCE")

ColvarDistance::ColvarDistance(const ActionOptions&ao):
Colvar(ao),
components(false)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  assert(atoms.size()==2);
  parseFlag("COMPONENTS",components);
  readActionAtomistic();
  checkRead();

  log.printf("  between atoms %d %d\n",atoms[0].serial(),atoms[1].serial());
//  if(pbc) log.printf("  using periodic boundary conditions\n");
//  else    log.printf("  without periodic boundary conditions\n");


  if(!components){

    addValueWithDerivatives("");
    getValue("")->setPeriodicity(false);

  }else{

    addValueWithDerivatives("x");
    getValue("x")->setPeriodicity(false);
    addValueWithDerivatives("y");
    getValue("y")->setPeriodicity(false);
    addValueWithDerivatives("z");
    getValue("z")->setPeriodicity(false);
  }
}


// calculator
void ColvarDistance::calculate(){

  Vector distance=getSeparation(0,1);
  // if(pbc){
  //   distance=pbcDistance(getPositions(0),getPositions(1));
  // } else {
  //   distance=delta(getPositions(0),getPositions(1));
  // }
  const double value=distance.modulo();
  const double invvalue=1.0/value;

  if(!components){

    setAtomsDerivatives(0,-invvalue*distance);
    setAtomsDerivatives(1,invvalue*distance);
    setBoxDerivatives  (-invvalue*Tensor(distance,distance));
    setValue           (value);

  }else{

    Value* valuex=getValue("x");
    Value* valuey=getValue("y");
    Value* valuez=getValue("z");

    setAtomsDerivatives (valuex,0,Vector(-1,0,0));
    setAtomsDerivatives (valuex,1,Vector(+1,0,0));
    setBoxDerivatives   (valuex,Tensor(distance,Vector(-1,0,0)));
    setValue            (valuex,distance[0]);

    setAtomsDerivatives (valuey,0,Vector(0,-1,0));
    setAtomsDerivatives (valuey,1,Vector(0,+1,0));
    setBoxDerivatives   (valuey,Tensor(distance,Vector(0,-1,0)));
    setValue            (valuey,distance[1]);

    setAtomsDerivatives (valuez,0,Vector(0,0,-1));
    setAtomsDerivatives (valuez,1,Vector(0,0,+1));
    setBoxDerivatives   (valuez,Tensor(distance,Vector(0,0,-1)));
    setValue            (valuez,distance[2]);
  };
}

}



