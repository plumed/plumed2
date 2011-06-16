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
DISTANCE ATOMS=x,y [COMPONENTS] [PBC]
If the COMPONENTS flag is present, the three components of the distance
can be accessed respectively as .x .y and .z.
If the PBC flag is present, distance is computed using periodic boundary conditions.

\par Example
The following input is printing the distance between atoms 3 and 5,
the distance between atoms 2 and 4 and its x component.
\verbatim
DISTANCE ATOMS=3,5             LABEL=d1
DISTANCE ATOMS=2,4 COMPONENTS  LABEL=d2
PRINT ARG=d1,d2,d2.x
\endverbatim

*/
//+ENDPLUMEDOC
   
class ColvarDistance : public Colvar {
  bool components;
  bool pbc;

public:
  ColvarDistance(const ActionOptions&);
// active methods:
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(ColvarDistance,"DISTANCE")

ColvarDistance::ColvarDistance(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
components(false),
pbc(true)
{
  vector<int> atoms;
  parseVector("ATOMS",atoms);
  assert(atoms.size()==2);
  parseFlag("COMPONENTS",components);
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  parseFlag("PBC",pbc);
  checkRead();

  log.printf("  between atoms %d %d\n",atoms[0],atoms[1]);
  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");

  addValueWithDerivatives("");

  if(components){
    addValueWithDerivatives("x");
    addValueWithDerivatives("y");
    addValueWithDerivatives("z");
  }

  requestAtoms(atoms);
}


// calculator
void ColvarDistance::calculate(){

  Vector distance;
  if(pbc){
    distance=pbcDistance(getPositions(0),getPositions(1));
  } else {
    distance=delta(getPositions(0),getPositions(1));
  }
  const double value=distance.modulo();
  const double invvalue=1.0/value;

  setAtomsDerivatives(0,-invvalue*distance);
  setAtomsDerivatives(1,invvalue*distance);
  setBoxDerivatives  (-invvalue*Tensor(distance,distance));
  setValue           (value);

  if(components){
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



