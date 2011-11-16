// This file is arse
/*#include "Colvar.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{
*/

//+PLUMEDOC COLVAR TEMPLATE
/**
This is just a template variable

*/
//+ENDPLUMEDOC

/*   
class ColvarTemplate : public Colvar {
  bool pbc;

public:
  ColvarTemplate(const ActionOptions&);
// active methods:
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(ColvarTemplate,"TEMPLATE")

ColvarTemplate::ColvarTemplate(const ActionOptions&ao):
Colvar(ao),
pbc(true)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  assert(atoms.size()==2);
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  parseFlag("PBC",pbc);
  checkRead();

  log.printf("  between atoms %d %d\n",atoms[0].serial(),atoms[1].serial());
  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");

  addValueWithDerivatives("");

  //requestAtoms(atoms);
}


// calculator
void ColvarTemplate::calculate(){

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
}

}*/



