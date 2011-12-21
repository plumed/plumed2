#include "Colvar.h"
#include "ActionRegister.h"
#include "Angle.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR ANGLE
/**
Calculate the angle between three atoms

\par Syntax
\verbatim
ANGLE ATOMS=x,y,z [PBC]
\endverbatim
If the PBC flag is present, distance is computed using periodic boundary conditions.

If *four* atoms are passed, the angle between the vector joining atoms 1-2
and the vector joining atoms 3-4 is computed

*/
//+ENDPLUMEDOC
   
class ColvarAngle : public Colvar {
  bool pbc;

public:
  ColvarAngle(const ActionOptions&);
// active methods:
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(ColvarAngle,"ANGLE")

ColvarAngle::ColvarAngle(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  parseFlag("PBC",pbc);
  checkRead();

  if(atoms.size()==3){
    log.printf("  between atoms %d %d %d\n",atoms[0].serial(),atoms[1].serial(),atoms[2].serial());
    atoms.resize(4);
    atoms[3]=atoms[2];
    atoms[2]=atoms[1];
  }else if(atoms.size()==4){
    log.printf("  between lines %d-%d and %d-%d\n",atoms[0].serial(),atoms[1].serial(),atoms[2].serial(),atoms[3].serial());
  }else assert(0);

  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");

  addValueWithDerivatives("");
  getValue("")->setPeriodicity(false);

  requestAtoms(atoms);
}

// calculator
void ColvarAngle::calculate(){

  Vector dij,dik;
  if(pbc){
    dij=pbcDistance(getPositions(2),getPositions(3));
    dik=pbcDistance(getPositions(1),getPositions(0));
  } else {
    dij=delta(getPositions(2),getPositions(3));
    dik=delta(getPositions(1),getPositions(0));
  }
  Vector ddij,ddik;
  Angle a;
  double angle=a.compute(dij,dik,ddij,ddik);
  setAtomsDerivatives(0,ddik);
  setAtomsDerivatives(1,-ddik);
  setAtomsDerivatives(2,-ddij);
  setAtomsDerivatives(3,ddij);
  setValue           (angle);
  setBoxDerivatives  (-(Tensor(dij,ddij)+Tensor(dik,ddik)));
}

}



