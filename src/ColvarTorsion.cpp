#include "Colvar.h"
#include "ActionRegister.h"
#include "Torsion.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR TORSION
/**
Calculate the torsion between four atoms

\par Syntax
\verbatim
TORSION ATOMS=x,y,z,t [PBC]
\endverbatim
If the PBC flag is present, distance is computed using periodic boundary conditions.

*/
//+ENDPLUMEDOC
   
class ColvarTorsion : public Colvar {
  bool pbc;

public:
  ColvarTorsion(const ActionOptions&);
// active methods:
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(ColvarTorsion,"TORSION")

ColvarTorsion::ColvarTorsion(const ActionOptions&ao):
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

  if(atoms.size()==4){
    log.printf("  between atoms %d %d %d %d\n",atoms[0].serial(),atoms[1].serial(),atoms[2].serial(),atoms[3].serial());
    atoms.resize(6);
    atoms[5]=atoms[3];
    atoms[4]=atoms[2];
    atoms[3]=atoms[2];
    atoms[2]=atoms[1];
  }else assert(0);

  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");

  addValueWithDerivatives("");
  getValue("")->setPeriodicity(true);
  getValue("")->setDomain(-pi,pi);

  requestAtoms(atoms);
}

// calculator
void ColvarTorsion::calculate(){

  Vector d0,d1,d2;
  if(pbc){
    d0=pbcDistance(getPositions(1),getPositions(0));
    d1=pbcDistance(getPositions(3),getPositions(2));
    d2=pbcDistance(getPositions(5),getPositions(4));
  } else {
    d0=delta(getPositions(1),getPositions(0));
    d1=delta(getPositions(3),getPositions(2));
    d2=delta(getPositions(5),getPositions(4));
  }
  Vector dd0,dd1,dd2;
  Torsion t;
  double torsion=t.compute(d0,d1,d2,dd0,dd1,dd2);
  setAtomsDerivatives(0,dd0);
  setAtomsDerivatives(1,-dd0);
  setAtomsDerivatives(2,dd1);
  setAtomsDerivatives(3,-dd1);
  setAtomsDerivatives(4,dd2);
  setAtomsDerivatives(5,-dd2);

  setValue           (torsion);
  setBoxDerivatives  (-(extProduct(d0,dd0)+extProduct(d1,dd1)+extProduct(d2,dd2)));
}

}



