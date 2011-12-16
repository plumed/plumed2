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

\attention
For this variable the derivatives are computed numerically, so they might not
be very efficient!
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

// I am too lazy now to implement derivatives !!!
  enforceNumericalDerivatives();
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
// this is the version without derivatives
  double torsion=t.compute(d0,d1,d2);
  setValue           (torsion);
//  setAtomsDerivatives(0,d0);
//  setAtomsDerivatives(1,-d0);
//  setAtomsDerivatives(2,d1);
//  setAtomsDerivatives(3,-d1);
//  setAtomsDerivatives(4,d2);
//  setAtomsDerivatives(5,-d2);
//  setBoxDerivatives  (-(Tensor(d0,dd0)+Tensor(d1,dd1)+Tensor(d2,dd2)));
}

}



