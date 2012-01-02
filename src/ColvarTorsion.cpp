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
Calculate the torsion between four atoms.
\par Syntax
\verbatim
TORSION ATOMS=a0,a1,a2,a3 [PBC]
\endverbatim
If the PBC flag is present, distance is computed using periodic boundary conditions.

Alternatively, compute the angle between two vectors projected on the plane
orthogonal to an axis, such as
\verbatim
TORSION V1=a0,a1 AXIS=a2,a3 V2=a4,a5
\endverbatim
Thus, two following variables are exactly the same one:
\verbatim
TORSION ATOMS=a0,a1,a2,a3
TORSION V1=a1,a0 AXIS=a1,a2 V2=a2,a3
\endverbatim

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
  vector<AtomNumber> atoms,v1,v2,axis;
  parseAtomList("ATOMS",atoms);
  parseAtomList("VECTOR1",v1);
  parseAtomList("VECTOR2",v2);
  parseAtomList("AXIS",axis);

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  parseFlag("PBC",pbc);
  checkRead();

  if(atoms.size()==4){
    assert(v1.size()==0 && v2.size()==0 && axis.size()==0);
    log.printf("  between atoms %d %d %d %d\n",atoms[0].serial(),atoms[1].serial(),atoms[2].serial(),atoms[3].serial());
    atoms.resize(6);
    atoms[5]=atoms[3];
    atoms[4]=atoms[2];
    atoms[3]=atoms[2];
    atoms[2]=atoms[1];
  }else if(atoms.size()==0){
    assert(v1.size()==2 && v2.size()==2 && axis.size()==2);
    log.printf("  between lines %d-%d and %d-%d, projected on the plane orthogonal to line %d-%d\n",
                v1[0].serial(),v1[1].serial(),v2[0].serial(),v2[1].serial(),axis[0].serial(),axis[1].serial());
    atoms.resize(6);
    atoms[0]=v1[1];
    atoms[1]=v1[0];
    atoms[2]=axis[0];
    atoms[3]=axis[1];
    atoms[4]=v2[0];
    atoms[5]=v2[1];
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



