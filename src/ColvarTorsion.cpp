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
Calculate the torsion between four atoms or alternatively use this command
to calculate the angle between two vectors projected on the plane
orthogonal to an axis. 

\par Examples

The torsional angle between atoms 1, 2, 3 and 4 can be calculated using:

\verbatim
TORSION ATOMS=1,2,3,4
\endverbatim

or alternatively using:

\verbatim
TORSION VECTOR1=2,1 AXIS=2,3 VECTOR2=3,4
\endverbatim

*/
//+ENDPLUMEDOC
   
class ColvarTorsion : public Colvar {
  bool pbc;

public:
  ColvarTorsion(const ActionOptions&);
// active methods:
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(ColvarTorsion,"TORSION")

void ColvarTorsion::registerKeywords(Keywords& keys){
   Colvar::registerKeywords( keys );
   keys.add("atoms","ATOMS","the four atoms involved in the torsional angle");
   keys.add("atoms","AXIS","two atoms that define an axis.  You can use this to find the angle in the plane perpendicular to the axis between the vectors specified using the VECTOR1 and VECTOR2 keywords."); 
   keys.add("atoms","VECTOR1","two atoms that define a vector.  You can use this in combination with VECTOR2 and AXIS");
   keys.add("atoms","VECTOR2","two atoms that define a vector.  You can use this in combination with VECTOR1 and AXIS");
}

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

  addValueWithDerivatives(); setPeriodic(-pi,pi);
  requestAtoms(atoms);
}

// calculator
void ColvarTorsion::calculate(){

  Vector d0,d1,d2;
  if(pbc){
    d0=pbcDistance(getPosition(1),getPosition(0));
    d1=pbcDistance(getPosition(3),getPosition(2));
    d2=pbcDistance(getPosition(5),getPosition(4));
  } else {
    d0=delta(getPosition(1),getPosition(0));
    d1=delta(getPosition(3),getPosition(2));
    d2=delta(getPosition(5),getPosition(4));
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



