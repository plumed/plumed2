#include "ColvarWithModifiers.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR ANGLE
/**
Calculates angle between sets of three atoms.  To calculate minimum angles, the number of angles less or more than a given value etc you should use this keyword
in conjuction with the colvar modifiers described in \ref Colvar.

\par Example

The following calculates the angle between the vectors connecting atoms 3 and 5 to atom 4 and stores the value on a1.value0.
\verbatim
ANGLE ATOMS=3,4,5 LABEL=a1
\endverbatim

The following calculates two angles.  That between the vectors connecting atoms 3 and 5 to atom 4 and that connecting atoms 6 and 8 to atom 7.  The values of these angles
are stored on a1.value0 and a1.value1 respectively.
\verbatim
ANGLE ATOMS1=3,4,5 ATOMS2=6,7,8 LABEL=a1
\endverbatim 

*/
//+ENDPLUMEDOC


class ColvarAngle : public ColvarWithModifiers {
public:
  ColvarAngle(const ActionOptions&);
// active methods:
  virtual double compute( const std::vector<unsigned>& indexes, std::vector<Vector>& derivatives, Tensor& virial );
};

PLUMED_REGISTER_ACTION(ColvarAngle,"ANGLE")

ColvarAngle::ColvarAngle(const ActionOptions&ao):
ColvarWithModifiers(ao)
{
  setNeighbourListStyle("none");
  registerKeyword(2,"ATOMS","calculate the angle defined by these three atoms.  To calculate multiple angles use ATOMS1, ATOMS2, ...");
  readActionAtomistic();
  int maxatoms=3; readAtomsKeyword(maxatoms);
  finishColvarSetup( 0, 0 );
}

double ColvarAngle::compute( const std::vector<unsigned>& indexes, std::vector<Vector>& derivatives, Tensor& virial ){
  assert( indexes.size()==3 && derivatives.size()==3 );
  virial.clear();

  Vector rji, rjk; 
  rji=getSeparation( indexes[1], indexes[0] ); rjk=getSeparation( indexes[1], indexes[2] );
  double caa=rji.modulo2();
  double cab=dot_product(rji,rjk); 
  double cbb=rjk.modulo2(); 
  
  double ccc = 1.0/sqrt(caa*cbb);
  double ac = cab*ccc; // cosine of teta
  double dVac = -ccc/sqrt(1.-ac*ac);

  for(unsigned ix=0;ix<3;ix++) {
    derivatives[0][ix] = dVac*(-cab/caa*rji[ix]+rjk[ix]);
    derivatives[1][ix] = -dVac*(-cab/caa*rji[ix]-cab/cbb*rjk[ix]+rji[ix]+rjk[ix]);
    derivatives[2][ix] = -dVac*(cab/cbb*rjk[ix]-rji[ix]);
  }
  virial = -1.0*( Tensor(derivatives[0],rji) + Tensor(derivatives[2],rjk) );
  return acos(ac);
}

}
