#include "ColvarDistinguishable.h"
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
   
class ColvarDistance : public ColvarDistinguishable {
public:
  ColvarDistance(const ActionOptions&);
// active methods:
  virtual double compute( const std::vector<unsigned>& indexes, std::vector<Vector>& derivatives, Tensor& virial );
};

PLUMED_REGISTER_ACTION(ColvarDistance,"DISTANCE")

ColvarDistance::ColvarDistance(const ActionOptions&ao):
ColvarDistinguishable(ao)
{
  allowKeyword("GROUP" );
  std::vector<double> domain( 2, 0.0 );
  readActionColvar( 2, domain );
  checkRead();
}

double ColvarDistance::compute( const std::vector<unsigned>& indexes, std::vector<Vector>& derivatives, Tensor& virial ){
  assert( indexes.size()==2 && derivatives.size()==2 );
  Vector distance=getSeparation( indexes[0], indexes[1] ); 
  const double value=distance.modulo();
  const double invvalue=1.0/value;
  derivatives[0]=-invvalue*distance;
  derivatives[1]=invvalue*distance;
  virial=-invvalue*Tensor(distance,distance);
  return value;  
}

}



