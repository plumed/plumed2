#include "ColvarDistinguishable.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR DISTANCE
/**
Calculate distances between atoms.  To calculate minimum distances, the number of distances less or more than a given value etc you should use this
keyword in conjuction with the colvar modifiers described in \ref Colvar. 

\par Example

The following calculates the distance between atoms 3 and 5 and stores the value on d1.value0. 
\verbatim
DISTANCE ATOMS=3,5 LABEL=d1
\endverbatim

The following calculates two distances.  That between atoms 3 and 5 and that between atoms 2 and 4.  These two values
are stored on d1.value0 and d1.value1
\verbatim
DISTANCE ATOMS1=3,5 ATOMS2=2,4 LABEL=d1
\endverbatim
 
The following calculates all the distance between the atoms in the group - i.e. the distances between atoms 3 and 4, between
3 and 5 and between 4 and 5.

\verbatim
DISTANCE GROUP=3,4,5 LABEL=d1
\endverbatim

Lastly, this calculates all the distances between the atoms in the two groups - i.e. the distances between 3 and 4 and 3 and 5.

\verbatim
DISTNACE GROUP1=3 GROUP2=4,5 LABEL=d1
\endverbatim

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



