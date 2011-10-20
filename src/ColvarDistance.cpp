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
  int component;
public:
  ColvarDistance(const ActionOptions&);
// active methods:
  virtual double compute( const std::vector<unsigned>& indexes, std::vector<Vector>& derivatives, Tensor& virial );
};

PLUMED_REGISTER_ACTION(ColvarDistance,"DISTANCE")

ColvarDistance::ColvarDistance(const ActionOptions&ao):
ColvarDistinguishable(ao),
component(-1)
{
  allowKeyword("GROUP" );
  registerKeyword(0,"COMPONENT","use this if you only want the X,Y or Z component of the distance");
  std::vector<double> domain( 2, 0.0 );
  readActionColvar( 2, domain );
  
  std::string lab="none"; parse("COMPONENT",lab);
  if( lab=="X" || lab=="x" ){
     component=0;
  } else if( lab=="Y" || lab=="y" ){
     component=1;
  } else if( lab=="Z" || lab=="z" ){
     component=2;
  } else if( lab!="none") {
      error( lab + " is not a valid argument for the COMPONENT keyword use X, Y or Z");
  } 
  checkRead();
}

double ColvarDistance::compute( const std::vector<unsigned>& indexes, std::vector<Vector>& derivatives, Tensor& virial ){
  assert( indexes.size()==2 && derivatives.size()==2 );
  Vector distance=getSeparation( indexes[0], indexes[1] ); 

  if ( component<0 ){
    const double value=distance.modulo();
    const double invvalue=1.0/value;
    derivatives[0]=-invvalue*distance;
    derivatives[1]=invvalue*distance;
    virial=-invvalue*Tensor(distance,distance);
    return value;
  } else {
    const double value=distance[component];
    const double invvalue=1.0/value;
    derivatives[0][component]=-1.0; 
    derivatives[1][component]=1.0; 
    virial=Tensor( distance,derivatives[0] );
    return value;
  }
}

}



