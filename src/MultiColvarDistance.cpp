#include "MultiColvar.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR MULTIDISTANCE
/**
Calculate the distances between one or many pairs of atoms.  You can then calculate functions of the distribution of
distances such as the minimum, the number less than a certain quantity and so on. 

\par Examples

The following input tells plumed to print the distance between atoms 3 and 5,
\verbatim
MULTIDISTANCE ATOMS=3,5 LABEL=d1
PRINT ARG=d1.*
\endverbatim
(See also \ref PRINT).

The following input tells plumed to print the distances between atoms 3 and 5 and between atoms 1 and 2
\verbatim
MULTIDISTANCE ATOMS1=3,5 ATOMS2=1,2 LABEL=d1
PRINT ARG=d1.*
\endverbatim
(See also \ref PRINT).

The following input tells plumed to calculate the distances between atoms 3 and 5 and between atoms 1 and 2
and then to calculate the number of these distances that are less than 0.1 nm.  The number of distances
less than 0.1nm is then printed to a file.
\verbatim
MULTIDISTANCE ATOMS1=3,5 ATOMS2=1,2 LABEL=d1 LESS_THAN=0.1
PRINT ARG=d1.lt0.1
\endverbatim
(See also \ref PRINT).

*/
//+ENDPLUMEDOC


class MultiColvarDistance : public MultiColvar {
public:
  static void registerKeywords( Keywords& keys );
  MultiColvarDistance(const ActionOptions&);
// active methods:
  virtual double compute( const std::vector<Vector>& pos, std::vector<Vector>& deriv, Tensor& virial );
};

PLUMED_REGISTER_ACTION(MultiColvarDistance,"MULTIDISTANCE")

void MultiColvarDistance::registerKeywords( Keywords& keys ){
  MultiColvar::registerKeywords( keys );
  MultiColvar::useNeighbourList("product",keys);
  keys.use("ATOMS"); keys.use("GROUP"); keys.use("GROUPA"); keys.use("GROUPB");
}

MultiColvarDistance::MultiColvarDistance(const ActionOptions&ao):
PLUMED_MULTICOLVAR_INIT(ao)
{
  // Read in the atoms
  int natoms=2; readAtoms( natoms );
  // Setup the neighbour list
  std::vector< std::pair<unsigned,unsigned> > pairs(1);
  pairs[0].first=0; pairs[0].second=1;
  createNeighbourList( pairs );
  // And check everything has been read in correctly
  checkRead();
}

double MultiColvarDistance::compute( const std::vector<Vector>& pos, std::vector<Vector>& deriv, Tensor& virial ){
   Vector distance; 
   distance=getSeparation( pos[0], pos[1] );
   const double value=distance.modulo();
   const double invvalue=1.0/value; 

   deriv[0]=-invvalue*distance;
   deriv[1]=invvalue*distance;
   virial=-invvalue*Tensor(distance,distance);
   return value;
}

}

