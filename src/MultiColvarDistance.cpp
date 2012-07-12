#include "MultiColvar.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC MCOLVAR DISTANCES
/*
Calculate the distances between one or many pairs of atoms.  You can then calculate functions of the distribution of
distances such as the minimum, the number less than a certain quantity and so on. 

\par Examples

The following input tells plumed to print the distance between atoms 3 and 5,
\verbatim
DISTANCES ATOMS=3,5 LABEL=d1
PRINT ARG=d1.*
\endverbatim
(See also \ref PRINT).

The following input tells plumed to print the distances between atoms 3 and 5 and between atoms 1 and 2
\verbatim
DISTANCES ATOMS1=3,5 ATOMS2=1,2 LABEL=d1
PRINT ARG=d1.*
\endverbatim
(See also \ref PRINT).

The following input tells plumed to calculate the distances between atoms 3 and 5 and between atoms 1 and 2
and then to calculate the number of these distances that are less than 0.1 nm.  The number of distances
less than 0.1nm is then printed to a file.
\verbatim
DISTANCES ATOMS1=3,5 ATOMS2=1,2 LABEL=d1 LESS_THAN=0.1
PRINT ARG=d1.lt0.1
\endverbatim
(See also \ref PRINT).

*/
//+ENDPLUMEDOC


class MultiColvarDistance : public MultiColvar {
private:
  double rcut;
public:
  static void registerKeywords( Keywords& keys );
  MultiColvarDistance(const ActionOptions&);
// active methods:
  virtual double compute( const unsigned& j, const std::vector<Vector>& pos, std::vector<Vector>& deriv, Tensor& virial );
/// Returns the number of coordinates of the field
  unsigned getNumberOfFieldDerivatives();
  bool isPeriodic(const unsigned nn){ return false; }
};

PLUMED_REGISTER_ACTION(MultiColvarDistance,"DISTANCES")

void MultiColvarDistance::registerKeywords( Keywords& keys ){
  MultiColvar::registerKeywords( keys );
  ActionWithDistribution::autoParallelize( keys );
  keys.add("optional","NL_CUTOFF","The cutoff for the neighbor list");
  keys.use("ATOMS"); keys.use("GROUP"); keys.use("GROUPA"); keys.use("GROUPB");
  keys.use("DISTRIBUTION");
}

MultiColvarDistance::MultiColvarDistance(const ActionOptions&ao):
PLUMED_MULTICOLVAR_INIT(ao),
rcut(-1)
{
  // Read in the atoms
  int natoms=2; readAtoms( natoms );
  // And setup the ActionWithDistribution
  requestDistribution();          
  // Read the cutoff for the neighbour list
  if( isTimeForNeighborListUpdate() ){
      parse("NL_CUTOFF",rcut);
      if( rcut>0 ) log.printf("  ignoring distances greater than %lf in neighbor list\n",rcut);
  }
  // And check everything has been read in correctly
  checkRead();
}

unsigned MultiColvarDistance::getNumberOfFieldDerivatives(){
  return 3*getNumberOfAtoms() + 9;
} 

double MultiColvarDistance::compute( const unsigned& j, const std::vector<Vector>& pos, std::vector<Vector>& deriv, Tensor& virial ){
   Vector distance; 
   distance=getSeparation( pos[0], pos[1] );
   const double value=distance.modulo();
   const double invvalue=1.0/value;

   // Check at neighbor list update time whether this distance is big
   if( isTimeForNeighborListUpdate() && rcut>0 ){
       if( value>rcut ){ stopCalculatingThisCV(); return 0.0; }
   } 

   // And finish the calculation
   deriv[0]=-invvalue*distance;
   deriv[1]=invvalue*distance;
   virial=-invvalue*Tensor(distance,distance);
   return value;
}

}

