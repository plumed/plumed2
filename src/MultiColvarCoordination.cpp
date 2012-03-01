#include "MultiColvar.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR COORDINATIONNUMBER
/**
Calculate the coordination numbers of atoms.  You can then calculate functions of the distribution of
coordination numbers such as the minimum, the number less than a certain quantity and so on. 

\par Examples

The following input tells plumed to calculate the coordination numbers of atoms 1-100 with themselves.
100 coordination numbers will be calculated in this example.
\verbatim
COORDINATIONNUMBER SPECIES=1-100 R_0=1.0
\endverbatim

The following input tells plumed to calculate how many atoms from 1-100 are within 3.0 of each of the atoms
from 101-110.  In the first 101 is the central atom, in the second 102 is the central atom and so on.
\verbatim
COORDINATIONNUMBER SPECIESA=101-110 SPECIESB=1-100 R_0=3.0
\endverbatim

*/
//+ENDPLUMEDOC


class MultiColvarCoordination : public MultiColvar {
private:
  SwitchingFunction switchingFunction;
public:
  static void registerKeywords( Keywords& keys );
  MultiColvarCoordination(const ActionOptions&);
// active methods:
  virtual double compute( const std::vector<Vector>& pos, std::vector<Vector>& deriv, Tensor& virial );
};

PLUMED_REGISTER_ACTION(MultiColvarCoordination,"COORDINATIONNUMBER")

void MultiColvarCoordination::registerKeywords( Keywords& keys ){
  MultiColvar::registerKeywords( keys );
  ActionWithDistribution::autoParallelize( keys );
  MultiColvar::useNeighbourList("sum",keys);
  keys.use("SPECIES"); keys.use("SPECIESA"); keys.use("SPECIESB");
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","12","The m parameter of the switching function ");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
}

MultiColvarCoordination::MultiColvarCoordination(const ActionOptions&ao):
PLUMED_MULTICOLVAR_INIT(ao)
{
  // Read in the switching function
  double r_0=-1.0, d_0; int nn, mm;
  parse("NN",nn); parse("MM",mm);
  parse("R_0",r_0); parse("D_0",d_0);
  if( r_0<0.0 ) error("you must set a value for R_0");
  switchingFunction.set(nn,mm,r_0,d_0);
  log.printf("  switching function parameters : r_0=%f, d_0=%f, nn=%d and mm=%d\n", r_0, d_0, nn, mm); 

  // Read in the atoms
  int natoms; readAtoms( natoms );

  // Setup the neighbour list
  std::vector< std::pair<unsigned,unsigned> > pairs(natoms-1);
  for(unsigned i=1;i<natoms;++i){ pairs[i-1].first=0; pairs[i-1].second=i; }
  createNeighbourList( pairs );
  // And check everything has been read in correctly
  checkRead();
}

double MultiColvarCoordination::compute( const std::vector<Vector>& pos, std::vector<Vector>& deriv, Tensor& virial ){
   double value=0, dfunc; Vector distance;
   for(unsigned i=1;i<pos.size();++i){
      distance=getSeparation( pos[0], pos[i] );
      value += switchingFunction.calculate( distance.modulo(), dfunc );  
      deriv[0] = deriv[0] + (-dfunc)*distance;
      deriv[i] = deriv[i] + dfunc*distance;
      virial = virial + (-dfunc)*Tensor(distance,distance);
   }
   return value;
}

}

