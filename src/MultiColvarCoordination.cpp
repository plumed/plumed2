#include "MultiColvar.h"
#include "NeighborList.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR COORDINATIONNUMBER
/**
Calculate the coordination numbers of atoms so that  you can then calculate functions of the distribution of
coordination numbers such as the minimum, the number less than a certain quantity and so on.   

\par Examples

The following input tells plumed to calculate the coordination numbers of atoms 1-100 with themselves.
The minimum coordination number is then calculated.
\verbatim
COORDINATIONNUMBER SPECIES=1-100 R_0=1.0 MIN=0.1
\endverbatim

The following input tells plumed to calculate how many atoms from 1-100 are within 3.0 of each of the atoms
from 101-110.  In the first 101 is the central atom, in the second 102 is the central atom and so on.  The 
number of coordination numbers more than 6 is then computed.
\verbatim
COORDINATIONNUMBER SPECIESA=101-110 SPECIESB=1-100 R_0=3.0 MORE_THAN=6.0
\endverbatim

*/
//+ENDPLUMEDOC


class MultiColvarCoordination : public MultiColvar {
private:
  NeighborList *nl;
  SwitchingFunction switchingFunction;
public:
  static void registerKeywords( Keywords& keys );
  MultiColvarCoordination(const ActionOptions&);
  ~MultiColvarCoordination();
// active methods:
  virtual double compute( const std::vector<Vector>& pos, std::vector<Vector>& deriv, Tensor& virial );
};

PLUMED_REGISTER_ACTION(MultiColvarCoordination,"COORDINATIONNUMBER")

void MultiColvarCoordination::registerKeywords( Keywords& keys ){
  MultiColvar::registerKeywords( keys );
  ActionWithDistribution::autoParallelize( keys );
  keys.use("SPECIES"); keys.use("SPECIESA"); keys.use("SPECIESB");
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","12","The m parameter of the switching function ");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","NL_CUTOFF","The cutoff for the neighbor list");
  keys.remove("AVERAGE");
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

  // Create the groups for the neighbor list
  std::vector<AtomNumber> ga_lista, gb_lista; AtomNumber aa;
  aa.setIndex(0); ga_lista.push_back(aa);
  for(unsigned i=1;i<natoms;++i){ aa.setIndex(i); gb_lista.push_back(aa); }

  // Setup the neighbor list
  double nl_cut=-1.0;
  if( isTimeForNeighborListUpdate() ) parse("NL_CUTOFF",nl_cut); 
  if(nl_cut>0.0){
     nl = new NeighborList(ga_lista,gb_lista,false,usesPbc(),getPbc(),nl_cut,0);
     log.printf("  ignoring distances greater than %lf in neighbor list\n",nl_cut); 
  } else {
     nl = new NeighborList(ga_lista,gb_lista,false,usesPbc(),getPbc());
  }

  // And check everything has been read in correctly
  checkRead();
}

MultiColvarCoordination::~MultiColvarCoordination(){
  delete nl;
}

double MultiColvarCoordination::compute( const std::vector<Vector>& pos, std::vector<Vector>& deriv, Tensor& virial ){
   double value=0, dfunc; Vector distance;

   // Update the neighbor list at neighbor list update time
   if( isTimeForNeighborListUpdate() ){ nl->update(pos); updateAtoms( nl->getFullAtomList() ); }

   // Calculate the coordination number
   for(unsigned i=0;i<nl->size();++i){
      unsigned i0=nl->getClosePair(i).first;
      unsigned i1=nl->getClosePair(i).second;

      distance=getSeparation( pos[i0], pos[i1] );
      value += switchingFunction.calculate( distance.modulo(), dfunc );  
      deriv[0] = deriv[0] + (-dfunc)*distance;
      deriv[i] = deriv[i] + dfunc*distance;
      virial = virial + (-dfunc)*Tensor(distance,distance);
   }

   return value;
}

}

