#include "SwitchingFunction.h"
#include "ColvarCoordination.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{ 

//+PLUMEDOC COLVAR COORDINATION
/**
This keyword can be used to calculate the coordination numbers for atoms in your system.  Minimum coordination numbers, maximum coordination numbers and so on can be calculated
using the colvar modifiers described in \ref Colvar.

\par Example

The following example instructs plumed to calculate coordination numbers for all the atoms in the group.  These coordination numbers count the number of atoms within the group that are within 0.3 nm of the central atom.  A neighbour list is used to make this calculation faster, this neighbour list is updated every 100 steps.
\verbatim
COORDINATION GROUP=1-100 R_0=0.3 NL_CUTOFF=0.5 UPDATE=100 
\endverbatim

The following example instructs plumed to calculate the coordination numbers between the two groups.  Each of these coordination numbers calculate the number of atoms from GROUP2 that are within 0.3 nm of one of the atoms in GROUP1.
\verbatim
COORDINATIoN GROUP1=1-50 GROUP2=51-300 NL_CUTOFF=0.5 UPDATE=100
\endverbatim 

*/
//+ENDPLUMEDOC

class ColvarCoordinationNumber : public ColvarCoordination {
  SwitchingFunction switchingFunction;
public:
  ColvarCoordinationNumber(const ActionOptions&);
// active methods:
  virtual double compute( const unsigned& nvecs, const Vector& cpos, const std::vector<Vector>& npos, std::vector<Vector>& derivatives, Tensor& vir ); 
};

PLUMED_REGISTER_ACTION(ColvarCoordinationNumber,"COORDINATION")

ColvarCoordinationNumber::ColvarCoordinationNumber(const ActionOptions&ao):
ColvarCoordination(ao)
{

  registerKeyword(1,"R_0", "the location of the inflection point for the switching function");
  registerKeyword(0,"D_0", "(default=0.0) for all r values less than this value the switching function equals 1 and is constant");
  registerKeyword(0,"NN", "(default=6) the exponent for the factor of (r-d_0)/r_0 that appears in the numerator of the switching function");
  registerKeyword(0,"MM", "(default=12) the exponent for the factor of (r-d_0)/r_0 that appears in the denominator of the switching function");
  std::vector<double> domain( 2, 0.0 );
  readActionColvar( 0, domain ); 
 
  int nn=6; parse("NN",nn); 
  int mm=12; parse("MM",mm);
  double r_0=-1; parse("R_0",r_0);
  double d_0=0.0; parse("D_0",d_0);
  switchingFunction.set(nn,mm,r_0,d_0);
  checkRead();
}

// calculator
double ColvarCoordinationNumber::compute( const unsigned& nvecs, const Vector& cpos, const std::vector<Vector>& npos, std::vector<Vector>& derivatives, Tensor& vir ) {
 assert( (nvecs+1)<=derivatives.size() && nvecs<=npos.size() );

 double ncoord=0, dfunc; Vector distance;
 for(unsigned n=0; n<nvecs;++n){
    distance=getSeparation( cpos, npos[n] );
    dfunc=0.; ncoord += switchingFunction.calculate(distance.modulo(), dfunc);
    derivatives[0] = derivatives[0] + (-dfunc)*distance;
    derivatives[n+1] = derivatives[n+1] + dfunc*distance;
    vir = vir+(-dfunc)*Tensor(distance,distance); 
 }
 return ncoord;
}

}
