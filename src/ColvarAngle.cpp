#include "ColvarDistinguishable.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

class ColvarAngle : public ColvarDistinguishable {
public:
  ColvarAngle(const ActionOptions&);
// active methods:
  virtual double compute( const std::vector<unsigned>& indexes, std::vector<Vector>& derivatives, Tensor& virial );
};

PLUMED_REGISTER_ACTION(ColvarAngle,"ANGLE")

ColvarAngle::ColvarAngle(const ActionOptions&ao):
ColvarDistinguishable(ao)
{
  std::vector<double> domain( 2, 0.0 );
  readActionColvar( 3, domain );
  checkRead();
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
