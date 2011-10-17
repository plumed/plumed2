#include "ColvarDistinguishable.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR TORSION 
/**
Calculates torsional angles defined by the positions of four atoms.  To calculate minimum torsions, the number of torsions less or more than a given value etc you should use this keyword
in conjuction with the colvar modifiers described in \ref Colvar.

\par Example

The following calculates the torsion between the plane containing atoms 3, 4 and 5 and the plane containing atoms 4, 5 and 6. The output is stored on t1.value0.
\verbatim
TORSION ATOMS=3,4,5,6 LABEL=t1
\endverbatim

The following calculates two torsional angle.  One involves atoms 3, 4, 5 and 6 while the other involves atoms 7, 8, 9 and 10.  The values of these angles
are stored on t1.value0 and t1.value1 respectively.
\verbatim
ANGLE ATOMS1=3,4,5,6 ATOMS2=7,8,9,10 LABEL=a1
\endverbatim 

*/
//+ENDPLUMEDOC

class ColvarTorsion : public ColvarDistinguishable {
  double pi,twopi;
public:
  ColvarTorsion(const ActionOptions&);
// active methods:
  virtual double compute( const std::vector<unsigned>& indexes, std::vector<Vector>& derivatives, Tensor& virial );
};

PLUMED_REGISTER_ACTION(ColvarTorsion,"TORSION")

ColvarTorsion::ColvarTorsion(const ActionOptions&ao):
ColvarDistinguishable(ao)
{
  std::vector<double> domain( 2, 0.0 );
  Tools::convert("+PI",pi); Tools::convert("+2PI",twopi);
  Tools::convert("-PI",domain[0]); Tools::convert("+PI",domain[1]); 
  readActionColvar( 4, domain );
  checkRead();
}

double ColvarTorsion::compute( const std::vector<unsigned>& indexes, std::vector<Vector>& derivatives, Tensor& virial ){
  assert( indexes.size()==4 && derivatives.size()==4 );
  virial.clear();

  Vector rij=getSeparation(indexes[0],indexes[1]);//  minimal_image(xi, xj, &mod_rij, r_ij);
  Vector rkj=getSeparation(indexes[2],indexes[1]);//  minimal_image(xk, xj, &mod_rkj, r_kj);
  Vector rkl=getSeparation(indexes[2],indexes[3]);//  minimal_image(xk, xl, &mod_rkl, r_kl);
  Vector rlj=getSeparation(indexes[3],indexes[1]);

  double n21=rkj.modulo(); double in21=1.0/(n21*n21);
  Vector vec1=cross_product(rij,rkj); double nv1=vec1.modulo(); double sc1=in21*dot_product(rij,rkj); //oprod(r_ij,r_kj,m);               
  Vector vec2=cross_product(rkj,rkl); double nv2=vec2.modulo(); double sc2=in21*dot_product(rkj,rkl); //oprod(r_kj,r_kl,n);
  double cos_phi=dot_product(vec1,vec2)/( nv1 * nv2 );  //  *cos_phi=cos_angle(m,n);

  for(unsigned i=0;i<3;++i) {
    derivatives[0][i]= -n21 * vec1[i] / (nv1*nv1);
    derivatives[3][i]=  n21 * vec2[i] / (nv2*nv2);
    derivatives[1][i]=  (sc1-1.0)*derivatives[0][i] - sc2*derivatives[3][i];
    derivatives[2][i]=  (sc2-1.0)*derivatives[3][i] - sc1*derivatives[0][i];
  }
  virial = Tensor(derivatives[0],rij) + Tensor(derivatives[2],rkj) + Tensor(derivatives[3],rlj);
  double psi=(dot_product(rij,vec2)<0.0)?-1.0*acos( cos_phi ):acos( cos_phi );

  if (psi >= pi) psi -= twopi;
  else if(psi < -pi) psi += twopi;
  return psi;
}

}
