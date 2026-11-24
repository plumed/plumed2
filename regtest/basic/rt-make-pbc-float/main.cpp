#include "plumed/tools/Pbc.h"
#include "plumed/tools/Tools.h"
#include "plumed/tools/Random.h"
#include <cmath>
#include <fstream>

using namespace PLMD;

// Tools::pbc and Pbc::distance in single precision must be as accurate as the
// single precision representation of their arguments allows.
// Until v2.10 Tools::pbc<float> had an error up to 3.8e-6 (in scaled coordinates)
// because of an offset added to its argument.
int main() {
  std::ofstream out("output");
  Random rnd;
  rnd.setSeed(-2026);
  double maxerr=0.0;
  for(unsigned i=0; i<100000; i++) {
    const float x=float(3.0*rnd.RandU01()-1.5);
    double d=Tools::pbc(x)-(double(x)-std::nearbyint(double(x)));
    d-=std::nearbyint(d); // +0.5 and -0.5 are equivalent at exact ties
    maxerr=std::max(maxerr,std::fabs(d));
  }
  out<<"Tools::pbc<float> exact: "<<(maxerr==0.0?"yes":"no")<<"\n";

  Pbc pbc;
  // a triclinic box
  pbc.setBox(Tensor(5.0,0.0,0.0, 1.5,4.5,0.0, 1.0,1.2,4.0));
  double maxdiff=0.0;
  for(unsigned i=0; i<100000; i++) {
    Vector a,b;
    for(unsigned k=0; k<3; k++) {
      a[k]=20.0*rnd.RandU01()-10.0;
      b[k]=20.0*rnd.RandU01()-10.0;
    }
    const Vector dd=pbc.distance(a,b);
    const auto df=pbc.distance(a.convert<float>(),b.convert<float>());
    for(unsigned k=0; k<3; k++) {
      maxdiff=std::max(maxdiff,std::fabs(dd[k]-df[k]));
    }
  }
  // positions up to 10 in single precision are accurate to ~5e-7
  out<<"Pbc::distance single vs double precision within 5e-6: "<<(maxdiff<5e-6?"yes":"no")<<"\n";
  return 0;
}
