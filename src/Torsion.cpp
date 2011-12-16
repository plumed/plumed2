#include "Torsion.h"

#include <cmath>

using namespace PLMD;

double Torsion::compute(const Vector& v1,const Vector& v2,const Vector& v3)const{
    const Vector nv2(v2*(1.0/v2.modulo()));
    const Vector a(crossProduct(nv2,v1));
    const Vector b(crossProduct(v3,nv2));
    const double cosangle=dotProduct(a,b);
    const double sinangle=dotProduct(crossProduct(a,b),nv2);
    return std::atan2(-sinangle,cosangle);
}


