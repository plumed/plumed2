#include "Angle.h"
#include "Tools.h"
#include <cmath>

using namespace PLMD;

double Angle::compute(const Vector& v1,const Vector& v2)const{
    return std::acos(dotProduct(v1,v2)/(v1.modulo()*v2.modulo()));
}

double Angle::compute(const Vector& v1,const Vector& v2,Vector& d1,Vector& d2)const{
    const double dp(dotProduct(v1,v2));
    const Vector& dp_dv1(v2);
    const Vector& dp_dv2(v1);
    const double sv1(v1.modulo2());
    const double sv2(v2.modulo2());
    const Vector dsv1_dv1(2*v1);
    const Vector dsv2_dv2(2*v2);
    const double nn(1.0/sqrt(sv1*sv2));
    const Vector dnn_dv1(-0.5*nn/sv1*dsv1_dv1);
    const Vector dnn_dv2(-0.5*nn/sv2*dsv2_dv2);

    const double dpnn(dp*nn);

    if(dpnn==1.0){
      d1=Vector(0.0,0.0,0.0);
      d2=Vector(0.0,0.0,0.0);
      return 0.0;
    }
    if(dpnn==-1.0){
      d1=Vector(0.0,0.0,0.0);
      d2=Vector(0.0,0.0,0.0);
      return pi;
    }
    const Vector ddpnn_dv1(dp*dnn_dv1+dp_dv1*nn);
    const Vector ddpnn_dv2(dp*dnn_dv2+dp_dv2*nn);

    const double x(-1.0/sqrt(1-dpnn*dpnn));

    d1=x*ddpnn_dv1;
    d2=x*ddpnn_dv2;


    return std::acos(dpnn);
}

