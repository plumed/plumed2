#include "Pbc.h"
#include "Tools.h"
#include "PlumedException.h"

using namespace PLMD;

Pbc::Pbc():
  type(unset)
{
  box.clear();
  invBox.clear();
}


void Pbc::setBox(const Tensor&b){
  box=b;
// UP TO NOW ONLY WORKS WITH ORTHOROMIBIC (should implement matrix inversion)
// detect type:
  const double epsilon=1e-14;

  if(box.determinant()<epsilon)return;

  bool cxy=false;
  bool cxz=false;
  bool cyz=false;
  if(box(0,1)*box(0,1)<epsilon && box(1,0)*box(1,0)<epsilon) cxy=true;
  if(box(0,2)*box(0,2)<epsilon && box(2,0)*box(2,0)<epsilon) cxz=true;
  if(box(1,2)*box(1,2)<epsilon && box(2,1)*box(2,1)<epsilon) cyz=true;

  invBox=box.inverse();

  if(cxy && cxz && cyz) type=orthorombic;
// NOT IMPLEMENTED YET; WILL FALL TO GENERIC
//  else if(cxy && cxz)   type=yz;
//  else if(cxy && cyz)   type=xz;
//  else if(cxz && cyz)   type=xy;
  else type=generic;
}

Vector Pbc::distance(const Vector&v1,const Vector&v2)const{
  Vector d=delta(v1,v2);
  if(type==unset){
  } else if(type==orthorombic) {
    for(int i=0;i<3;i++) d[i]=Tools::pbc(d[i]*invBox(i,i))*box(i,i);
  } else if(type==generic) {
// first reduce:
   Vector s=matmul(invBox.transpose(),d);
   for(int i=0;i<3;i++) s[i]=Tools::pbc(s[i]);
   d=matmul(box.transpose(),s);
// then check all the neighbors:
    Vector a1(box(0,0),box(0,1),box(0,2));
    Vector a2(box(1,0),box(1,1),box(1,2));
    Vector a3(box(2,0),box(2,1),box(2,2));
    Vector best(d);
    double lbest=best.modulo2();
    Vector trial;
    double ltrial;
    for(int i=-1;i<=1;i++) for(int j=-1;j<=1;j++) for(int k=-1;k<=1;k++) {
      trial=d+i*a1+j*a2+k*a3;
      ltrial=trial.modulo2();
      if(ltrial<lbest){
        best=trial;
        lbest=ltrial;
      }
    }
    d=best;
  } else plumed_merror("unknown pbc type");
  return d;
}

Vector Pbc::realToScaled(const Vector&d)const{
  return matmul(invBox.transpose(),d);
}

Vector Pbc::scaledToReal(const Vector&d)const{
  return matmul(box.transpose(),d);
}

bool Pbc::isOrthorombic()const{
  return type==orthorombic;
}

const Tensor& Pbc::getBox()const{
  return box;
}


