#include "plumed/tools/Pbc.h"
#include "plumed/tools/Random.h"
#include "plumed/tools/Exception.h"
#include "plumed/tools/Stopwatch.h"
#include <iostream>
#include <fstream>

using namespace PLMD;

int run(Stopwatch & sw,int boxtype,double* av_nshifts=NULL){
  Random r;
  int failures=0;
  r.setSeed(-20);
  int nshifts=0;
  int nbox=128;
  int nvec=512;
  for(int i=0;i<nbox;i++){
// random matrix with some zero element
  Tensor box;
  for(int j=0;j<3;j++) for(int k=0;k<3;k++) if(r.U01()>0.2){
    box[j][k]=2.0*r.U01()-1.0;
  }
  switch(boxtype){
    case 0:
// cubic
      for(int j=0;j<3;j++) for(int k=0;k<3;k++) if(j!=k) box[j][k]=0.0;
      for(int j=1;j<3;j++) box[j][j]=box[0][0];
      break;
    case 1:
// orthorombic
      for(int j=0;j<3;j++) for(int k=0;k<3;k++) if(j!=k) box[j][k]=0.0;
      break;
    case 2:
// hexagonal
      {
      int perm=r.U01()*100;
      Vector a;
      a(0)=r.U01()*2-2; a(1)=0.0;a(2)=0.0;
      double d=r.U01()*2-2;
      Vector b(0.0,d,0.0);
      Vector c(0.0,0.5*d,sqrt(3.0)*d*0.5);
      box.setRow((perm+0)%3,a);
      box.setRow((perm+1)%3,b);
      box.setRow((perm+2)%3,c);
      }
      break;
    case 3:
// bcc
      {
      int perm=r.U01()*100;
      double d=r.U01()*2-2;
      Vector a(d,d,d);
      Vector b(d,-d,d);
      Vector c(d,d,-d);
      box.setRow((perm+0)%3,a);
      box.setRow((perm+1)%3,b);
      box.setRow((perm+2)%3,c);
      }
      break;
    case 4:
// fcc
      {
      int perm=r.U01()*100;
      double d=r.U01()*2-2;
      Vector a(d,d,0);
      Vector b(d,0,d);
      Vector c(0,d,d);
      box.setRow((perm+0)%3,a);
      box.setRow((perm+1)%3,b);
      box.setRow((perm+2)%3,c);
      }
      break;
    default:
// triclinic
      break;
    }

    Pbc pbc;
    pbc.setBox(box);
    std::vector<Vector> vv(nvec);
    std::vector<Vector> result(nvec);
    std::vector<Vector> result_apply(nvec);
    std::vector<Vector> reference(nvec);
    for(int j=0;j<nvec;j++){
// random vector
      Vector v(r.U01()-0.5,r.U01()-0.5,r.U01()-0.5);
      v*=5;
// set some component to zero
      for(int j=0;j<3;j++) if(r.U01()>0.2) v(j)=0.0;
      vv[j]=v;
    }

// fast version
    sw.start("Fast");
    for(int j=0;j<nvec;j++) result[j]=pbc.distance(Vector(0,0,0),vv[j],&nshifts);
    sw.stop("Fast");

// apply version (in-place, vector operation)
    sw.start("Vector");
    result_apply=vv;
    pbc.apply(result_apply);
    sw.stop("Vector");

    reference=result;
// full search around that
    sw.start("Full");
    for(int j=0;j<nvec;j++) pbc.fullSearch(reference[j]);
    sw.stop("Full");

// compare
    for(int j=0;j<nvec;j++) if(std::fabs(modulo2(result[j])-modulo2(result_apply[j]))>1e-15) failures++;
    for(int j=0;j<nvec;j++) if(std::fabs(modulo2(result[j])-modulo2(reference[j]))>1e-15) failures++;
  }
  if(av_nshifts) *av_nshifts=double(nshifts)/double(nbox*nvec);
  return failures;
}


int main(){
  Stopwatch sw;
  sw.start();
  std::ofstream ofs("logfile");
  ofs<<std::fixed;
  for(unsigned type=0;type<6;type++){
    double nsh;
    int err=run(sw,type,&nsh);
    ofs<<"Box type "<<type<<"\n";
    ofs<<"Failures "<<err<<"\n";
    ofs.precision(1);
    ofs<<"Shifts   "<<nsh<<"\n\n";
  }
  sw.stop();
  std::cout<<sw;
  return 0;
}
