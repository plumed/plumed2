#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <string.h>

int main(int argc,char*argv[]){
  unsigned nsym=0;
  for(int i=1;i<argc;i++) if(!strcmp(argv[i],"-h")){
    printf("Usage: wham -n NSYM\n");
    printf("Then pass a file with NSYM columns\n");
    return 0;
  }
  for(int i=1;i<argc;i++) if(!strcmp(argv[i],"-n")) sscanf(argv[i+1],"%u",&nsym);
  std::cerr<<"nsym "<<nsym<<"\n";

// array with exponentials of V
  std::vector<double> expv;
// array with partition functions
  std::vector<double> Z(nsym,1.0);
// array with inverse of partition functions
  std::vector<double> invZ(nsym,1.0);
  for(;;){
    double d;
    int x=std::scanf("%lf",&d);
    if(x==EOF) break;
    expv.push_back(d);
  }
  int nframes=expv.size()/nsym;
// array with weights:
  std::vector<double> weight(nframes,1.0);
  double mm=+1e10;
  for(unsigned i=0;i<nframes*nsym;i++) if(expv[i]<mm) mm=expv[i];
  for(unsigned i=0;i<nframes*nsym;i++) expv[i]-=mm;
  for(unsigned i=0;i<nframes*nsym;i++) expv[i]=std::exp(-expv[i]);

  for(unsigned it=0;it<1000;it++){

// store Z
    std::vector<double> Zold=Z;

// recompute weights
    for(unsigned i=0;i<nframes;i++){
      double ew=0.0;
      for(unsigned j=0;j<nsym;j++) ew+=expv[i*nsym+j]*invZ[j];
      weight[i]=1.0/ew;
    }

// normalize weights
    {
      double norm=0.0;
      for(unsigned i=0;i<nframes;i++) norm+=weight[i];
      norm=1.0/norm;
      for(unsigned i=0;i<nframes;i++) weight[i]*=norm;
    }
   
// recompute Z
    Z.assign(nsym,0.0);
    for(unsigned i=0;i<nframes;i++) for(unsigned j=0;j<nsym;j++) Z[j]+=weight[i]*expv[i*nsym+j];

// normalize Z
    {
      double norm=0.0;
      for(unsigned j=0;j<nsym;j++) norm+=Z[j];
      norm=1.0/norm;
      for(unsigned j=0;j<nsym;j++) Z[j]*=norm;
    }
// compute inverse Z
    for(unsigned j=0;j<nsym;j++) invZ[j]=1.0/Z[j];

// compute change in log Z
    double eps=0.0; for(unsigned i=0;i<Z.size();i++){
      double d=std::log(Z[i]/Zold[i]);
      eps+=d*d;
    }

// some log
    if(it%10==0){
      std::cerr<<"iteration "<<it<<"\n";
      std::cerr<<"Z: ";
      for(unsigned j=0;j<nsym;j++) std::cerr<<" "<<Z[j];
      std::cerr<<"\n";
      std::cerr<<"eps: "<<eps<<"\n";
    }

// possibly stop
    if(eps<1e-8)break;
  }
  for(unsigned i=0;i<nframes;i++) std::cout<<" "<<weight[i]<<"\n";
  return 0;
}
