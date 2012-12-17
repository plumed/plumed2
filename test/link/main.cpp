#include "../../src/wrapper/Plumed.h"
#include <cstdio>


int ddoit(){
  fprintf(stderr,"Callback\n");
  return 10;
}

int main(int argc,char**argv){
  PLMD::Plumed plumed;
  bool inst=PLMD::Plumed::installed();
  const char* input="plumed.dat";

  if(argc>1 && argv[1]) input=argv[1];

  if(inst){

    int check;
    plumed.cmd("checkAction DISTANCE",&check);
    printf("Checking DISTANCE: %i\n",check);
    plumed.cmd("checkAction DISTANCE1",&check);
    printf("Checking DISTANCE1: %i\n",check);

//   {
//     plumed_function_holder xxx={(plumed_function_pointer)ddoit};
//     plumed.cmd("doit",&xxx);
//   }

    plumed.cmd("setMDEngine","specific");
    plumed.cmd("setLog",stdout);
    int n=20;
    plumed.cmd("setNatoms",&n);
    plumed.cmd("setPlumedDat",input);
    plumed.cmd("init");
  }
  return 0;
}
