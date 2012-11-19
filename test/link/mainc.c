#include "../../src/wrapper/Plumed.h"
#include <stdlib.h>

int main(int argc,char** argv){
  plumed p;
  int natoms=3;

  if(plumed_installed()){
    plumed_gcreate();
    plumed_gcmd("setMDEngine","ACCode");
    plumed_gcmd("setNatoms",&natoms);
    plumed_gcmd("init",NULL);
    plumed_gcmd("read",NULL);
    plumed_gfinalize();

    p=plumed_create();
    plumed_cmd(p,"setNatoms",&natoms);
    plumed_cmd(p,"init",NULL);
    plumed_cmd(p,"read",NULL);
    plumed_finalize(p);
  }
  return 0;
}
