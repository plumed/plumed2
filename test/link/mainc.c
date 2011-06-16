#include "Plumed.h"
#include <stdlib.h>

int main(int argc,char** argv){
  int check;
  plumed p;
  int natoms=3;

  plumed_installed(&check);
  if(check==1){
    plumed_g_create();
    plumed_g_cmd("setMDEngine","ACCode");
    plumed_g_cmd("setNatoms",&natoms);
    plumed_g_cmd("init",NULL);
    plumed_g_cmd("read",NULL);
    plumed_g_finalize();

    p=plumed_create();
    plumed_cmd(p,"setNatoms",&natoms);
    plumed_cmd(p,"init",NULL);
    plumed_cmd(p,"read",NULL);
    plumed_finalize(p);
  }
  return 0;
}
