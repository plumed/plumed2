#include "plumed/wrapper/Plumed.h"
#include <stdlib.h>
#include <stdio.h>

int main(int argc,char*argv[]) {
  plumed p;
  double* positions;
  double* masses;
  double* forces;
  double box[3][3];
  double virial[3][3];
  int ene,stopflag,itmp;
  double bias;
  unsigned i,j;
  FILE* log;

  plumed_gcreate();
  plumed_gcmd("setNatoms",1*10);
  plumed_gcmd("init",0);
  plumed_gfinalize();

  p=plumed_create();

  positions=(double*)malloc(3*10*sizeof(double));
  forces=(double*)malloc(3*10*sizeof(double));
  masses=(double*)malloc(10*sizeof(double));
  for(i=0;i<10;i++) {
    positions[3*i+0]=1.0+10*i;
    positions[3*i+1]=2.0+10*i;
    positions[3*i+2]=3.0+10*i;
  }

  for(i=0;i<3*10;i++) forces[i]=0.0;
  for(i=0;i<10;i++) masses[i]=1.0;
  for(i=0;i<3;i++) for(j=0;j<3;j++) box[i][j]=0.0;
  for(i=0;i<3;i++) for(j=0;j<3;j++) virial[i][j]=0.0;

  stopflag=0;
  ene=1;
  bias=-10;

  plumed_cmd(p,"setNatoms",1*10);;

  log=fopen("testfile","w");

  plumed_cmd(p,"init");
  plumed_cmd(p,"setStopFlag",&stopflag);
  plumed_cmd(p,"readInputLine","p: POSITION ATOM=2");
  plumed_cmd(p,"readInputLine","g: GYRATION ATOMS=@allatoms");
  plumed_cmd(p,"readInputLine","r: RESTRAINT ARG=g AT=0 KAPPA=3");
  plumed_cmd(p,"readInputLine","COMMITTOR ARG=p.x STRIDE=1 BASIN_LL1=0 BASIN_UL1=30");
  plumed_cmd(p,"readInputLine","PRINT ARG=p.*,r.* FILE=testme2");
  itmp=1;
  plumed_cmd(p,"setStep",&itmp);
  plumed_cmd(p,"setPositions",positions,30);
  plumed_cmd(p,"setForces",forces,30);
  plumed_cmd(p,"setMasses",masses,10);
  plumed_cmd(p,"setBox",&box[0][0],9);
  plumed_cmd(p,"setVirial",&virial[0][0],9);
  fprintf(log,"stopflag should be 0: %d\n",stopflag);
  fprintf(log,"isEnergyNeeded should be 1: %d\n",ene);
  plumed_cmd(p,"isEnergyNeeded",&ene);
  fprintf(log,"isEnergyNeeded should be 0: %d\n",ene);
  fprintf(log,"stopflag should be 0: %d\n",stopflag);
  plumed_cmd(p,"calc");
  fprintf(log,"stopflag should be 1: %d\n",stopflag);
  plumed_cmd(p,"getBias",&bias);
  fprintf(log,"bias: %lf\n",bias);

  fclose(log);
  plumed_finalize(p);

  free(positions);
  free(forces);
  free(masses);

  return 0;
}
