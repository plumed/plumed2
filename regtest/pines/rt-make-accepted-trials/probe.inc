// Mode 0 local control, 1 same-step foreign trial, 2 fresh foreign control,
// 3 extra prepareDependencies before the ordinary prepareCalc.
#define __PLUMED_WRAPPER_LIBCXX11 1
#define __PLUMED_WRAPPER_LIBCXX17 1
#include "plumed/wrapper/Plumed.h"
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
int runProbe(int argc,char**argv) {
  if(argc!=3) {
    return 2;
  }
  std::string action=argv[1];
  int mode=std::stoi(argv[2]);
  PLMD::Plumed p;
  p.cmd("setNatoms",3);
  p.cmd("setTimestep",.001);
  p.cmd("setKbT",2.5);
  p.cmd("setLogFile","plumed.log");
  p.cmd("init");
  std::ofstream ref("ref.pdb");
  ref<<"CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1           1\n"
     <<"ATOM      1  AR  ARG A   1       0.000   0.000   0.000  1.00  0.00\n"
     <<"ATOM      2  AR  ARG A   2      10.000   0.000   0.000  1.00  0.00\n"
     <<"ATOM      3  AR  ARG A   3      30.000   0.000   0.000  1.00  0.00\nEND\n";
  ref.close();
  if(action=="piv") {
    p.cmd("readInputLine","c: PIV PRECISION=100000 REF_FILE=ref.pdb PIVATOMS=1 ATOMTYPES=AR ONLYDIRECT SFACTOR=1 SORT=1 SWITCH1={RATIONAL R_0=0.7 MM=12 NN=4}");
    p.cmd("readInputLine","a: RESTRAINT ARG=c AT=0 KAPPA=1");
  } else if(action=="pines") {
    p.cmd("readInputLine","c: PINES REF_FILE=ref.pdb N_BLOCKS=1 BLOCK1={G1={ATOMID=0} G2={ATOMID={1 2}} SIZE=1 BUFFER=0 PL_REFRESH=1} SWITCH1={RATIONAL R_0=0.7 MM=12 NN=4}");
    p.cmd("readInputLine","a: RESTRAINT ARG=c.ELEMENT-0 AT=0 KAPPA=1");
  } else if(action=="emmi") {
    std::ofstream atoms("atoms.pdb");
    atoms<<"ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00\n"
         <<"ATOM      2  CA  ALA A   2      20.000   0.000   0.000  1.00  0.00\n"
         <<"ATOM      3  CA  ALA A   3      40.000   0.000   0.000  1.00  0.00\nEND\n";
    atoms.close();
    std::ofstream gmm("gmm.dat");
    gmm<<"#! FIELDS Id Weight Mean_0 Mean_1 Mean_2 Cov_00 Cov_01 Cov_02 Cov_11 Cov_12 Cov_22 Beta\n"
       <<"0 0.5 1.8 0 0 0.01 0 0 0.01 0 0.01 0\n"
       <<"1 0.5 3.8 0 0 0.01 0 0 0.01 0 0.01 0\n";
    gmm.close();
    p.cmd("readInputLine","MOLINFO STRUCTURE=atoms.pdb");
    p.cmd("readInputLine","c: EMMI NO_AVER NOPBC TEMP=300 NL_STRIDE=10 NL_CUTOFF=0.1 ATOMS=1-3 GMM_FILE=gmm.dat SIGMA_MIN=0.01 RESOLUTION=0.1 NOISETYPE=MARGINAL");
    p.cmd("readInputLine","a: BIASVALUE ARG=c.scoreb");
  } else if(action=="funcpath") {
    std::ofstream path("path.dat");
    path<<"0 0.5\n1 2\n2 4\n";
    path.close();
    p.cmd("readInputLine","d: DISTANCE ATOMS=1,2 NOPBC");
    p.cmd("readInputLine","c: FUNCPATHGENERAL ARG=d COEFFICIENTS=1 REFERENCE=path.dat LAMBDA=2 NEIGH_SIZE=1 NEIGH_STRIDE=0.002");
    p.cmd("readInputLine","a: RESTRAINT ARG=c.z AT=0 KAPPA=1");
  } else {
    return 3;
  }
  std::array<double,9>x{{0,0,0,2,0,0,4,0,0}},f{},v{},box{{10,0,0,0,10,0,0,0,10}};
  std::array<double,3>mass{{1,1,1}};
  auto bind=[&](int step,bool foreign) {
    x[3]=foreign?5.0:(step<3?3.99:4.01);
    f.fill(0);
    v.fill(0);
    p.cmd("setStep",step);
    p.cmd("setPositions",x.data());
    p.cmd("setMasses",mass.data());
    p.cmd("setBox",box.data());
    p.cmd("setForces",f.data());
    p.cmd("setVirial",v.data());
  };
  std::cout<<std::setprecision(17);
  for(int step=0; step<5; step++) {
    bind(step,mode==2 || ((mode==5 || mode==7) && step>2));
    if(mode==3) {
      p.cmd("prepareDependencies");
    }
    p.cmd("calc");
    double bias=0;
    p.cmd("getBias",&bias);
    std::cout<<"normal "<<step<<" "<<bias<<" "<<f[3]<<"\n";
    if((mode==1 || mode==5) && step==2) {
      bind(step,true);
      p.cmd("setExchangeStep",1);
      p.cmd("prepareCalc");
      p.cmd("performCalcNoForces");
      p.cmd("getBias",&bias);
      std::cout<<"trial "<<step<<" "<<bias<<"\n";
      p.cmd("setExchangeStep",0);
    }
  }
  return 0;
}
