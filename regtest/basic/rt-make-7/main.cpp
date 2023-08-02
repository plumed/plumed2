#include "plumed/wrapper/Plumed.h"
#include "plumed/tools/Communicator.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdio>

using namespace PLMD;

class A{
  const int a;
public:
  operator int(){return a;};
  A(const int a):
    a(a)
  {}
};

// short test to see if MPI surrogate function do what they should
void small_test_mpi() {
  Communicator comm;
  Vector x0(1,2,3),x=x0;
  Vector y0(4,5,6),y=y0;
  
  plumed_assert(x[0]==x0[0] && x[1]==x0[1] && x[2]==x0[2]);
  plumed_assert(y[0]==y0[0] && y[1]==y0[1] && y[2]==y0[2]);
  comm.Sum(x);
  plumed_assert(x[0]==x0[0] && x[1]==x0[1] && x[2]==x0[2]);
  plumed_assert(y[0]==y0[0] && y[1]==y0[1] && y[2]==y0[2]);
  comm.Allgather(x,y);
  plumed_assert(x[0]==x0[0] && x[1]==x0[1] && x[2]==x0[2]);
  plumed_assert(y[0]==x0[0] && y[1]==x0[1] && y[2]==x0[2]);
  std::vector<int> count(1,3);
  std::vector<int> displ(1,0);
  x=x0;
  y=y0;
  comm.Allgatherv(y,x,count.data(),displ.data());
  plumed_assert(x[0]==y0[0] && x[1]==y0[1] && x[2]==y0[2]);
  plumed_assert(y[0]==y0[0] && y[1]==y0[1] && y[2]==y0[2]);
}

template<typename T>
void test_convert() {
  Plumed plumed;
  int size=sizeof(T);
  T t=0.0;
  plumed.cmd("setRealPrecision",size);
  plumed.cmd("convert cos(0.0)",&t);
  plumed_assert(t==1.0);
}

void test_checkAction() {
  Plumed plumed;
  int i=10;
  plumed.cmd("checkAction DISTANCE",&i);
  plumed_assert(i==1);
  i=10;
  plumed.cmd("checkAction ANY_NONEXISTING_NAME",&i);
  plumed_assert(i==0);
}

void test_cl() {
  {
    Plumed plumed;
    int ret;
    plumed.cmd("CLTool setArgvLine","plumed --help");
    ret=10;
    plumed.cmd("CLTool run",&ret);
    plumed_assert(ret==0);
  }
  {
    Plumed plumed;
    auto fp=std::fopen("tmp_out","w");
    int ret;
    plumed.cmd("CLTool setArgvLine","plumed --help");
    plumed.cmd("CLTool setOut",fp);
    ret=10;
    plumed.cmd("CLTool run",&ret);
    plumed_assert(ret==0);
    std::fclose(fp);

    std::ifstream check("tmp_out");
    std::string line;
    auto found=false;
    while(std::getline(check,line)) {
      if(line=="Commands:") found=true;
    }
    plumed_assert(found);
  }
}

void test_xyz() {
  Plumed p;
  auto natoms=4;
  std::vector<double> posx(natoms),posy(natoms),posz(natoms);
  std::vector<double> forx(natoms,0.0),fory(natoms,0.0),forz(natoms,0.0);
  std::vector<double> masses(natoms);
  for(auto i=0; i<posx.size(); i++) masses[i]=i+1;
  for(auto i=0; i<posx.size(); i++) posx[i]=10*i+0;
  for(auto i=0; i<posy.size(); i++) posy[i]=10*i+1;
  for(auto i=0; i<posz.size(); i++) posz[i]=10*i+2;
  double cell[9];
  double virial[9];
  for(auto i=0; i<9; i++) cell[i]=0.0;
  for(auto i=0; i<9; i++) virial[i]=0.0;
  cell[0]=100.0;
  cell[4]=100.0;
  cell[8]=100.0;
  int stopflag=0;
  A a(natoms);
  p.cmd("setNatoms",a);
  p.cmd("init",nullptr); // Test this: https://github.com/plumed/plumed2/issues/705
  try {
    p.cmd("setStopFlag",stopflag);
    plumed_error()<<"should have crashed";
  } catch (PLMD::Plumed::ExceptionTypeError & e) {
  }
  p.cmd("setStopFlag",&stopflag);
  p.cmd("readInputLine","DUMPATOMS ATOMS=@mdatoms FILE=test_xyz.xyz");
  p.cmd("readInputLine","c: COM ATOMS=@mdatoms");
  p.cmd("readInputLine","p: POSITION ATOM=c");
  p.cmd("readInputLine","RESTRAINT ARG=p.x,p.y,p.z AT=0.0,0.0,0.0 KAPPA=0.0,0.0,0.0 SLOPE=1.0,2.0,3.0");
  p.cmd("setBox",cell);
  p.cmd("setStep",0);
  p.cmd("setVirial",virial);
  p.cmd("setMasses",masses.data());
  p.cmd("setPositionsX",posx.data());
  p.cmd("setPositionsY",posy.data());
  p.cmd("setPositionsZ",posz.data());
  p.cmd("setForcesX",forx.data());
  p.cmd("setForcesY",fory.data());
  p.cmd("setForcesZ",forz.data());
  p.cmd("calc");
  std::ofstream ofs("test_xyz.forces");
  for(auto i=0; i<natoms; i++) ofs<<forx[i]<<" "<<fory[i]<<" "<<forz[i]<<"\n";
}

int main(){

  small_test_mpi();

  test_convert<double>();
  test_convert<float>();

  test_cl();

  Plumed* plumed=new Plumed;

  unsigned natoms=10;

  std::vector<double> positions(3*natoms,0.0);
  for(unsigned i=0;i<natoms;i++) positions[i]=i/10.0;
  std::vector<double> masses(natoms,1.0);
  std::vector<double> forces(3*natoms,0.0);
  std::vector<double> box(9,0.0);
  std::vector<double> virial(9,0.0);

  try {
    double dnatoms=natoms;
    plumed->cmd("setNatoms",dnatoms);
    plumed_error() << "should have failed with a typecheck error";
  } catch(PLMD::Plumed::ExceptionTypeError & e) {
  }
  plumed->cmd("setNatoms",natoms);
  plumed->cmd("setLogFile","test.log");
  plumed->cmd("init");
  plumed->cmd("readInputLine","UNITS LENGTH=A");

  plumed->cmd("readInputLine","d: TORSION ATOMS=1,1,1,1");
  plumed->cmd("clear"); // this is to test the clear command

  plumed->cmd("readInputLine","d: DISTANCE ATOMS=1,2");
  plumed->cmd("readInputLine","d1: DISTANCE ATOMS={1 2}"); // check if braces are parsed correctly
  plumed->cmd("readInputLine","PRINT ARG=d,d1 FILE=COLVAR");
  plumed->cmd("readInputLines","RESTRAINT ARG=d AT=0 KAPPA=1\n"
                               "METAD ARG=d PACE=1 SIGMA=1 HEIGHT=0 FILE=H1\n"
                               "METAD ARG=d PACE=2 SIGMA=1 HEIGHT=0 FILE=H2"
             );

  std::ofstream ofs("output");

  for(int step=0;step<10;step++){
    plumed->cmd("setStep",&step);
    plumed->cmd("setPositions",&positions[0],{natoms,3});
    plumed->cmd("setBox",&box[0],{3,3});
    plumed->cmd("setForces",&forces[0],forces.size());
    plumed->cmd("setVirial",&virial[0],9);
    plumed->cmd("setMasses",&masses[0],masses.size());
// first compute using modified positions:
    positions[0]=0.05;
    plumed->cmd("prepareCalc");
    plumed->cmd("performCalcNoUpdate");
    positions[0]=0;
    double bias=0;
    plumed->cmd("getBias",&bias,1);
    ofs<<bias<<"\n";
// first compute using regular positions:
    plumed->cmd("prepareCalc");
    plumed->cmd("performCalcNoUpdate");
    plumed->cmd("getBias",&bias,1);
    ofs<<bias<<"\n";
// hills should only be added at regular positions:
    plumed->cmd("update");
  }

  test_checkAction();

  delete plumed;

  test_xyz();

  return 0;
}
