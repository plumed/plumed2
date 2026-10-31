#ifdef __PLUMED_LIBCXX11
// In order to correctly catch the thrown C++11 exceptions,
// we notify the Plumed wrapper that those exceptions are recognized by the compiler.
#define __PLUMED_WRAPPER_LIBCXX11 1
#endif
#include "mpi.h"
#include "plumed/wrapper/Plumed.h"
#include <fstream>
#include <string>
#include <vector>
#include <unistd.h>

// With MPI, only the first process opens the output files (the others open /dev/null).
// The check that the same file is not opened twice (issue #802) must still stop
// all the processes, also when the same file is reached through a different path,
// otherwise the processes that did not stop would hang in the next collective operation.

using namespace PLMD;

void test_line(std::ostream & ofs,Plumed & p,const std::string & arg,bool shouldFail,const std::string & shown="") {
  ofs<<"readInputLine "<<(shown.empty()?arg:shown)<<std::endl;
  try {
    p.cmd("readInputLine",arg.c_str());
    ofs<<(shouldFail?"+++ !!!! uncatched !!!!":"+++ ok")<<std::endl;
  } catch(Plumed::Exception&e) {
    std::string msg(e.what());
    ofs<<(shouldFail?"+++ catched":"+++ !!!! unexpected exception !!!!")<<std::endl;
    if(msg.find("is already open for writing")!=std::string::npos) {
      ofs<<"+++ message ok"<<std::endl;
    }
  }
}

int main(int argc,char**argv) {
  MPI_Init(&argc,&argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  std::ofstream ofs("output"+std::to_string(rank));

  std::vector<char> cwd(4096);
  if(!getcwd(cwd.data(),cwd.size())) {
    return 1;
  }
  const std::string abspath=std::string(cwd.data())+"/COLVAR";
  {
    Plumed plumed;
    MPI_Comm comm=MPI_COMM_WORLD;
    plumed.cmd("setMPIComm",&comm);
    int natoms=10;
    plumed.cmd("setNatoms",&natoms);
    plumed.cmd("setLogFile","test.log");
    plumed.cmd("init");

    test_line(ofs,plumed,"d: DISTANCE ATOMS=1,2",false);
    test_line(ofs,plumed,"p1: PRINT ARG=d FILE=COLVAR",false);
    test_line(ofs,plumed,"p2: PRINT ARG=d FILE=COLVAR",true);
    test_line(ofs,plumed,"p3: PRINT ARG=d FILE=./COLVAR",true);
    test_line(ofs,plumed,"p4: PRINT ARG=d FILE="+abspath,true,"p4: PRINT ARG=d FILE=/absolute/path/to/COLVAR");
    test_line(ofs,plumed,"p5: PRINT ARG=d FILE=COLVAR2",false);

    // a step, to check that the files are written by the first process only once
    std::vector<double> pos(3*natoms,0.0), forces(3*natoms,0.0), masses(natoms,1.0);
    for(int i=0; i<3*natoms; i++) {
      pos[i]=0.1*i;
    }
    double box[9]= {10,0,0,0,10,0,0,0,10};
    double virial[9]= {0};
    int step=0;
    plumed.cmd("setStep",&step);
    plumed.cmd("setPositions",pos.data(),3*natoms);
    plumed.cmd("setForces",forces.data(),3*natoms);
    plumed.cmd("setMasses",masses.data(),natoms);
    plumed.cmd("setBox",box,9);
    plumed.cmd("setVirial",virial,9);
    plumed.cmd("calc");
  }
  ofs<<"done"<<std::endl;
  MPI_Finalize();
  return 0;
}
