#ifdef __PLUMED_LIBCXX11
// In order to correctly catch the thrown C++11 exceptions,
// we notify the Plumed wrapper that those exceptions are recognized by the compiler.
#define __PLUMED_WRAPPER_LIBCXX11 1
#endif
#include "plumed/wrapper/Plumed.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <unistd.h>
#include <sys/stat.h>

bool exists(const std::string & path) {
  struct stat st;
  return stat(path.c_str(),&st)==0;
}

// Check that the same file cannot be opened for writing twice
// by actions of the same plumed object (see issue #802)

using namespace PLMD;

void test_line(std::ostream & ofs,Plumed & p,const std::string & arg,bool shouldFail,const std::string & shown="") {
  ofs<<"readInputLine "<<(shown.empty()?arg:shown)<<std::endl;
  try {
    p.cmd("readInputLine",arg.c_str());
    ofs<<(shouldFail?"+++ !!!! uncatched !!!!":"+++ ok")<<std::endl;
  } catch(Plumed::Exception&e) {
    std::string msg(e.what());
    ofs<<(shouldFail?"+++ catched":"+++ !!!! unexpected exception !!!!")<<std::endl;
    if(msg.find("is already open for writing")!=std::string::npos) ofs<<"+++ message ok"<<std::endl;
  }
}

int main() {

  std::ofstream ofs("output");

  std::vector<char> cwd(4096);
  if(!getcwd(cwd.data(),cwd.size())) return 1;
  std::string abspath=std::string(cwd.data())+"/COLVAR";

  {
    Plumed plumed;
    int natoms=10;
    plumed.cmd("setNatoms",&natoms);
    plumed.cmd("setLogFile","test.log");
    plumed.cmd("init");

    test_line(ofs,plumed,"d: DISTANCE ATOMS=1,2",false);
    test_line(ofs,plumed,"p1: PRINT ARG=d FILE=COLVAR",false);
// same name
    test_line(ofs,plumed,"p2: PRINT ARG=d FILE=COLVAR",true);
// same file, different name
    test_line(ofs,plumed,"p3: PRINT ARG=d FILE=./COLVAR",true);
    test_line(ofs,plumed,"p4: PRINT ARG=d FILE="+abspath,true,"p4: PRINT ARG=d FILE=/absolute/path/to/COLVAR");
// same file in append mode
    test_line(ofs,plumed,"p5: PRINT ARG=d FILE=COLVAR RESTART=YES",true);
// a different file is fine
    test_line(ofs,plumed,"p6: PRINT ARG=d FILE=COLVAR2",false);
// /dev/null can be shared
    test_line(ofs,plumed,"p7: PRINT ARG=d FILE=/dev/null",false);
    test_line(ofs,plumed,"p8: PRINT ARG=d FILE=/dev/null",false);
// two METAD with the default FILE=HILLS
    test_line(ofs,plumed,"m1: METAD ARG=d SIGMA=0.1 HEIGHT=1.0 PACE=2",false);
    test_line(ofs,plumed,"m2: METAD ARG=d SIGMA=0.1 HEIGHT=1.0 PACE=2",true);
    test_line(ofs,plumed,"m3: METAD ARG=d SIGMA=0.1 HEIGHT=1.0 PACE=2 FILE=HILLS2",false);
// the check is done before backing up, so the files written by p1 and m1 should not have been moved
    ofs<<"bck.0.COLVAR exists: "<<exists("bck.0.COLVAR")<<std::endl;
    ofs<<"bck.0.HILLS exists: "<<exists("bck.0.HILLS")<<std::endl;
  }

  {
// a different plumed object is not affected by files opened by the previous one,
// which are now closed
    Plumed plumed;
    int natoms=10;
    plumed.cmd("setNatoms",&natoms);
    plumed.cmd("setLogFile","test2.log");
    plumed.cmd("init");
    test_line(ofs,plumed,"d: DISTANCE ATOMS=1,2",false);
    test_line(ofs,plumed,"p1: PRINT ARG=d FILE=COLVAR",false);
  }

  return 0;
}
