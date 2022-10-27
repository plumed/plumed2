#include "plumed/wrapper/Plumed.h"
#include <fstream>

using namespace PLMD;

int main(){

  std::ofstream ofs("log");

// test if constructor/destructor are thread safe
  {
    ofs<<"Testing if constructor and destructor are thread safe ..."<<std::endl;
    unsigned nouter=10;
    unsigned ninner=10000;
    for(unsigned i=0;i<nouter;i++) {
#pragma omp parallel
      for(unsigned j=0;j<ninner;j++) PLMD::Plumed p;
    }
    ofs<<"OK"<<std::endl;
  }

  return 0;
}

