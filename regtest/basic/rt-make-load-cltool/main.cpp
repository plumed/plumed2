#include "plumed/wrapper/Plumed.h"
#include <vector>
#include <thread>
#include <cstdio>

using namespace PLMD;

// i is thread number (starting from 0)
void fun(Plumed p,unsigned i) {
  char buffer[1000];
  std::snprintf(buffer,1000,"output_threads%d",i);
  auto fp=std::fopen(buffer,"w");
  int err=0;
  try {
    const char* argv[10];
    int argc=0;
    argv[argc++]="plumed";
    
    // we load different versions in different threads
    // all of them override the standard kt cltool
    if(i%3==1) {
      argv[argc++]="--load";
      argv[argc++]="./kT10." __PLUMED_SOEXT;
    } else if(i%3==2) {
      argv[argc++]="--load";
      argv[argc++]="./kT20." __PLUMED_SOEXT;
    }

    argv[argc++]="kt";
    argv[argc++]="--temp";
    argv[argc++]="300.68090634510198801675"; // enough digits to get 2.5 in the output
    p.cmd("CLTool setArgc",&argc);
    p.cmd("CLTool setArgv",argv);
    p.cmd("CLTool setOut",fp);
    p.cmd("CLTool run",&err);
    std::fclose(fp);
  } catch(...) {
    std::fclose(fp);
    throw;
  }
  if(err!=0) throw std::runtime_error("an error happened");
}

int main() {
  // threads
  {
    unsigned nthreads=16;
    std::vector<std::thread> threads;
    // plumed objects are declared outside
    std::vector<Plumed> pl(nthreads);
    for(unsigned j=0;j<nthreads;j++) threads.emplace_back(fun,pl[j],j);
    // threads are joined before destruction
    for(unsigned j=0;j<nthreads;j++) threads[j].join();

    // this is done to make sure that all plumed objects survive till the end
    // and make the number of unique so/dylib loaded reproducible
  }

  return 0;
}

