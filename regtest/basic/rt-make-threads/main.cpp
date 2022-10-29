#include "plumed/wrapper/Plumed.h"
#include <fstream>
#include <thread>
#include <functional>
#include <vector>

using namespace PLMD;

void run(std::ostream & os, const std::string & name,std::function<void(void)> f,unsigned nthreads=4,unsigned nrepeats=10){
  os<<"Test "<<name<<" with OpenMP...\n";
  {
    for(unsigned i=0;i<nrepeats;i++) {
#pragma omp parallel num_threads(nthreads)
      f();
    }
  }
  os<<"OK\n";
  os<<"Test "<<name<<" with C++11 threads...\n";
  for(unsigned i=0;i<nrepeats;i++) {
    std::vector<std::thread> threads;
    for(unsigned j=0;j<nthreads;j++) threads.emplace_back(std::thread(f));
    for(unsigned j=0;j<nthreads;j++) threads[j].join();
  }
  os<<"OK\n";
}

int main(){

  std::ofstream ofs("log");

  run(ofs,"constructor and destructor",[&](){for(unsigned j=0;j<10000;j++) PLMD::Plumed q;});

  PLMD::Plumed p;
  run(ofs,"reference counter",[&](){for(unsigned j=0;j<100000;j++) PLMD::Plumed q=p;});

  return 0;
}

