#include "plumed/wrapper/Plumed.h"
#include "plumed/tools/Exception.h"
#include <fstream>
#include <thread>
#include <functional>
#include <vector>
#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace PLMD;

void run(std::ostream & os, const std::string & name,std::function<void(int)> f,unsigned nthreads=4,unsigned nrepeats=10){
  os<<"Test "<<name<<" with OpenMP...\n";
  {
    for(unsigned i=0;i<nrepeats;i++) {
#if defined(_OPENMP)
#pragma omp parallel num_threads(nthreads)
      f(omp_get_thread_num());
#else
      f(0);
#endif
    }
  }
  os<<"OK"<<std::endl;
  os<<"Test "<<name<<" with C++11 threads...\n";
  for(unsigned i=0;i<nrepeats;i++) {
    std::vector<std::thread> threads;
    for(unsigned j=0;j<nthreads;j++) threads.emplace_back(std::thread(f,j));
    for(unsigned j=0;j<nthreads;j++) threads[j].join();
  }
  os<<"OK"<<std::endl;
}

int main(){
  try { // make sure "log" file is flushed

  std::ofstream ofs("log");

  run(ofs,"constructor and destructor",[&](int n){for(unsigned j=0;j<10000;j++) PLMD::Plumed q;});

  PLMD::Plumed p;
  run(ofs,"reference counter",[&](int n){for(unsigned j=0;j<100000;j++) PLMD::Plumed q=p;});

  run(ofs,"driver",[&](int n){
    char buffer[1024];
    { // remove output, if present
      std::sprintf(buffer,"test%d.xyz",n);
      std::remove(buffer);
    }
    { // write plumed file
      std::sprintf(buffer,"plumed%d.dat",n);
      auto fp=std::fopen(buffer,"w");
      std::fprintf(fp,"DUMPATOMS ATOMS=1 FILE=test%d.xyz PRECISION=2\n",n);
      std::fclose(fp);
    }
    { // run
      std::sprintf(buffer,"plumed driver --mf_xtc traj.xtc --plumed plumed%d.dat",n);
      PLMD::Plumed q;
      q.cmd("CLTool setArgvLine",buffer);
      int err;
      q.cmd("CLTool run",&err);
      plumed_assert(err==0);
    }
    { // compare with test_reference.xyz
      std::sprintf(buffer,"test%d.xyz",n);
      std::fstream file1(buffer);
      std::fstream file2("test_reference.xyz");
      while(!file2.eof())
      {
        char string1[1024];
        char string2[1024];
        file1.getline(string1,1024);
        file2.getline(string2,1024);
        if(std::strcmp(string1,string2)) {
          plumed_error()<<"Different lines:\n"<<string1<<"\n"<<string2<<"\n";
        }
      }
    }
  },4,10);

  } catch(...) {
    throw;
  }

  return 0;
}

