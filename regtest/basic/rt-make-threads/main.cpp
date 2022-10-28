#include "plumed/wrapper/Plumed.h"
#include <fstream>
#include <thread>

using namespace PLMD;

int main(){

  std::ofstream ofs("log");

// test if constructor/destructor are thread safe
  {
    ofs<<"Testing if constructor and destructor are thread safe with OpenMP ..."<<std::endl;
    unsigned nouter=10;
    unsigned ninner=10000;
    for(unsigned i=0;i<nouter;i++) {
#pragma omp parallel
      for(unsigned j=0;j<ninner;j++) PLMD::Plumed p;
    }
    ofs<<"OK"<<std::endl;
  }

// test if reference counter is thread safe
  {
    ofs<<"Testing if reference counter is thread safe with OpenMP ..."<<std::endl;
    unsigned nouter=10;
    unsigned ninner=100000;
    PLMD::Plumed p;
    for(unsigned i=0;i<nouter;i++) {
#pragma omp parallel
      for(unsigned j=0;j<ninner;j++) PLMD::Plumed q=p;
    }
    ofs<<"OK"<<std::endl;
  }

/// test if constructor/destructor are thread safe
  {
    ofs<<"Testing if constructor and destructor are thread safe with C++11 threads ..."<<std::endl;
    unsigned nouter=10;
    unsigned ninner=10000;
    for(unsigned i=0;i<nouter;i++) {
      std::thread thread1([&](){for(unsigned j=0;j<ninner;j++) PLMD::Plumed p;});
      std::thread thread2([&](){for(unsigned j=0;j<ninner;j++) PLMD::Plumed p;});
      std::thread thread3([&](){for(unsigned j=0;j<ninner;j++) PLMD::Plumed p;});
      std::thread thread4([&](){for(unsigned j=0;j<ninner;j++) PLMD::Plumed p;});
      std::thread thread5([&](){for(unsigned j=0;j<ninner;j++) PLMD::Plumed p;});
      std::thread thread6([&](){for(unsigned j=0;j<ninner;j++) PLMD::Plumed p;});
      std::thread thread7([&](){for(unsigned j=0;j<ninner;j++) PLMD::Plumed p;});
      std::thread thread8([&](){for(unsigned j=0;j<ninner;j++) PLMD::Plumed p;});
      thread1.join();
      thread2.join();
      thread3.join();
      thread4.join();
      thread5.join();
      thread6.join();
      thread7.join();
      thread8.join();
    }
    ofs<<"OK"<<std::endl;
  }

// test if reference counter is thread safe
  {
    ofs<<"Testing if reference counter is thread safe with C++11 threads ..."<<std::endl;
    unsigned nouter=10;
    unsigned ninner=100000;
    PLMD::Plumed p;
    for(unsigned i=0;i<nouter;i++) {
      std::thread thread1([&](){for(unsigned j=0;j<ninner;j++) PLMD::Plumed q=p;});
      std::thread thread2([&](){for(unsigned j=0;j<ninner;j++) PLMD::Plumed q=p;});
      std::thread thread3([&](){for(unsigned j=0;j<ninner;j++) PLMD::Plumed q=p;});
      std::thread thread4([&](){for(unsigned j=0;j<ninner;j++) PLMD::Plumed q=p;});
      std::thread thread5([&](){for(unsigned j=0;j<ninner;j++) PLMD::Plumed q=p;});
      std::thread thread6([&](){for(unsigned j=0;j<ninner;j++) PLMD::Plumed q=p;});
      std::thread thread7([&](){for(unsigned j=0;j<ninner;j++) PLMD::Plumed q=p;});
      std::thread thread8([&](){for(unsigned j=0;j<ninner;j++) PLMD::Plumed q=p;});
      thread1.join();
      thread2.join();
      thread3.join();
      thread4.join();
      thread5.join();
      thread6.join();
      thread7.join();
      thread8.join();
    }
    ofs<<"OK"<<std::endl;
  }

  return 0;
}

