
#include "OpenMP.h"
#include "Tools.h"
#include <cstdlib>

namespace PLMD{


void OpenMP::init(){
  if(std::getenv("PLUMED_NUM_THREADS")) Tools::convert(std::getenv("PLUMED_NUM_THREADS"),numThreads);
  if(std::getenv("PLUMED_CACHELINE_SIZE")) Tools::convert(std::getenv("PLUMED_CACHELINE_SIZE"),cachelineSize);
  initialized=true;
}

bool OpenMP::initialized=false;

unsigned OpenMP::numThreads=1;

unsigned OpenMP::cachelineSize=512;

}



