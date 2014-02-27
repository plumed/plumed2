
#include "OpenMP.h"
#include "Tools.h"
#include <cstdlib>

namespace PLMD{

unsigned OpenMP::getCachelineSize(){
  static unsigned cachelineSize=512;
  if(std::getenv("PLUMED_CACHELINE_SIZE")) Tools::convert(std::getenv("PLUMED_CACHELINE_SIZE"),cachelineSize);
  return cachelineSize;
}

unsigned OpenMP::getNumThreads(){
  static unsigned numThreads=1;
  if(std::getenv("PLUMED_NUM_THREADS")) Tools::convert(std::getenv("PLUMED_NUM_THREADS"),numThreads);
  return numThreads;
}


}



