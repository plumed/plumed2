#ifndef __PLUMED_ExchangePatterns_h
#define __PLUMED_ExchangePatterns_h

#include "Random.h"

namespace PLMD {
  class ExchangePatterns {
    Random random;
public:
  void setSeed(int seed);
  void getList(int *ind, int nrepl);
};
}
#endif
