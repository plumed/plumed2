#ifndef __PLUMED_ExchangePatterns_h
#define __PLUMED_ExchangePatterns_h

#include "Random.h"

namespace PLMD {
  class ExchangePatterns {
    int    PatternFlag;
    Random random;
public:
  enum PatternFlags { NONE, RANDOM, NEIGHBOR, TOTAL };
  void setSeed(int seed);
  void getList(int *ind, int nrepl);
  void setFlag(const int);
  void getFlag(int&);
};
}
#endif
