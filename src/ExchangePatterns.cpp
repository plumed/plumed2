#include "ExchangePatterns.h"

using namespace std;

namespace PLMD{

void ExchangePatterns::setFlag(const int flag){
  PatternFlag=flag;
}

void ExchangePatterns::getFlag(int &flag){
  PatternFlag=flag;
}

void ExchangePatterns::setSeed(int seed)
{
  random.setSeed(seed);
}

void ExchangePatterns::getList(int *ind, int nrepl)
{
  /* in principle here we can add a switch(patter) case in order to get a list of exchanges dependent on a specific pattern */
  switch(PatternFlag)
  {
    case RANDOM:
      for(int i=0;i<nrepl;i++) {
        int stat=1;
        while(stat) {
          stat=0;
          ind[i] = random.RandU01()*nrepl;
          for(int j=0;j<i;j++) if(ind[i]==ind[j]) stat=1;
        }
      }
      break;
    case NEIGHBOR:
      for(int i=0;i<nrepl;i++) ind[i]=i; 
      break; 
  }
}

}
