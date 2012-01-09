#ifndef __PLUMED_Grex_h
#define __PLUMED_Grex_h

#include "WithCmd.h"
#include <string>
#include <vector>

namespace PLMD{

class PlumedMain;
class Atoms;
class PlumedCommunicator;

class GREX:
  public WithCmd
{
  bool initialized;
  PlumedCommunicator& intracomm;
  PlumedCommunicator& intercomm;
  PlumedMain& plumedMain;
  Atoms&      atoms;
  int partner;
  double localDeltaBias;
  double foreignDeltaBias;
  std::vector<double> allDeltaBias;
  std::string buffer;
  int myreplica;
public:
  GREX(PlumedMain&);
  ~GREX();
  void cmd(const std::string&key,void*val=NULL);
  void calculate();
  void savePositions();
};

}

#endif
