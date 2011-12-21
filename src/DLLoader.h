#ifndef __PLUMED_DLLoader_h
#define __PLUMED_DLLoader_h

#include <stack>
#include <string>

namespace PLMD{

class DLLoader{
  std::stack<void*> handles;
  std::string lastError;
public:
  ~DLLoader();
  void* load(const std::string&);
  const std::string & error();
  static bool installed();
};

}

#endif
