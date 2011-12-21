#ifndef __PLUMED_DLLoader_h
#define __PLUMED_DLLoader_h

#include <stack>
#include <string>

namespace PLMD{

/// Class taking care of dynamic loading.
/// It contains wrappers to the dlopen() routine.
/// It is designed so that when an object of this class goes
/// out of scope all the libraries loaded by it are unloaded. In this
/// manner, loaded libraries are automatically unloaded at the end of
/// execution. Libraries are loaded with RTDL_LOCAL option, which
/// means that they are not accessible from outside. Still, if they
/// contain self-registering classes, they will register themselves
/// to the ActionRegister object.
class DLLoader{
  std::stack<void*> handles;
  std::string lastError;
public:
/// Cleanup
  ~DLLoader();
/// Load a library, returning its handle
  void* load(const std::string&);
/// Returns the last error in dynamic loader
  const std::string & error();
/// Returs true if the dynamic loader is available (on some systems it may not).
  static bool installed();
};

}

#endif
