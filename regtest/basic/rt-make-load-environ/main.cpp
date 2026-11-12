#include "plumed/wrapper/Plumed.h"
#include <fstream>
#include <string>
#include <exception>

#ifdef __linux__
// Reproduce an executable's copy-relocated environment symbols. libc startup
// initializes the interposed __environ, while a deep-bound DSO can see libc's
// original, empty environ. The runtime wrapper enables RTLD_DEEPBIND by default.
extern "C" {
  char** __environ=nullptr;
  extern char** environ __attribute__((alias("__environ")));
}
#endif

int main() {
  std::ofstream result("result");
  try {
    PLMD::Plumed p;
    p.cmd("setLogFile","test.log");
    p.cmd("init");
    bool failed=false;
    try {
      // The mock compiler records its environment and deliberately fails.
      p.cmd("readInputLine","LOAD FILE=Probe.cpp");
    } catch(const std::exception&) {
      failed=true;
    }
    std::string inherited;
    std::getline(std::ifstream("child-environment"),inherited);
    if(failed && inherited=="inherited through a deep-bound kernel") {
      result<<"OK\n";
    } else {
      result<<"compiler child lost the host environment\n";
    }
  } catch(const std::exception& e) {
    result<<e.what()<<'\n';
  }
}
