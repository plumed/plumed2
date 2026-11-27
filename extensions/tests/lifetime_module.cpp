#include "config/Config.h"
#include <cstdio>
#include <cstdlib>

#ifdef LIFETIME_ROOT
#include "core/Action.h"
#include "core/ActionRegister.h"

namespace PLMD {
class LifetimeAction:public Action {
public:
  static void registerKeywords(Keywords& keys) {
    Action::registerKeywords(keys);
  }
  explicit LifetimeAction(const ActionOptions& options):Action(options) {
    checkRead();
  }
  ~LifetimeAction() override {
    // The automatic registry view must still be valid during member teardown.
    if(!actionRegister().check("LIFETIME_PROBE")) {
      std::abort();
    }
  }
  void calculate() override {}
  void apply() override {}
};
PLUMED_REGISTER_ACTION(LifetimeAction,"LIFETIME_PROBE")
}

// Force a native dependency, independently of ModuleLoader's opening order.
extern "C" const char* plumed_module_lifetime_dep();
extern "C" const void* lifetime_dependency() {
  return reinterpret_cast<const void*>(&plumed_module_lifetime_dep);
}
#endif

namespace {
struct Marker {
  ~Marker() {
    if(const char* path=std::getenv("PLUMED_TEST_LIFETIME_EVENTS")) {
      if(auto file=std::fopen(path,"a")) {
        std::fprintf(file,"%s\n",LIFETIME_NAME);
        std::fclose(file);
      }
    }
  }
} marker;
}
