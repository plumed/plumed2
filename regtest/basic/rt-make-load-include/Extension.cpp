#include "core/ActionAnyorder.h"
#include "core/ActionRegister.h"
#include "value.h"
#include "extra.h"

#if MKLIB_INCLUDE_FIXTURE != 7
#error Existing PLUMED_MKLIB_CFLAGS must be preserved
#endif

namespace PLMD {
class CheckInclude : public ActionAnyorder {
public:
  static void registerKeywords(Keywords& keys) {
    ActionAnyorder::registerKeywords(keys);
    keys.add("compulsory","EXPECTED","expected header value");
  }
  explicit CheckInclude(const ActionOptions& ao): Action(ao), ActionAnyorder(ao) {
    int expected=0;
    parse("EXPECTED",expected);
    checkRead();
    if(expected!=INCLUDE_VALUE+EXTRA_VALUE) {
      error("Wrong include path or stale library");
    }
  }
};
PLUMED_REGISTER_ACTION(CheckInclude,"CHECK_INCLUDE")
}
