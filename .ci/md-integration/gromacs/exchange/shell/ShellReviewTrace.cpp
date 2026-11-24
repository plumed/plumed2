#include "function/Function.h"
#include "core/ActionRegister.h"
#include <cstdio>
namespace PLMD {
namespace function {
class ShellReviewTrace : public Function {
public:
  explicit ShellReviewTrace(const ActionOptions& options): Action(options), Function(options) {
    addValueWithDerivatives();
    setNotPeriodic();
    checkRead();
  }
  static void registerKeywords(Keywords& keys) {
    Function::registerKeywords(keys);
  }
  void calculate() override {
    std::fprintf(stderr,"SHELL_REVIEW step=%lld auxiliary=%d d=%.17g min=%.17g bias=%.17g\n",static_cast<long long>(getStep()),int(isCalculationOnly()),getArgument(0),getArgument(1),getArgument(2));
    setValue(0);
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      setDerivative(i,0);
    }
  }
};
PLUMED_REGISTER_ACTION(ShellReviewTrace,"SHELLREVIEWTRACE")
}
}
