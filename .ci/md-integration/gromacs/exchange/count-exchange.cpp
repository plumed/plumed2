// Private LOAD fixture: count forward evaluations without changing the potential.
#include "core/ActionRegister.h"
#include "function/Function.h"

namespace PLMD {
class CountExchange : public function::Function {
public:
  static void registerKeywords(Keywords& keys) {
    function::Function::registerKeywords(keys);
  }
  explicit CountExchange(const ActionOptions& ao):Action(ao),function::Function(ao) {
    addValueWithDerivatives();
    setNotPeriodic();
    checkRead();
  }
  void calculate() override {
    setValue(getArgument(0));
    setDerivative(0,1.0);
    log.printf("EXCHANGE_COUNT step=%lld exchange=%d auxiliary=%d value=%.17g\n",
               static_cast<long long>(getStep()),int(getExchangeStep()),
               int(getExchangeStep()),getArgument(0));
  }
};
PLUMED_REGISTER_ACTION(CountExchange,"COUNT_EXCHANGE")
}
