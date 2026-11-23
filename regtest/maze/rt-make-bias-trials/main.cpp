#include "probe.inc"
int main() {
  const auto control=evaluate("maze-bias",0);
  for(const int mode : {
        1,2,3
      }) {
    const auto trial=evaluate("maze-bias",mode);
    if(trial.normal!=control.normal || trial.history!=control.history) {
      throw std::runtime_error("MAZE bias history changed on trial");
    }
  }
  if(evaluate("maze-no-print",0).normal!=control.normal) {
    throw std::runtime_error("MAZE optimizer scheduling depends on PRINT");
  }
  std::ofstream("output")<<"MAZE bias state and explicit optimizer dependency: OK\n";
}
