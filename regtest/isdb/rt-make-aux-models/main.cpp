#include "probe.inc"
int main() {
  for(const std::string action : {
        "rescale","caliber","shadow"
      }) {
    const auto control=evaluate(action,0);
    for(const int mode : {
          1,2,3
        }) {
      const auto trial=evaluate(action,mode);
      if(trial.normal!=control.normal || trial.history!=control.history || trial.biasFile!=control.biasFile) {
        throw std::runtime_error(action+": auxiliary model state changed");
      }
    }
  }
  std::ofstream("output")<<"RESCALE deposition/output, CALIBER regression, SHADOW reference: OK\n";
}
