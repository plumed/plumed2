#include "probe.inc"
int main() {
  // Discover statically linked modules before entering per-probe directories.
  PLMD::Plumed bootstrap;
  const auto control=evaluate("voronoi",0);
  for(const int mode : {
        1,2,6
      }) {
    if(evaluate("voronoi",mode).normal!=control.normal) {
      throw std::runtime_error("Voronoi trial or duplicate-prepare state leaked");
    }
  }
  std::ofstream("output")<<"Voronoi trial geometry and pending rebuild: OK\n";
}
