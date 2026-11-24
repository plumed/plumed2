#include "probe.inc"
int main() {
  // Discover statically linked modules before entering per-probe directories.
  PLMD::Plumed bootstrap;
  const auto control=evaluate("read",0);
  for(int mode : {
        1,2,3,6
      }) {
    if(evaluate("read",mode).normal!=control.normal) {
      throw std::runtime_error("READ ordinary history changed");
    }
  }
  // The third ordinary row has r=3. An own trial must use that same value.
  // Inspect the captured trial explicitly rather than only rejection invariance.
  std::filesystem::create_directory("own-energy");
  std::filesystem::current_path("own-energy");
  std::ostringstream capture;
  auto* saved=std::cout.rdbuf(capture.rdbuf());
  char name[]="probe", action[]="read", mode[]="2";
  char* args[]= {name,action,mode};
  runProbe(3,args);
  std::cout.rdbuf(saved);
  std::istringstream in(capture.str());
  std::string row;
  double ordinary=0,trial=0;
  while(std::getline(in,row)) {
    std::istringstream fields(row);
    std::string kind;
    int step;
    double bias;
    fields>>kind>>step>>bias;
    if(kind=="normal" && step==2) {
      ordinary=bias;
    }
    if(kind=="trial") {
      trial=bias;
    }
  }
  if(ordinary!=trial || ordinary==0) {
    throw std::runtime_error("READ consumed next-time value during own trial");
  }
  std::filesystem::current_path("..");
  std::ofstream("output")<<"READ preserves the current record during same-step auxiliary evaluation: OK\n";
}
