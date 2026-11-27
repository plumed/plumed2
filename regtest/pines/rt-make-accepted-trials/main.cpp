#include "probe.inc"
#include <sstream>
#include <stdexcept>
#include <vector>
std::vector<std::string> evaluate(int mode) {
  std::ostringstream capture;
  auto* saved=std::cout.rdbuf(capture.rdbuf());
  std::string number=std::to_string(mode);
  char name[]="probe", action[]="pines";
  char* args[]= {name,action,number.data()};
  runProbe(3,args);
  std::cout.rdbuf(saved);
  std::istringstream rows(capture.str());
  std::vector<std::string> result;
  std::string row;
  while(std::getline(rows,row))
    if(row.rfind("normal ",0)==0) {
      result.push_back(row);
    }
  if(result.size()!=5) {
    throw std::runtime_error("missing frames");
  }
  return result;
}
int main() {
  if(evaluate(0)!=evaluate(1)) {
    throw std::runtime_error("rejected trial changed PINES");
  }
  const auto fresh=evaluate(2), accepted=evaluate(5), jump=evaluate(7);
  for(int i=3; i<5; ++i)
    if(accepted[i]!=fresh[i] || accepted[i]!=jump[i]) {
      throw std::runtime_error("accepted swap missed nearest pair");
    }
  std::ofstream("output")<<"PINES refresh=1 accepted and rejected trials agree with current-coordinate selection: OK\n";
}
