#include "probe.inc"
#include <cmath>
#include <sstream>
#include <stdexcept>

struct Result {
  std::vector<std::array<double,3>> normal;
  double trial=0.;
  double trialForce=0.;
};

Result evaluate(std::string action,int mode) {
  std::ostringstream capture;
  struct Restore {
    std::streambuf* saved;
    ~Restore() {
      std::cout.rdbuf(saved);
    }
  } restore {std::cout.rdbuf(capture.rdbuf())};
  std::string number=std::to_string(mode);
  char program[]="probe";
  char* argv[]= {program,action.data(),number.data()};
  if(runProbe(3,argv)!=0) {
    throw std::runtime_error("probe failed");
  }
  Result result;
  std::istringstream in(capture.str());
  std::string row;
  while(std::getline(in,row)) {
    std::istringstream values(row);
    std::string kind;
    int step;
    values>>kind>>step;
    if(kind=="normal") {
      double bias=0.,force=0.;
      values>>bias>>force;
      result.normal.push_back({double(step),bias,force});
    } else if(kind=="trial") {
      values>>result.trial>>result.trialForce;
    }
  }
  if(result.normal.size()!=5) {
    throw std::runtime_error("missing normal samples");
  }
  return result;
}

int main() {
  for(const std::string action : {
        "coordination","piv","pines","emmi","funcpath-full"
      }) {
    const auto control=evaluate(action,0),trial=evaluate(action,1),foreign=evaluate(action,2),duplicate=evaluate(action,3);
    for(int mode : {
          4,5
        }) {
      const auto early=evaluate(action,mode);
      if(early.normal!=control.normal || !std::isfinite(early.trial) ||
          std::abs(early.trial-foreign.normal[2][1])>1.e-12 ||
          !std::isfinite(early.trialForce) ||
          std::abs(early.trialForce-(mode==4 ? foreign.normal[2][2] : 0.))>1.e-12) {
        throw std::runtime_error(action+": early auxiliary context lost geometry, forces or normal state");
      }
    }
    if(control.normal!=trial.normal || control.normal!=duplicate.normal ||
        std::abs(trial.trial-foreign.normal[2][1])>1.e-12) {
      throw std::runtime_error(action+": trial geometry or normal-state restoration failed");
    }
  }
  if(evaluate("coordination",0).normal!=evaluate("coordination",6).normal) {
    throw std::runtime_error("duplicate preparation lost the post-trial list rebuild");
  }
  bool pivRejected=false;
  try {
    evaluate("piv-cadence",1);
  } catch(const std::exception& error) {
    pivRejected=std::string(error.what()).find("require UPDATEPIV=1")!=std::string::npos;
  }
  if(!pivRejected) {
    throw std::runtime_error("unsafe held PIV trial value was accepted");
  }
  bool pinesRejected=false;
  try {
    evaluate("pines-cadence",1);
  } catch(const std::exception& error) {
    pinesRejected=std::string(error.what()).find("require PL_REFRESH=1")!=std::string::npos;
  }
  if(!pinesRejected) {
    throw std::runtime_error("unsafe reduced PINES trial refresh cadence was accepted");
  }
  if(evaluate("funcpath",0).normal!=evaluate("funcpath",3).normal) {
    throw std::runtime_error("duplicate path preparation changed the normal bias");
  }
  bool rejected=false;
  try {
    evaluate("funcpath",1);
  } catch(const std::exception& error) {
    rejected=std::string(error.what()).find("disable NEIGH_SIZE")!=std::string::npos;
  }
  if(!rejected) {
    throw std::runtime_error("unsafe reduced-list path trial was accepted");
  }
  std::ofstream("output")<<"PIV, PINES and EMMI trial geometry and local restoration: OK\n"
                         <<"Path duplicate preparation and full-list trials: OK\n"
                         <<"Unsafe reduced-list path trials rejected: OK\n";
}
