#include "plumed/wrapper/Plumed.h"
#include <fstream>
#include <map>
#include <string>
#include <vector>

// Exact restart of RANDOM_EXCHANGES: the generator of the exchange lists is stored in
// the checkpoint. The MD code asks for a list of exchanges at some steps (as GROMACS does
// at the exchange steps, getExchangesList). A reference run goes from step 0 to 39 and
// asks at the steps 5, 10, ..., 35; an interrupted run stops at step 24 after having
// written the checkpoint of step 20 (doCheckPoint at the start of that step); the restarted
// run, which starts again from step 20, must give the same lists at the steps 20 ... 35.

constexpr int nreplicas=4;

static bool ask(int step) {
  return step%5==0 && step>0;
}

// returns the list asked at each step
static std::map<int,std::vector<int>> run(int first,int last,bool restart,bool signal) {
  std::map<int,std::vector<int>> lists;
  PLMD::Plumed p;
  int natoms=2;
  double timestep=0.001;
  p.cmd("setNatoms",natoms);
  p.cmd("setMDEngine","test");
  p.cmd("setTimestep",timestep);
  p.cmd("setLogFile",restart?"log_restart":(signal?"log_part":"log_reference"));
  if(restart) {
    int r=1;
    p.cmd("setRestart",&r);
  }
  p.cmd("init");
  p.cmd("readInputLines",
        "CHECKPOINT FILE=plumed.cpt\n"
        "RANDOM_EXCHANGES SEED=7\n"
        "d: DISTANCE ATOMS=1,2\n"
        "RESTRAINT ARG=d AT=0.5 KAPPA=1.0\n");
  int nrep=nreplicas;
  p.cmd("setNumberOfReplicas",nrep);
  for(int step=first; step<last; ++step) {
    p.cmd("setStep",step);
    if(signal && step==20) {
      int one=1;
      p.cmd("doCheckPoint",one);
    }
    double positions[2][3]= {{0.0,0.0,0.0},{0.3+0.01*step,0.1,0.0}};
    double forces[2][3]= {{0.0}};
    double masses[2]= {1.0,1.0};
    double charges[2]= {0.0,0.0};
    double box[3][3]= {{2.5,0,0},{0,2.5,0},{0,0,2.5}};
    double virial[3][3]= {{0}};
    p.cmd("setPositions",&positions[0][0], {2,3});
    p.cmd("setForces",&forces[0][0], {2,3});
    p.cmd("setMasses",&masses[0], {2});
    p.cmd("setCharges",&charges[0], {2});
    p.cmd("setBox",&box[0][0], {3,3});
    p.cmd("setVirial",&virial[0][0], {3,3});
    p.cmd("calc");
    if(ask(step)) {
      std::vector<int> list(nreplicas);
      p.cmd("getExchangesList",list.data(), {static_cast<std::size_t>(nreplicas)});
      lists[step]=list;
    }
  }
  return lists;
}

int main() {
  std::ofstream out("output");
  const auto reference=run(0,40,false,true);
  run(0,25,false,true);
  const auto restarted=run(20,40,true,false);
  bool same=true;
  for(const auto & l : restarted) {
    out<<"step "<<l.first<<":";
    for(int i : l.second) {
      out<<" "<<i;
    }
    out<<"\n";
    same=same && reference.at(l.first)==l.second;
  }
  out<<"restarted lists identical to the reference: "<<(same?"YES":"NO")<<"\n";
  // without the restart the generator would start again from its seed
  bool differs=false;
  const auto fresh=run(20,40,false,false);
  for(const auto & l : fresh) {
    differs=differs || reference.at(l.first)!=l.second;
  }
  out<<"a run without the checkpoint gives different lists: "<<(differs?"YES":"NO")<<"\n";
  return 0;
}
