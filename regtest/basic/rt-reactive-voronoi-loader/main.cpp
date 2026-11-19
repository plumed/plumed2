#include "plumed/wrapper/Plumed.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>

namespace {

const std::string backend="./backend." __PLUMED_SOEXT;
const std::string legacy="./legacy." __PLUMED_SOEXT;
const std::string neighbors=" NLIST NL_CUTOFF=1 NL_STRIDE=10 NL_SKIN=0.4";
const std::string input="VORONOI_COORDINATION CENTERS=1 ASSIGNED=2 REFERENCE=0 KAPPA=5";

void check(bool pass, const std::string& message) {
  if(!pass) {
    throw std::runtime_error(message);
  }
}

struct System {
  PLMD::Plumed plumed;
  std::array<double,6> positions{{0,0,0,0.5,0,0}};
  std::array<double,6> forces{};
  std::array<double,2> masses{{1,1}};
  std::array<double,9> box{{4,0,0,0,4,0,0,0,4}}, virial{};
  double value=0.0;
  std::string file;

  explicit System(const std::string& label, int device=-1): file(label+".log") {
    plumed.cmd("setRealPrecision",8);
    plumed.cmd("setNatoms",2);
    plumed.cmd("setTimestep",0.001);
    plumed.cmd("setLogFile",file.c_str());
    if(device>=0) {
      plumed.cmd("setGpuDeviceId",device);
    }
    plumed.cmd("init");
  }

  void read(const std::string& line) {
    plumed.cmd("readInputLine",line.c_str());
  }

  void declare(const std::string& options) {
    read("q: "+input+options);
    read("bias: RESTRAINT ARG=q AT=0 KAPPA=2");
    long rank=-1;
    plumed.cmd("getDataRank q",&rank);
    plumed.cmd("setMemoryForData q",&value);
  }

  void evaluate(int step) {
    forces.fill(0.0);
    virial.fill(0.0);
    plumed.cmd("setStep",step);
    plumed.cmd("setPositions",positions.data());
    plumed.cmd("setForces",forces.data());
    plumed.cmd("setMasses",masses.data());
    plumed.cmd("setBox",box.data());
    plumed.cmd("setVirial",virial.data());
    plumed.cmd("calc");
    double bias=0.0;
    plumed.cmd("getBias",&bias);
    const auto zero=[](double x) {
      return std::abs(x)<1e-12;
    };
    check(std::abs(value-1.0)<1e-12 && std::abs(bias-1.0)<1e-12 &&
          std::all_of(forces.begin(),forces.end(),zero) &&
          std::all_of(virial.begin(),virial.end(),zero),file+": wrong value or derivatives");
  }
};

std::string log(const std::string& label) {
  std::ifstream in(label+".log");
  return {std::istreambuf_iterator<char>(in),std::istreambuf_iterator<char>()};
}

void positive(const std::string& label, const std::string& options,
              const std::string& selection, int device=-1, bool shear=false) {
  {
    System system(label,device);
    if(shear) {
      system.box[3]=0.5;
    }
    system.declare(options);
    system.evaluate(0);
    system.positions[3]=0.6;
    system.evaluate(1);
    system.positions[3]=0.8;
    system.evaluate(2);
  }
  check(log(label).find(selection)!=std::string::npos,label+": wrong backend or fallback");
}

void rejects(const std::string& label, const std::string& options,
             const std::string& reason, double separation=0.5, int device=-1) {
  try {
    System system(label,device);
    system.declare(options);
    system.positions[3]=separation;
    system.evaluate(0);
  } catch(const std::exception& error) {
    check(std::string(error.what()).find(reason)!=std::string::npos,
          label+": unexpected error: "+error.what());
    return;
  }
  throw std::runtime_error(label+": expected an error");
}

void exactNeighborSettings() {
  {
    System system("exact-settings");
    system.read("a: "+input+neighbors);
    system.read("b: "+input+neighbors);
    system.read("c: "+input+" NLIST NL_CUTOFF=1.0000000000000002 NL_STRIDE=10 NL_SKIN=0.4");
    system.read("d: "+input+" NLIST NL_CUTOFF=1 NL_STRIDE=10 NL_SKIN=0.4000000000000001");
  }
  const std::string report=log("exact-settings");
  const std::string sharing="sharing reactive Voronoi geometry with ";
  const auto first=report.find(sharing);
  check(first!=std::string::npos && report.find(sharing,first+sharing.size())==std::string::npos,
        "only the actions with identical neighbor settings should share");
}

}

int main() {
  try {
    exactNeighborSettings();
    const std::string cuda=" BACKEND=CUDA CUDA_LIBRARY="+backend;
    const std::string automatic=" BACKEND=AUTO CUDA_LIBRARY="+backend;
    const std::string old=" CUDA_LIBRARY="+legacy;
    positive("full",cuda,"selected CUDA (device 0",0);
    positive("buffered",automatic+neighbors,"selected CUDA (device 0; NLIST");
    positive("legacy"," BACKEND=CUDA"+old,"selected CUDA (current device");
    positive("legacy-skin"," BACKEND=AUTO"+old+neighbors,"does not support NLIST");
    positive("absent"," BACKEND=AUTO CUDA_LIBRARY=./missing-voronoi-backend","AUTO initialization fallback");
    positive("sheared",automatic,"selected CPU (triclinic", -1,true);
    rejects("legacy-device"," BACKEND=CUDA"+old,"needs the v2 factory",0.5,0);
    rejects("invalid-device",cuda,"context creation failed",0.5,1);
    rejects("coincident",automatic,"CUDA reactive Voronoi forward failed",0.0);
    rejects("outside",automatic+neighbors,"no CENTER inside the candidate list",1.2);
    std::ofstream("output")<<"PASS reactive Voronoi runtime backend loading and fallback\n";
  } catch(const std::exception& error) {
    std::cerr<<error.what()<<'\n';
    return 1;
  }
}
