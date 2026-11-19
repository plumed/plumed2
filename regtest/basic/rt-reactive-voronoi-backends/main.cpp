#include "plumed/wrapper/Plumed.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

struct Case {
  std::string name;
  std::string action;
  std::string options;
  double kappa=5.0;
  bool triclinic=false;
};

struct System {
  PLMD::Plumed plumed;
  std::vector<double> positions{
    0.2,0.3,0.4, 1.3,0.6,0.2, 2.2,1.7,0.8, 3.0,2.4,1.9,
    0.4,0.3,0.4, 0.1,0.5,0.4, 0.2,0.2,0.6,
    1.4,0.7,0.3, 2.3,1.8,0.9, 3.1,2.5,2.0
  };
  std::vector<double> forces=std::vector<double>(30,0.0);
  std::vector<double> masses=std::vector<double>(10,1.0);
  std::array<double,9> box{{7,0,0,0,8,0,0,0,9}};
  std::array<double,9> virial{};
  std::array<double,4> values{};
  double energy=0.0;

  System(const Case& specification, bool shared, const std::string& mode,
         const std::string& backend, unsigned identifier, bool legacy=false) {
    int precision=8;
    int atoms=10;
    double timestep=0.0005;
    plumed.cmd("setRealPrecision",&precision);
    plumed.cmd("setNatoms",&atoms);
    plumed.cmd("setTimestep",&timestep);
    const std::string log="system-"+std::to_string(identifier)+".log";
    plumed.cmd("setLogFile",log.c_str());
    plumed.cmd("init");
    if(specification.triclinic) {
      box[3]=0.8;
      box[6]=0.4;
      box[7]=0.7;
      // A lattice translation exercises the image used by Cartesian positions.
      positions[3]+=0.8;
      positions[4]+=8.0;
    }
    const std::string control=(shared || legacy ? "" : " NO_SHARE")+mode+backend;
    const std::string atomsInput=" CENTERS=1,2,3,4 ASSIGNED=5,6,7,8,9,10";
    const std::string reference=" REFERENCE=1.5,1.2,1.1,2.2";
    const std::string kappa=" KAPPA="+std::to_string(specification.kappa);
    auto read=[&](const std::string& line) {
      plumed.cmd("readInputLine",line.c_str());
    };
    read("a: "+specification.action+atomsInput+reference+kappa+
         " "+specification.options+control);
    read("b: VORONOI_COORDINATION"+atomsInput+reference+kappa+
         " POWER=2 SELECT=1,3 COEFFICIENTS=1.1,0.7"+control);
    read("c: VORONOI_COORDINATION"+atomsInput+reference+
         " KAPPA=8 POWER=2 SIGN=NEGATIVE"+control);
    read("d: "+specification.action+
         " CENTERS=3,1,4,2 ASSIGNED=10,8,5,9,7,6"
         " REFERENCE=1.1,1.5,2.2,1.2"+kappa+
         " "+specification.options+control);
    read("bias: RESTRAINT ARG=a,b,c,d AT=0.3,-0.1,0.2,0.3 KAPPA=3,2,1,3");
    for(unsigned index=0; index<values.size(); ++index) {
      const std::string label(1,static_cast<char>('a'+index));
      long rank=-1;
      plumed.cmd("getDataRank "+label,&rank);
      plumed.cmd("setMemoryForData "+label,&values[index]);
    }
  }

  void evaluate(int step=7) {
    std::fill(forces.begin(),forces.end(),0.0);
    virial.fill(0.0);
    plumed.cmd("setStep",&step);
    plumed.cmd("setPositions",positions.data());
    plumed.cmd("setMasses",masses.data());
    plumed.cmd("setBox",box.data());
    plumed.cmd("setForces",forces.data());
    plumed.cmd("setVirial",virial.data());
    plumed.cmd("calc");
    plumed.cmd("getBias",&energy);
  }
};

int main(int argc,char** argv) {
  const bool legacy=argc>1 && std::string(argv[1])=="--legacy";
  const std::string backend=argc>1 && !legacy ? argv[1] : "";
  if(argc>2 && std::string(argv[2])=="--timing") {
    Case specification{"tiny-cv","VORONOI_COORDINATION","POWER=2"};
    System system(specification,true," NOPBC",backend,0);
    for(int step=0; step<50; ++step) {
      system.evaluate(step);
    }
    constexpr int timedCalls=1000;
    for(int repeat=0; repeat<3; ++repeat) {
      const auto start=std::chrono::steady_clock::now();
      for(int step=0; step<timedCalls; ++step) {
        system.positions[12]+=1e-9;
        system.evaluate(50+step);
      }
      const double seconds=std::chrono::duration<double>(
                             std::chrono::steady_clock::now()-start).count();
      std::cout<<"tiny_cv_ms_per_call "<<seconds*1000.0/timedCalls<<"\n";
    }
    return 0;
  }
  const bool strictCuda=backend.find("BACKEND=CUDA")!=std::string::npos;
  std::vector<Case> cases{
    {"coord-power1","VORONOI_COORDINATION","POWER=1"},
    {"coord-power2","VORONOI_COORDINATION","POWER=2"},
    {"coord-positive","VORONOI_COORDINATION","POWER=1 SIGN=POSITIVE"},
    {"coord-negative","VORONOI_COORDINATION","POWER=2 SIGN=NEGATIVE"},
    {"coord-select","VORONOI_COORDINATION","POWER=2 SELECT=2,4 COEFFICIENTS=1.3,-0.4"},
    {"coord-diffuse","VORONOI_COORDINATION","POWER=2",0.5},
    {"coord-sharp","VORONOI_COORDINATION","POWER=2",100.0},
    {"distance-cross","VORONOI_DISTANCE","GROUP1=1,3 GROUP2=2,4"},
    {"distance-within","VORONOI_DISTANCE","GROUP1=1,2,3,4"},
    {"position-x","VORONOI_POSITION","AXIS=X ORIGIN=0.35"},
    {"position-y","VORONOI_POSITION","AXIS=Y ORIGIN=-0.7"},
    {"position-z","VORONOI_POSITION","AXIS=Z ORIGIN=1.2"},
    {"position-positive","VORONOI_POSITION","AXIS=X ORIGIN=0.35 SIGN=POSITIVE NORMALIZE"},
    {"position-negative","VORONOI_POSITION","AXIS=Y ORIGIN=-0.7 SIGN=NEGATIVE NORMALIZE"},
    {"position-positive-raw","VORONOI_POSITION","AXIS=X ORIGIN=0.35 SIGN=POSITIVE"},
    {"position-negative-raw","VORONOI_POSITION","AXIS=Y ORIGIN=-0.7 SIGN=NEGATIVE"},
    {"position-absolute-raw","VORONOI_POSITION","AXIS=Z ORIGIN=1.2 ABSOLUTE"},
    {"position-absolute","VORONOI_POSITION","AXIS=Z ORIGIN=1.2 ABSOLUTE NORMALIZE"},
    {"position-subset","VORONOI_POSITION","AXIS=X ORIGIN=-0.4 SELECT=1,4 NORMALIZE"},
    {"position-triclinic-y","VORONOI_POSITION","AXIS=Y ORIGIN=-2.7 NORMALIZE",5.0,true}
  };
  double maximumValues=0.0;
  double maximumForces=0.0;
  double maximumVirial=0.0;
  double maximumFiniteForce=0.0;
  double maximumFiniteVirial=0.0;
  unsigned evaluations=0;
  unsigned identifier=0;
  std::ofstream metrics("metrics");
  metrics<<std::setprecision(17);
  for(const auto& specification: cases) {
    for(const std::string mode: {
          "", " NOPBC",
          " NLIST NL_CUTOFF=1.85 NL_STRIDE=3 NL_SKIN=0.4",
          " NOPBC NLIST NL_CUTOFF=1.85 NL_STRIDE=1"
        }) {
      if(strictCuda && specification.triclinic && mode.find("NOPBC")==std::string::npos) {
        continue;
      }
      System candidate(specification,true,mode,backend,identifier++,legacy);
      System reference(specification,false,mode,"",identifier++,legacy);
      for(unsigned frame=0; frame<6; ++frame) {
        candidate.positions[12]+=0.003;
        if(mode.find("NL_SKIN")!=std::string::npos && frame==1) {
          // Bring a buffered-out center into the cutoff on the same step.
          candidate.positions[12]+=0.7;
          candidate.positions[13]+=0.8;
        }
        candidate.box[0]+=0.002;
        reference.positions=candidate.positions;
        reference.box=candidate.box;
        const int step=frame<4 ? 7 : 7+3*frame;
        candidate.evaluate(step);
        reference.evaluate(step);
        for(unsigned component=0; component<4; ++component) {
          maximumValues=std::max(maximumValues,
                                 std::abs(candidate.values[component]-reference.values[component]));
        }
        maximumValues=std::max(maximumValues,std::abs(candidate.energy-reference.energy));
        maximumValues=std::max(maximumValues,std::abs(candidate.values[0]-candidate.values[3]));
        for(unsigned component=0; component<30; ++component) {
          maximumForces=std::max(maximumForces,
                                 std::abs(candidate.forces[component]-reference.forces[component]));
        }
        for(unsigned component=0; component<9; ++component) {
          maximumVirial=std::max(maximumVirial,
                                 std::abs(candidate.virial[component]-reference.virial[component]));
        }
        metrics<<evaluations++<<" "<<candidate.energy;
        for(double value: candidate.values) {
          metrics<<" "<<value;
        }
        for(double force: candidate.forces) {
          metrics<<" "<<force;
        }
        for(double value: candidate.virial) {
          metrics<<" "<<value;
        }
        metrics<<"\n";
      }
      candidate.evaluate();
      const auto initialPositions=candidate.positions;
      const auto initialBox=candidate.box;
      const auto initialForces=candidate.forces;
      const auto initialVirial=candidate.virial;
      const double increment=1e-6;
      for(unsigned component=0; component<30; ++component) {
        candidate.positions=initialPositions;
        candidate.positions[component]+=increment;
        candidate.evaluate();
        const double plus=candidate.energy;
        candidate.positions[component]-=2*increment;
        candidate.evaluate();
        const double minus=candidate.energy;
        maximumFiniteForce=std::max(maximumFiniteForce,
                                    std::abs(initialForces[component]+(plus-minus)/(2*increment)));
      }
      candidate.positions=initialPositions;
      for(unsigned first=0; first<3; ++first) {
        for(unsigned second=0; second<3; ++second) {
          // A strict CUDA PBC action rejects the sheared boxes used by this FD.
          if(strictCuda &&
              mode.find("NOPBC")==std::string::npos && first!=second) {
            continue;
          }
          double energies[2];
          for(unsigned side=0; side<2; ++side) {
            const double strain=side==0 ? increment : -increment;
            candidate.positions=initialPositions;
            candidate.box=initialBox;
            for(unsigned atom=0; atom<10; ++atom) {
              candidate.positions[3*atom+second]+=strain*initialPositions[3*atom+first];
            }
            for(unsigned row=0; row<3; ++row) {
              candidate.box[3*row+second]+=strain*initialBox[3*row+first];
            }
            candidate.evaluate();
            energies[side]=candidate.energy;
          }
          // The PLUMED command interface returns the positive strain derivative.
          const double error=std::abs(initialVirial[3*first+second]-
                                      (energies[0]-energies[1])/(2*increment));
          maximumFiniteVirial=std::max(maximumFiniteVirial,error);
          if(error>2e-5) {
            std::cerr<<"virial FD "<<specification.name<<mode<<" "
                     <<first<<" "<<second<<" "<<error<<" analytic "<<initialVirial[3*first+second]<<" numeric "<<(energies[0]-energies[1])/(2*increment)<<"\n";
          }
        }
      }
    }
  }
  std::cout<<std::setprecision(12)
           <<"max_value_energy_error "<<maximumValues
           <<"\nmax_force_error "<<maximumForces
           <<"\nmax_virial_error "<<maximumVirial
           <<"\nmax_force_finite_difference_error "<<maximumFiniteForce
           <<"\nmax_virial_finite_difference_error "<<maximumFiniteVirial
           <<"\nevaluations "<<evaluations<<"\n";
  const bool pass=maximumValues<2e-9 && maximumForces<2e-8 &&
                  maximumVirial<2e-8 && maximumFiniteForce<2e-5 &&
                  maximumFiniteVirial<2e-5;
  std::ofstream output("output");
  output<<(pass ? "PASS" : "FAIL")
        <<" reactive Voronoi variants, sharing, forces and cell derivatives\n";
  return pass ? 0 : 1;
}
