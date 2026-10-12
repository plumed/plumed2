#include "plumed/wrapper/Plumed.h"
#include <mpi.h>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>

struct Harness {
  PLMD::Plumed p;
  std::array<double,9> x{{0,0,0,1,0,0,2,0,0}}, f{}, box{{10,0,0,0,10,0,0,0,10}}, virial{};
  std::array<double,3> masses{{1,1,1}};
  Harness(MPI_Comm intra, MPI_Comm inter, int local_rank, bool weighted,
          bool central, int moment, int power, bool use_moment) {
    const int natoms=3;
    const double timestep=0.01, kbt=2.49433863;
    p.cmd("GREX setMPIIntracomm",&intra);
    if(local_rank==0) {
      p.cmd("GREX setMPIIntercomm",&inter);
    }
    p.cmd("GREX init");
    p.cmd("setMPIComm",&intra);
    p.cmd("setNatoms",&natoms);
    p.cmd("setTimestep",&timestep);
    p.cmd("setKbT",&kbt);
    p.cmd("init");
    p.cmd("readInputLine","d: DISTANCE ATOMS=1,2 NOPBC");
    p.cmd("readInputLine","e: DISTANCE ATOMS=1,3 NOPBC");
    std::string input="a: ENSEMBLE ARG="+std::string(weighted ? "d,e REWEIGHT TEMP=300" : "d");
    input+=" MOMENT="+std::to_string(moment);
    if(central) {
      input+=" CENTRAL";
    }
    if(power!=1) {
      input+=" POWER="+std::to_string(power);
    }
    p.cmd("readInputLine",input.c_str());
    p.cmd("readInputLine",use_moment ? "b: BIASVALUE ARG=a.d_m" : "b: BIASVALUE ARG=a.d");
  }
  double sample(double value, double energy) {
    const int step=0;
    x[3]=value;
    x[6]=energy;
    f.fill(0);
    virial.fill(0);
    p.cmd("setStep",&step);
    p.cmd("setPositions",x.data());
    p.cmd("setMasses",masses.data());
    p.cmd("setBox",box.data());
    p.cmd("setForces",f.data());
    p.cmd("setVirial",virial.data());
    p.cmd("calc");
    double bias=0;
    p.cmd("getBias",&bias);
    return bias;
  }
};

int main(int argc, char** argv) {
  MPI_Init(&argc,&argv);
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  const int beads=argc>1 ? std::atoi(argv[1]) : 2;
  if(beads<2 || size%beads!=0) {
    MPI_Abort(MPI_COMM_WORLD,2);
  }
  const int bead=rank/(size/beads), local_rank=rank%(size/beads);
  MPI_Comm intra,inter;
  MPI_Comm_split(MPI_COMM_WORLD,bead,local_rank,&intra);
  MPI_Comm_split(MPI_COMM_WORLD,local_rank,bead,&inter);
  std::ofstream output,details;
  if(rank==0) {
    output.open("output");
    details.open("details");
  }
  bool all_pass=true;
  for(bool weighted : {
        false,true
      }) {
    for(bool central : {
          false,true
        }) {
      for(int moment : {
            2,3
          }) {
        for(int power : {
              1,2
            }) {
          for(bool use_moment : {
                false,true
              }) {
            Harness h(intra,inter,local_rank,weighted,central,moment,power,use_moment);
            const double value=0.8+0.3*bead, energy=1.0+0.7*bead, eps=1e-6;
            h.sample(value,energy);
            const double value_force=h.f[3], energy_force=h.f[6];
            bool passed=true;
            for(int changed=0; changed<beads; ++changed) {
              const double delta=bead==changed ? eps : 0;
              const double plus=h.sample(value+delta,energy);
              const double minus=h.sample(value-delta,energy);
              if(bead==changed) {
                passed=passed && std::abs(value_force+(plus-minus)/(2*eps))<2e-7;
              }
              if(weighted) {
                const double eplus=h.sample(value,energy+delta);
                const double eminus=h.sample(value,energy-delta);
                if(bead==changed) {
                  passed=passed && std::abs(energy_force+(eplus-eminus)/(2*eps))<2e-7;
                }
              }
            }
            int local=passed ? 1 : 0, global=0;
            MPI_Allreduce(&local,&global,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
            all_pass=all_pass && global;
            if(rank==0) {
              output << (weighted ? "weighted" : "uniform") << " "
                     << (central ? "central" : "raw") << " moment=" << moment
                     << " power=" << power << " " << (use_moment ? "moment" : "mean")
                     << " " << (global ? "PASS" : "FAIL") << "\n";
              details << "value_force=" << value_force << " energy_force=" << energy_force << "\n";
            }
          }
        }
      }
    }
  }
  MPI_Comm_free(&intra);
  MPI_Comm_free(&inter);
  MPI_Finalize();
  return all_pass ? 0 : 1;
}
