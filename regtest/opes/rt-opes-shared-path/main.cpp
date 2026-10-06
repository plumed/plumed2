#include "plumed/wrapper/Plumed.h"
#include <mpi.h>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// A configurable number of synchronized beads, with one or more ranks per bead. The bias values
// used by the oracle are sampled before update at the deposition coordinates.
struct Harness {
  PLMD::Plumed p;
  std::array<double,6> x{{0,0,0,1,0,0}}, f{};
  std::array<double,9> box{{10,0,0,0,10,0,0,0,10}}, virial{};
  std::array<double,2> masses{{1,1}};
  Harness(MPI_Comm intra, MPI_Comm inter, int local_rank,
          const std::string& name, bool shared, bool wall,
          const std::string& restart="", const std::string& extra="",
          const std::string& sigma="0.25") {
    int natoms=2;
    double timestep=0.01, kbt=2.49433863;
    p.cmd("GREX setMPIIntracomm",&intra);
    if(local_rank==0) {
      p.cmd("GREX setMPIIntercomm",&inter);
    }
    p.cmd("GREX init");
    p.cmd("setMPIComm",&intra);
    p.cmd("setNatoms",&natoms);
    p.cmd("setTimestep",&timestep);
    p.cmd("setKbT",&kbt);
    const int restarting=restart.empty() ? 0 : 1;
    p.cmd("setRestart",&restarting);
    p.cmd("init");
    p.cmd("readInputLine","d: DISTANCE ATOMS=1,2 NOPBC");
    if(wall) {
      p.cmd("readInputLine","w: UPPER_WALLS ARG=d AT=1.1 KAPPA=3");
    }
    std::string input="b: OPES_METAD ARG=d PACE=1 TEMP=300 BARRIER=10 SIGMA="+sigma+" FIXED_SIGMA WALKERS_MPI FMT=%24.17g FILE="+name+".kernels STATE_WFILE="+name+".state STATE_WSTRIDE=1";
    if(shared) {
      input+=" WALKERS_SHARED_BIAS";
    }
    if(wall) {
      input+=" EXTRA_BIAS=w.bias";
    }
    if(!restart.empty()) {
      input+=" STATE_RFILE="+restart;
    }
    input+=extra;
    p.cmd("readInputLine",input.c_str());
  }
  double sample(int step, double distance, bool update) {
    x[3]=distance;
    f.fill(0);
    virial.fill(0);
    p.cmd("setStep",&step);
    p.cmd("setPositions",x.data());
    p.cmd("setMasses",masses.data());
    p.cmd("setBox",box.data());
    p.cmd("setForces",f.data());
    p.cmd("setVirial",virial.data());
    p.cmd("prepareCalc");
    p.cmd("performCalcNoUpdate");
    double bias=0;
    p.cmd("getBias",&bias);
    if(update) {
      p.cmd("update");
    }
    return bias;
  }
};

int main(int argc, char** argv) {
  MPI_Init(&argc,&argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  const int beads=argc>1 ? std::atoi(argv[1]) : 2;
  if(beads<2 || size%beads!=0) {
    MPI_Abort(MPI_COMM_WORLD,2);
  }
  const int bead=rank/(size/beads), local_rank=rank%(size/beads);
  const double bead_fraction=double(bead)/(beads-1);
  MPI_Comm intra, inter;
  MPI_Comm_split(MPI_COMM_WORLD,bead,local_rank,&intra);
  MPI_Comm_split(MPI_COMM_WORLD,local_rank,bead,&inter);
  std::ofstream output;
  if(rank==0) {
    output.open("output");
  }
  bool all_pass=true;
  auto report=[&](const std::string& label, bool okay) {
    int pass=okay ? 1 : 0, global=0;
    MPI_Allreduce(&pass,&global,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    all_pass=all_pass && global;
    if(rank==0) {
      output << label << " " << (global ? "PASS" : "FAIL") << "\n";
    }
  };
  for(bool wall : {
        false,true
      }) {
    for(bool shared : {
          false,true
        }) {
      const std::string name=std::string(shared ? "shared" : "local")+(wall ? "-wall" : "-bare");
      std::vector<std::vector<double>> expected;
      bool fd_pass=true;
      {
        Harness h(intra,inter,local_rank,name,shared,wall);
        for(int step=0; step<5; ++step) {
          const double energy=h.sample(step,0.8+0.35*bead_fraction+0.05*step,true);
          std::vector<double> energies(beads);
          if(local_rank==0) {
            MPI_Allgather(&energy,1,MPI_DOUBLE,energies.data(),1,MPI_DOUBLE,inter);
          }
          MPI_Bcast(energies.data(),beads,MPI_DOUBLE,0,intra);
          if(step>0) {
            if(shared) {
              double mean=0;
              for(double energy : energies) {
                mean+=energy/beads;
              }
              energies.assign(beads,mean);
            }
            for(double& e : energies) {
              e/=2.49433863; // PLUMED kBT at 300 K, kJ/mol.
            }
            expected.push_back(energies);
          }
        }
        // Freeze the field and verify F_dyn = -P*d(mean_b B_b)/dx_b.
        const double distance=1.02+0.3*bead_fraction, eps=1e-6;
        h.sample(5,distance,false);
        const double force=h.f[3];
        const double plus=h.sample(5,distance+eps,false);
        const double minus=h.sample(5,distance-eps,false);
        fd_pass=std::isfinite(force) && std::abs(force+(plus-minus)/(2*eps))<2e-7;
      }
      MPI_Barrier(MPI_COMM_WORLD);
      bool weights_pass=true;
      if(rank==0) {
        std::ifstream in(name+".kernels");
        std::string line;
        unsigned row=0;
        while(std::getline(in,line)) {
          if(line.empty() || line[0]=='#') {
            continue;
          }
          double time,center,sigma,height,logweight;
          std::istringstream parser(line);
          weights_pass=weights_pass && bool(parser>>time>>center>>sigma>>height>>logweight);
          weights_pass=weights_pass && row<expected.size()*beads;
          if(row<expected.size()*beads) {
            weights_pass=weights_pass && std::abs(logweight-expected[row/beads][row%beads])<1e-9;
          }
          ++row;
        }
        weights_pass=weights_pass && row==expected.size()*beads;
      }
      report(name+" weights",weights_pass);
      report(name+" frozen-force",fd_pass);
      bool same_restart=true;
      try {
        Harness h(intra,inter,local_rank,name+"-resume",shared,wall,name+".state");
        same_restart=std::isfinite(h.sample(5,1.1+0.2*bead_fraction,true));
      } catch(const std::exception&) {
        same_restart=false;
      }
      report(name+" same-mode-restart",same_restart);
      bool rejected=false;
      try {
        Harness h(intra,inter,local_rank,name+"-reject",!shared,wall,name+".state");
      } catch(const std::exception& e) {
        rejected=std::string(e.what()).find("weight mode mismatch")!=std::string::npos;
      }
      report(name+" changed-mode-rejected",rejected);
    }
  }
  bool mismatched_flag=false;
  try {
    Harness h(intra,inter,local_rank,"mismatched-flags",bead==0,false);
  } catch(const std::exception& e) {
    mismatched_flag=std::string(e.what()).find("must agree on every walker")!=std::string::npos;
  }
  report("mismatched-flags-rejected",mismatched_flag);
  if(rank==0) {
    std::ifstream source("shared-bare.state");
    std::ofstream changed("wrong-beads.state");
    std::string line;
    while(std::getline(source,line)) {
      if(line.find("#! SET shared_path_walkers")!=std::string::npos) {
        line="#! SET shared_path_walkers "+std::to_string(beads+1);
      }
      changed << line << "\n";
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  bool wrong_beads=false;
  try {
    Harness h(intra,inter,local_rank,"wrong-beads-restart",true,false,"wrong-beads.state");
  } catch(const std::exception& e) {
    wrong_beads=std::string(e.what()).find("bead count mismatch")!=std::string::npos;
  }
  report("changed-bead-count-rejected",wrong_beads);
  for(const std::string option : {
        "adaptive", "excluded"
      }) {
    bool rejected=false;
    try {
      const std::string extra=option=="excluded" && bead==0 ? " EXCLUDED_REGION=d" : "";
      const std::string sigma=option=="adaptive" && bead==0 ? "ADAPTIVE" : "0.25";
      Harness h(intra,inter,local_rank,"unsupported-"+option,true,false,"",extra,sigma);
    } catch(const std::exception& e) {
      rejected=std::string(e.what()).find("requires explicit SIGMA and no EXCLUDED_REGION")!=std::string::npos;
    }
    report("single-walker-"+option+"-rejected",rejected);
  }
  MPI_Comm_free(&intra);
  MPI_Comm_free(&inter);
  MPI_Finalize();
  return all_pass ? 0 : 1;
}
