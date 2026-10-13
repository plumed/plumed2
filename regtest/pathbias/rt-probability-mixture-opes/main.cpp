#include "plumed/wrapper/Plumed.h"
#include <mpi.h>
#include <array>
#include <algorithm>
#include <iterator>
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
          const std::string& sigma="0.25", int beads=2, bool mixture=false) {
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
    if(mixture) {
      for(const auto& line:std::vector<std::string> {
      "ell: CUSTOM ARG=b.bias FUNC=-x/2.49433863 PERIODIC=NO",
      "lm: PATH_LOGMEANEXP ARG=ell EXPECTED_REPLICAS="+std::to_string(beads),
        "va: CUSTOM ARG=lm FUNC=-2.49433863*x PERIODIC=NO",
        "correction: CUSTOM ARG=va,b.bias FUNC=x-y PERIODIC=NO",
        "apply: BIASVALUE ARG=correction"
      })
      p.cmd("readInputLine",line.c_str());
    }
    p.cmd("readInputLine",("PRINT ARG=b.bias,b.neff,b.zed,b.nker STRIDE=1 FMT=%24.17g FILE="+name+".cv").c_str());
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

std::string contents(const std::string& file) {
  std::ifstream in(file);
  return std::string(std::istreambuf_iterator<char>(in),std::istreambuf_iterator<char>());
}

int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  if(size!=4) {
    MPI_Abort(MPI_COMM_WORLD,2);
  }
  bool passed=true;
  std::ofstream output;
  if(rank==0) {
    output.open("output");
  }
  auto record=[&](const std::string& name,bool good) {
    int local=good?1:0,all=0;
    MPI_Allreduce(&local,&all,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    passed=passed && all;
    if(rank==0) {
      output << name << " " << (all?"PASS":"FAIL") << "\n";
    }
  };
  for(int beads: {
        2,4
      }) {
    int local=rank%(size/beads),bead=rank/(size/beads);
    MPI_Comm intra,inter;
    MPI_Comm_split(MPI_COMM_WORLD,bead,local,&intra);
    MPI_Comm_split(MPI_COMM_WORLD,local,bead,&inter);
    std::string label="P"+std::to_string(beads),train=label+"-train",frozen=label+"-frozen";
    {
      Harness h(intra,inter,local,train,true,false);
      for(int step=0; step<7; ++step) {
        h.sample(step,0.8+0.1*bead+0.03*step,true);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    std::string original=rank==0?contents(train+".state"):"";
    bool good=true;
    {
      Harness direct(intra,inter,local,label+"-direct",true,false,train+".state"," UPDATE_UNTIL=0");
      Harness mixed(intra,inter,local,frozen,true,false,train+".state"," UPDATE_UNTIL=0","0.25",beads,true);
      const double x=1.02+0.07*bead,kbt=2.49433863;
      double vb=direct.sample(0,x,true),fb=direct.f[3];
      std::vector<double> energies(beads);
      if(local==0) {
        MPI_Allgather(&vb,1,MPI_DOUBLE,energies.data(),1,MPI_DOUBLE,inter);
      }
      MPI_Bcast(energies.data(),beads,MPI_DOUBLE,0,intra);
      double minimum=*std::min_element(energies.begin(),energies.end()),sum=0;
      for(double e:energies) {
        sum+=std::exp(-(e-minimum)/kbt);
      }
      const double va=minimum-kbt*std::log(sum/beads);
      const double alpha=std::exp(-(vb-minimum)/kbt)/sum;
      const double total=mixed.sample(0,x,true),force=mixed.f[3];
      double virial=local==0?mixed.virial[0]:0;
      MPI_Allreduce(MPI_IN_PLACE,&virial,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      good=std::abs(total-va)<2e-11 && std::abs(force-alpha*fb)<2e-9 && std::abs(fb)>1e-3;
      int step=1;
      for(int changed=0; changed<beads; ++changed)
        for(double h: {
              1e-4,1e-5,1e-6
            }) {
          double ep=mixed.sample(step++,x+(changed==bead?h:0),true);
          double em=mixed.sample(step++,x-(changed==bead?h:0),true);
          if(changed==bead) {
            good=good && std::abs(force+(ep-em)/(2*h))<1e-8+1e-6*std::abs(force);
          }
        }
      const double ep=mixed.sample(step++,x*1.00001,true);
      const double em=mixed.sample(step++,x*0.99999,true);
      good=good && std::abs(virial-(ep-em)/0.00002)<1e-8+1e-6*std::abs(virial);
      good=good && std::abs(mixed.sample(step++,x,true)-total)<1e-12;
      good=good && std::abs(direct.sample(step++,x,true)-vb)<1e-12;
    }
    record(label+"_energy_force_virial",good);
    MPI_Barrier(MPI_COMM_WORLD);
    good=true;
    if(rank==0) {
      good=!original.empty() && contents(train+".state")==original;
      for(int b=0; b<beads; ++b) {
        std::ifstream in(frozen+"."+std::to_string(b)+".cv");
        std::string line;
        std::array<double,3> initial{};
        unsigned count=0;
        while(std::getline(in,line)) {
          if(line.empty() || line[0]=='#') {
            continue;
          }
          double time,bias;
          std::array<double,3> diagnostic{};
          std::istringstream row(line);
          good=good && bool(row>>time>>bias>>diagnostic[0]>>diagnostic[1]>>diagnostic[2]);
          if(count==0) {
            initial=diagnostic;
          } else {
            good=good && diagnostic==initial;
          }
          ++count;
        }
        good=good && count>10 && initial[2]>0;
      }
    }
    record(label+"_frozen_state",good);
    MPI_Comm_free(&intra);
    MPI_Comm_free(&inter);
  }
  MPI_Finalize();
  return passed?0:1;
}
