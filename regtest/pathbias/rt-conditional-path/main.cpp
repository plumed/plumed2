#include "plumed/wrapper/Plumed.h"
#include <mpi.h>
#include <array>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

struct Harness {
  PLMD::Plumed p;
  std::array<double,6> x{{0,0,0,1,0,0}}, f{};
  std::array<double,9> box{{20,0,0,0,20,0,0,0,20}}, virial{};
  std::array<double,2> masses{{1,1}};
  Harness(MPI_Comm intra, MPI_Comm inter, int local, double coupling,
          const std::string& loga="log(0.1+x)",
          const std::string& logm="-1+0.2*x+0.03*x*x") {
    const int natoms=2;
    const double dt=0.01;
    p.cmd("GREX setMPIIntracomm",&intra);
    if(local==0) {
      p.cmd("GREX setMPIIntercomm",&inter);
    }
    p.cmd("GREX init");
    p.cmd("setMPIComm",&intra);
    p.cmd("setNatoms",&natoms);
    p.cmd("setTimestep",&dt);
    p.cmd("init");
    for(const std::string& line : std::vector<std::string> {
    "d: DISTANCE ATOMS=1,2 COMPONENTS NOPBC",
    "mean: ENSEMBLE ARG=d.x",
    "c: CUSTOM ARG=mean.d.x FUNC=x*x PERIODIC=NO",
    "h: CUSTOM ARG=d.x FUNC=exp(-x*x/2) PERIODIC=NO",
    "fraction: ENSEMBLE ARG=h",
    "loga: CUSTOM ARG=fraction.h FUNC="+loga+" PERIODIC=NO",
    "logm: CUSTOM ARG=c FUNC="+logm+" PERIODIC=NO",
    "ratio: CONDITIONAL_PATH ARG=c,loga,logm COUPLING="+std::to_string(coupling)+" LOWER=0 UPPER=4",
      "energy: CUSTOM ARG=c,ratio FUNC=0.7*x*x-2.3*y PERIODIC=NO",
      "bias: BIASVALUE ARG=energy"
    })
    p.cmd("readInputLine",line.c_str());
  }
  double sample(double value) {
    const int step=0;
    x[3]=value;
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

double oracle(const std::vector<double>& x,double coupling) {
  double mean=0, score=0.1;
  for(double v : x) {
    mean+=v/x.size();
    score+=std::exp(-v*v/2)/x.size();
  }
  const double c=mean*mean;
  const double normalizer=std::exp(-1+0.2*c+0.03*c*c);
  return 0.7*c*c-2.3*std::log(1-coupling+coupling*score/normalizer);
}

int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  if(size!=4) {
    MPI_Abort(MPI_COMM_WORLD,2);
  }
  std::ofstream output;
  if(rank==0) {
    output.open("output");
  }
  bool passed=true;
  auto record=[&](const std::string& name,bool local) {
    int value=local?1:0, all=0;
    MPI_Allreduce(&value,&all,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    passed=passed && all;
    if(rank==0) {
      output << name << " " << (all?"PASS":"FAIL") << "\n";
    }
  };
  for(int beads : {
        2,4
      }) {
    const int local=rank%(size/beads), bead=rank/(size/beads);
    MPI_Comm intra,inter;
    MPI_Comm_split(MPI_COMM_WORLD,bead,local,&intra);
    MPI_Comm_split(MPI_COMM_WORLD,local,bead,&inter);
    std::vector<double> values;
    for(int b=0; b<beads; ++b) {
      values.push_back(-0.4+0.45*b);
    }
    for(double coupling : {
          0.0,0.5,0.9
        }) {
      Harness h(intra,inter,local,coupling);
      const double energy=h.sample(values[bead]);
      const double force=h.f[3];
      bool good=std::abs(energy-oracle(values,coupling))<2e-12;
      good=good && std::abs(h.f[0]+force)<2e-12;
      const double step=1e-5;
      for(int changed=0; changed<beads; ++changed) {
        auto plus=values,minus=values;
        plus[changed]+=step;
        minus[changed]-=step;
        const double ep=h.sample(plus[bead]), em=h.sample(minus[bead]);
        if(bead==changed) {
          good=good && std::abs(force+(ep-em)/(2*step))<2e-8;
          good=good && std::abs(force+(oracle(plus,coupling)-oracle(minus,coupling))/(2*step))<2e-8;
        }
      }
      record("beads="+std::to_string(beads)+" coupling="+std::to_string(coupling)+" value-force",good);
      const double shifted=h.sample(values[(bead+1)%beads]);
      record("cyclic-energy",std::abs(shifted-energy)<2e-12);
      double permuted_force=0;
      MPI_Sendrecv(&force,1,MPI_DOUBLE,(bead+beads-1)%beads,42,
                   &permuted_force,1,MPI_DOUBLE,(bead+1)%beads,42,inter,MPI_STATUS_IGNORE);
      record("cyclic-force",std::abs(h.f[3]-permuted_force)<2e-12);
      Harness restarted(intra,inter,local,coupling);
      record("stateless-reconstruction",std::abs(restarted.sample(values[bead])-energy)<2e-12
             && std::abs(restarted.f[3]-force)<2e-12);
    }
    for(const std::string extreme : {
          "1000","-1000"
        }) {
      Harness h(intra,inter,local,0.5,extreme+"+0*x","0*x");
      double mean=0;
      for(double v : values) {
        mean+=v/beads;
      }
      const double expected=0.7*std::pow(mean,4)-2.3*(extreme=="1000"?1000-std::log(2):-std::log(2));
      record("extreme-"+extreme,std::abs(h.sample(values[bead])-expected)<2e-10);
    }
    {
      Harness h(intra,inter,local,0.5);
      bool rejected=false;
      try {
        h.sample(3.0);
      } catch(const std::exception&) {
        rejected=true;
      }
      record("domain-rejected",rejected);
    }
    for(double invalid : {
          -0.1,1.0
          }) {
      bool rejected=false;
      try {
        Harness h(intra,inter,local,invalid);
      } catch(const std::exception&) {
        rejected=true;
      }
      record("coupling-rejected",rejected);
    }
    MPI_Comm_free(&intra);
    MPI_Comm_free(&inter);
  }
  MPI_Finalize();
  return passed?0:1;
}
