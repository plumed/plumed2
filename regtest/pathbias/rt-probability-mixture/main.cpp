#include "plumed/wrapper/Plumed.h"
#include <mpi.h>
#include <algorithm>
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
  Harness(MPI_Comm intra,MPI_Comm inter,int local,int beads,double eta,
          const std::string& field="0.4*x*x+0.2*x*x*x*x",int expected=0) {
    int natoms=2;
    double dt=0.01;
    p.cmd("GREX setMPIIntracomm",&intra);
    if(local==0) {
      p.cmd("GREX setMPIIntercomm",&inter);
    }
    p.cmd("GREX init");
    p.cmd("setMPIComm",&intra);
    p.cmd("setNatoms",&natoms);
    p.cmd("setTimestep",&dt);
    p.cmd("init");
    std::vector<std::string> lines{
      "d: DISTANCE ATOMS=1,2 COMPONENTS NOPBC",
      "v: CUSTOM ARG=d.x FUNC="+field+" PERIODIC=NO",
      "ell: CUSTOM ARG=v FUNC=-x/2.3 PERIODIC=NO",
      "lm: PATH_LOGMEANEXP ARG=ell EXPECTED_REPLICAS="+std::to_string(expected?expected:beads),
      "va: CUSTOM ARG=lm FUNC=-2.3*x PERIODIC=NO"
    };
    if(eta>=0) {
      lines.push_back("mean: ENSEMBLE ARG=d.x");
      lines.push_back("vc: CUSTOM ARG=mean.d.x FUNC=0.7*x*x PERIODIC=NO");
      lines.push_back("total: PROBABILITY_MIX ARG=vc,va KBT=2.3 COUPLING="+std::to_string(eta)+" LOG_NORMALIZER=0.2");
      lines.push_back("bias: BIASVALUE ARG=total");
    } else {
      lines.push_back("bias: BIASVALUE ARG=va");
    }
    for(const auto& line:lines) {
      p.cmd("readInputLine",line.c_str());
    }
  }
  double sample(double value,int step=0) {
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

double oracle(const std::vector<double>& x,double eta) {
  long double sum=0,mean=0;
  for(double v:x) {
    sum+=std::exp(-(0.4L*v*v+0.2L*v*v*v*v)/2.3L)/x.size();
    mean+=v/static_cast<long double>(x.size());
  }
  if(eta<0) {
    return -2.3L*std::log(sum);
  }
  return -2.3L*std::log((1-eta)*std::exp(-0.7L*mean*mean/2.3L)+eta*sum/std::exp(0.2L));
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
  auto record=[&](const std::string& name,bool good) {
    int local=good?1:0,all=0;
    MPI_Allreduce(&local,&all,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    passed=passed && all;
    if(rank==0) {
      output << name << " " << (all?"PASS":"FAIL") << "\n";
    }
  };
  for(int beads: {
        1,2,4
      }) {
    int local=rank%(size/beads),bead=rank/(size/beads);
    MPI_Comm intra,inter;
    MPI_Comm_split(MPI_COMM_WORLD,bead,local,&intra);
    MPI_Comm_split(MPI_COMM_WORLD,local,bead,&inter);
    std::vector<double> values;
    for(int b=0; b<beads; ++b) {
      values.push_back(-0.4+0.45*b);
    }
    for(double eta: {
          -1.0,0.0,0.5,0.9
          }) {
      Harness h(intra,inter,local,beads,eta);
      double energy=h.sample(values[bead]),force=h.f[3],virial=h.virial[0];
      bool good=std::abs(energy-oracle(values,eta))<2e-12 && std::abs(h.f[0]+force)<2e-12;
      for(int changed=0; changed<beads; ++changed) {
        for(double step: {
              1e-4,1e-5,1e-6
            }) {
          auto plus=values,minus=values;
          plus[changed]+=step;
          minus[changed]-=step;
          double ep=h.sample(plus[bead]),em=h.sample(minus[bead]);
          if(bead==changed) {
            good=good && std::abs(force+(ep-em)/(2*step))<1e-8+1e-6*std::abs(force);
          }
        }
      }
      double globalVirial=local==0?virial:0;
      MPI_Allreduce(MPI_IN_PLACE,&globalVirial,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      auto plus=values,minus=values;
      for(auto& v:plus) {
        v*=1.00001;
      }
      for(auto& v:minus) {
        v*=0.99999;
      }
      good=good && std::abs(globalVirial-(oracle(plus,eta)-oracle(minus,eta))/0.00002)<2e-8;
      Harness rebuilt(intra,inter,local,beads,eta);
      good=good && std::abs(rebuilt.sample(values[bead])-energy)<2e-12;
      auto permuted=values;
      std::rotate(permuted.begin(),permuted.begin()+1,permuted.end());
      good=good && std::abs(h.sample(permuted[bead])-energy)<2e-12;
      record("P"+std::to_string(beads)+"_eta"+std::to_string(eta),good);
    }
    {
      Harness h(intra,inter,local,beads,-1,"2300*x");
      const double value=bead==0?-1:1;
      const double energy=h.sample(value);
      const double expected=beads==1?-2300:-2300+2.3*std::log(beads);
      record("P"+std::to_string(beads)+"_extreme",std::isfinite(h.f[3]) && std::abs(energy-expected)<1e-10);
    }
    bool caught=false;
    try {
      Harness h(intra,inter,local,beads,-1,"x",beads+1);
    } catch(...) {
      caught=true;
    }
    record("P"+std::to_string(beads)+"_wrong_count",caught);
    caught=false;
    try {
      Harness h(intra,inter,local,beads,-1,"log(x)");
      h.sample(bead==0?-1:1);
    } catch(...) {
      caught=true;
    }
    record("P"+std::to_string(beads)+"_nonfinite",caught);
    if(beads>1) {
      caught=false;
      try {
        Harness h(intra,inter,local,beads,-1);
        h.sample(1,bead);
      } catch(...) {
        caught=true;
      }
      record("P"+std::to_string(beads)+"_step_mismatch",caught);
    }
    MPI_Comm_free(&intra);
    MPI_Comm_free(&inter);
  }
  MPI_Finalize();
  return passed?0:1;
}
