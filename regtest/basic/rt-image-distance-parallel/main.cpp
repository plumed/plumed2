#include "test-tools.h"
#include "mpi-test.h"
#include "plumed/tools/OpenMP.h"
#include <iomanip>

namespace {
void same(MPI_Comm world,double value,const X& force,const B& virial,
          double expected,const X& expectedForce,const B& expectedVirial) {
  const auto identical=[](double a,double b) {
    return std::memcmp(&a,&b,sizeof(double))==0;
  };
  int ok=identical(value,expected) && force.size()==expectedForce.size();
  for(unsigned i=0; i<force.size(); ++i)
    for(unsigned k=0; k<3; ++k) {
      ok=ok && identical(force[i][k],expectedForce[i][k]);
    }
  for(unsigned k=0; k<9; ++k) {
    ok=ok && identical(virial[k],expectedVirial[k]);
  }
  int all=0;
  MPI_Allreduce(&ok,&all,1,MPI_INT,MPI_MIN,world);
  require(all,"serial/parallel value, forces and virial must agree bitwise");
}

int run(MPI_Comm world,int argc,char** argv) {
  int rank=0,size=0;
  MPI_Comm_rank(world,&rank);
  MPI_Comm_size(world,&size);
  const auto setup=[&](PLMD::Plumed& p) {
    p.cmd("setMPIComm",&world);
  };
  const std::string logfile="rank-"+std::to_string(rank)+".log";
  const B rect{{4,0,0,0,5,0,0,0,6}};
  const B cubic{{4,0,0,0,4,0,0,0,4}};
  const B skew{{4,0,0,1.3,3.7,0,-0.8,0.7,5.1}};
  const B adversarial{{4,0,0,12,0.7,0,0,0,6}};
  std::mt19937 rng(356);
  std::uniform_real_distribution<double> uniform(-0.7,0.7);
  std::ofstream out;
  if(rank==0) {
    out.open("output");
  }
  if(argc>1 && std::string(argv[1])=="--benchmark") {
    const std::string mode=argc>2 ? argv[2] : "parallel";
    const bool reference=std::getenv("IMAGE_DISTANCE_REFERENCE_BUILTIN");
    require(mode!="reference" || reference,"use benchmark.py to build the reference benchmark");
    for(unsigned n: {
          100u,1000u,5000u
        }) {
      // Reset for every size and mode: all processes and runs use identical input.
      rng.seed(356+n);
      X x(n);
      for(auto& point:x)
        for(double& v:point) {
          v=uniform(rng);
        }
      const unsigned frames=n<200 ? 1000 : n<2000 ? 100 : 8;
      Harness h(n,"",false,false,mode=="serial" ? "SERIAL" : "",setup,logfile,
                mode=="reference" ? "IMAGE_DISTANCE_REFERENCE" : "IMAGE_DISTANCE");
      double distance=0;
      for(unsigned warmup=0; warmup<3; ++warmup) {
        distance=h.sample(x,rect);
      }
      std::array<double,3> timings{};
      for(double& ms:timings) {
        MPI_Barrier(world);
        const auto start=std::chrono::steady_clock::now();
        for(unsigned frame=0; frame<frames; ++frame) {
          h.sample(x,rect);
        }
        const double elapsed=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count()/frames;
        MPI_Allreduce(&elapsed,&ms,1,MPI_DOUBLE,MPI_MAX,world);
      }
      std::sort(timings.begin(),timings.end());
      if(rank==0) {
        std::cout<<mode<<","<<size<<","<<PLMD::OpenMP::getNumThreads()<<","<<n<<","<<frames<<","<<timings[1]<<","<<std::setprecision(17)<<distance<<"\n";
      }
    }
    return 0;
  }
  auto compare=[&](const X& x,const B& box,const std::string& message) {
    Harness serial(x.size(),"",false,false,"SERIAL",setup,logfile);
    Harness parallel(x.size(),"",false,false,"",setup,logfile);
    for(unsigned repeat=0; repeat<3; ++repeat) {
      const double expected=serial.sample(x,box);
      const X f=serial.forces;
      const B v=serial.virial;
      const double actual=parallel.sample(x,box);
      same(world,actual,parallel.forces,parallel.virial,expected,f,v);
    }
    if(rank==0) {
      out<<message<<" PASS\n";
    }
  };
  compare(X{{{0.2,0.3,0.4}}},rect,"one atom and empty worker ranks");
  compare(X{{{0,0,0}},{{12.4,0.3,0.2}}},rect,"extended group and empty worker ranks");
  compare(X{{{0.2,0.3,0.4}}},adversarial,"skew self image");
  X compact(513);
  for(auto& point:compact)
    for(double& v:point) {
      v=uniform(rng);
    }
  compare(compact,skew,"large triclinic work sharing");
  X tied(513,V{{0,0,0}});
  tied[33]= {{1,1,0}};
  tied[199]= {{1,1,0}};
  compare(tied,cubic,"ties across rank and thread partitions");
  X translated=compact;
  for(auto& point:translated)
    for(int k=0; k<3; ++k) {
      point[k]+=2*skew[k]-3*skew[3+k]+0.27;
    }
  compare(translated,skew,"translated triclinic work sharing");
  // Only rank zero owns the sole distinct pair. All ranks must receive the error
  // before throwing; otherwise a worker could hang in the winner collective.
  int rejected=0,total=0;
  try {
    Harness h(2,"",false,false,"",setup,logfile);
    h.sample(X{{{0,0,0}},{{1e12,0,0}}},rect);
  } catch(const std::exception&) {
    rejected=1;
  }
  MPI_Allreduce(&rejected,&total,1,MPI_INT,MPI_SUM,world);
  require(total==size,"all ranks must reject unrepresentable search bounds");
  if(rank==0) {
    out<<"collective invalid-input propagation PASS\n";
  }
  rejected=0;
  try {
    compact[129][0]=1e12;
    Harness h(compact.size(),"",false,false,"",setup,logfile);
    h.sample(compact,rect);
  } catch(const std::exception&) {
    rejected=1;
  }
  MPI_Allreduce(&rejected,&total,1,MPI_INT,MPI_SUM,world);
  require(total==size,"all ranks must reject OpenMP-worker bounds errors");
  if(rank==0) {
    out<<"threaded invalid-input propagation PASS\n";
  }
  return 0;
}
}

int main(int argc,char** argv) {
  // Initialize the shared OpenMP setting before internal thread-MPI launches its
  // rank threads. Thread count remains fixed throughout a run; test 1/2/4 threads
  // in separate invocations using PLUMED_NUM_THREADS.
  PLMD::OpenMP::getNumThreads();
  return PLMD::mpi_test::run(argc,argv,4,run);
}
