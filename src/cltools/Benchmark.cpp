/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2019-2023 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "CLTool.h"
#include "core/CLToolRegister.h"
#include "tools/Communicator.h"
#include "tools/Tools.h"
#include "config/Config.h"
#include "tools/PlumedHandle.h"
#include "tools/Stopwatch.h"
#include "tools/Log.h"
#include "tools/DLLoader.h"
#include "tools/Random.h"
#include "tools/TrajectoryParser.h"
#include "tools/AtomDistribution.h"

#include <cstdio>
#include <string>
#include <vector>
#include <fstream>
#include <atomic>
#include <csignal>
#include <cstdio>
#include <random>
#include <algorithm>
#include <chrono>
#include <string_view>
#include <optional>

namespace PLMD {
namespace cltools {

//+PLUMEDOC TOOLS benchmark
/*
benchmark is a lightweight reimplementation of [driver](driver.md) that can be used to run benchmark calculations

The main difference between [driver](driver.md) and benchmark is that benchmark generates a trajectory in memory rather than reading a
trajectory from a file. This approach is better for timing the overhead of the plumed library.  If you do similar benchmarking with driver
the timings you get are dominated by the time spent doing the I/O operations that are required to read the trajectory.

## Basic usage

If you want to use benchmark you first create a sample `plumed.dat` file for testing. For example:

```plumed
WHOLEMOLECULES ENTITY0=1-10000
p: POSITION ATOM=1
RESTRAINT ARG=p.x KAPPA=1 AT=0
```

You can then run this benchmark using the following command:

```plumed
plumed benchmark
```

Notice, that benchmark will read an input file called `plumed.dat` by default.  You can specify a different name for you PLUMED input file
by using the `--plumed` flag.

## Running with a different PLUMED version

If you want to run a benchmark against a previous plumed version in a controlled setting you can do so by using the command:

```plumed
plumed-runtime benchmark --kernel /path/to/lib/libplumedKernel.so
```

If you use this command the version of PLUMED that is in your environment calls the version of the library that is specified using the
`--kernel` flag.  Running the benchmark in this way ensures that you are running in a controlled setting, where systematic errors
in the comparison are minimized.

!!! warning "using plumed-runtime"

    You use the `plumed-runtime` executable here to avoid conflicts between different
    plumed versions. You will find the `plumed-runtime` executable in your path if you are using the non installed version of plumed,
    and in `$prefix/lib/plumed` if you installed plumed in $prefix,.

## Comparing multiple versions

The best way to compare two versions of plumed on the same input is to pass multiple colon-separated kernels as shown below:

```plumed
plumed-runtime benchmark --kernel /path/to/lib/libplumedKernel.so:/path2/to/lib/libplumedKernel.so:this
```

Here `this` means the kernel of the version with which you are running the benchmark. This comparison runs the three
instances simultaneously (alternating them) so that systematic differences in the load of your machine will affect them
to the same extent.

In case the different versions require modified plumed.dat files, or if you simply want to compare
two different plumed input files that compute the same thing, you can also use multiple plumed input files:

```plumed
plumed-runtime benchmark --kernel /path/to/lib/libplumedKernel.so:this --plumed plumed1.dat:plumed2.dat
```

Similarly, you might want to run two different inputs using the same kernel by using an input like this:

```plumed
plumed-runtime benchmark --plumed plumed1.dat:plumed2.dat
```

## Profiling

If you want to attach a profiler to the process on the fly, you might find it convenient to use `--nsteps -1`.
This options ensures that the simulation runs forever unless interrupted with CTRL-C. When interrupted, the result of the timers
should be displayed anyway.
You can also set a maximum time for the calculating by using the `--maxtime` flag.

If you run a profiler when testing multiple PLUMED versions it can be difficult to determine which function is from
each version. We therefore recommended you recompile separate PLUMED instances with a separate C++ namespace (`-DPLMD=PLUMED_version_1`)
so that you will be able to distinguish them. In addition, compiling with `CXXFLAGS="-g -O3"` will make the profiling
report more complete and will likely highlight lines of code that are particularly computationally demanding.

## MPI runs

You can also run a benchmark that emulates a domain decomposition if plumed has been compiled with MPI
and you run with `mpirun` and a command like the one shown below:

```plumed
mpirun -np 4 plumed-runtime benchmark
```

If you load separate PLUMED instances as discussed above, they should all be compiled against the same MPI version.
Notice that when using MPI signals (CTRL-C) might not work.

Since some of the data transfer could happen asynchronously, you might want to use the `--sleep` option
to simulate a lag between the `prepareCalc` and `performCalc` actions. This part of the calculation will not contribute
to the output timings, but will obviously slow down your test.

## Output

In the output you will see the usual reports about timings produced by the internal
timers of the tested plumed instances.

In addition, this tool monitors the timing externally, with some slightly different criterion:

- First, the initialization (construction of the input) will be shown with a separate timer,
  as well as the timing for the first step.
- Second, the timer corresponding to the calculation will be split in three parts, reporting
  execution of the first 20% (warm-up) and the next two blocks of 40% each.
- Finally, you might notice some discrepancies because some of the actions that are usually
  not expensive are not included in the internal timers. The external timer will
  thus provide a better estimate of the total elapsed time that includes everything.

The internal timers are still useful to monitor what happens at the different stages
of the calculattion.  If you want more detailed information you can also use a
[DEBUG](DEBUG.md) action with the `DETAILED_TIMERS`, to determine how much time is spnt in each action.

When you run multiple version, a comparative analisys of the time spent within PLUMED in the various
instances will be done.  For each PLUMED instance you run, this analysis shows the ratio between the total time each PLUMED instance ran for and the total time the first
PLUMED instance ran for. In other words, the first time that the first PLUMED instance ran for is used as the basis for comparisons. Errors on these estimates of the timings
are calculated using bootstrapping and the warm-up phase is discarded in the analysis.

*/
//+ENDPLUMEDOC

// We use an anonymous namespace here to avoid clashes with variables
// declared in other parts of the code
namespace {

//this is a sugar for changing idea faster about the rng
using generator = std::mt19937;

std::atomic<bool> signalReceived(false);

class SignalHandlerGuard {
public:
  SignalHandlerGuard(int signal, void (*newHandler)(int)) : signal_(signal) {
    // Store the current handler before setting the new one
    prevHandler_ = std::signal(signal, newHandler);
    if (prevHandler_ == SIG_ERR) {
      throw std::runtime_error("Failed to set signal handler");
    }
  }

  ~SignalHandlerGuard() {
    // Restore the previous handler upon destruction
    std::signal(signal_, prevHandler_);
  }

  // Delete copy constructor and assignment operator to prevent copying
  SignalHandlerGuard(const SignalHandlerGuard&) = delete;
  SignalHandlerGuard& operator=(const SignalHandlerGuard&) = delete;

private:
  int signal_;
  void (*prevHandler_)(int);
};

extern "C" void signalHandler(int signal) {
  if (signal == SIGINT) {
    signalReceived.store(true);
    fprintf(stderr, "Signal interrupt received\n");
  }
  if (signal == SIGTERM) {
    signalReceived.store(true);
    fprintf(stderr, "Signal termination received\n");
  }
}

/// This base class contains members that are movable with default operations
struct KernelBase {
  std::string path;
  std::string plumed_dat;
  PlumedHandle handle;
  Stopwatch stopwatch;
  std::vector<long long int> timings;
  double comparative_timing=-1.0;
  double comparative_timing_error=-1.0;
  KernelBase(const std::string & path_,const std::string & plumed_dat_, Log* log_):
    path(path_),
    plumed_dat(plumed_dat_),
    handle([&]() {
    if(path_=="this") {
      return PlumedHandle();
    } else {
      return PlumedHandle::dlopen(path.c_str());
    }
  }()),
  stopwatch(*log_) {
  }
};

/// Local structure handling a kernel and the related timers.
/// This structure specifically contain the Log, which needs special treatment
/// in move semantics
struct Kernel :
  public KernelBase {
  Log* log=nullptr;
  Kernel(const std::string & path_,const std::string & the_plumed_dat, Log* log_):
    KernelBase(path_,the_plumed_dat,log_),
    log(log_) {
  }

  ~Kernel() {
    if(log) {
      (*log)<<"\n";
      (*log)        <<"Kernel:      "<<path<<"\n";
      (*log)        <<"Input:       "<<plumed_dat<<"\n";
      if(comparative_timing>0.0) {
        (*log).printf("Comparative: %.3f +- %.3f\n",comparative_timing,comparative_timing_error);
      }
    }
  }

  Kernel(Kernel && other) noexcept:
    KernelBase(std::move(other)),
    log(other.log) {
    other.log=nullptr; // ensure no log is done in the moved away object
  }

  Kernel & operator=(Kernel && other) noexcept {
    if(this != &other) {
      KernelBase::operator=(std::move(other));
      log=other.log;
      other.log=nullptr; // ensure no log is done in the moved away object
    }
    return *this;
  }
};


class Benchmark:
  public CLTool {
public:
  static void registerKeywords( Keywords& keys );
  explicit Benchmark(const CLToolOptions& co );
  int main(FILE* in, FILE*out,Communicator& pc) override;

  std::string description()const override {
    return "run a calculation with a fixed trajectory to find bottlenecks in PLUMED";
  }

  //this does the parsing
  std::optional<std::unique_ptr<AtomDistribution>> parseAtomDistribution(Log& log) {
    {
      std::string trajectoryFile="";
      int nn=0;
      std::string trajectory_fmt="";
      for (const auto & trj_type : TrajectoryParser::trajectoryOptions()) {
        std::string tmp;
        parse("--i"+trj_type, tmp);
        if (tmp.length()>0) {
          log << "Using --i"<<trj_type<<"=" << tmp << "\n";
          trajectory_fmt=trj_type;
          ++nn;
          trajectoryFile=tmp;
        }
      }
      bool use_molfile=false;
#ifdef __PLUMED_HAS_MOLFILE_PLUGINS
      {
        auto plugins_names=TrajectoryParser::getMolfilePluginsnames() ;
        for(unsigned i=0; i<plugins_names.size(); i++) {
          std::string molfile_key="--mf_"+plugins_names[i];
          std::string traj_molfile;
          parse(molfile_key,traj_molfile);
          if(traj_molfile.length()>0) {
            ++nn;
            log << "Using --mf_"<<plugins_names[i]<<"=" << traj_molfile << "\n";
            trajectoryFile=traj_molfile;
            trajectory_fmt=plugins_names[i];
            use_molfile=true;
          }
        }
      }
#endif
      if(nn>1) {
        std::fprintf(stderr,"ERROR: cannot provide more than one trajectory file\n");
        //let the "main"
        return std::nullopt;
      }
      if (nn==1) {
        return std::make_unique<fileTraj>(
                 trajectory_fmt,
                 trajectoryFile,
                 use_molfile,
                 -1
               );
      }
    }
    std::string atomicDistr;
    parse("--atom-distribution",atomicDistr);
    if(atomicDistr != "") {
      log << "Using --atom-distribution=" << atomicDistr << "\n";
      return getAtomDistribution(atomicDistr);
    }
    return std::nullopt;
  }

//parse and evenually decorate the AtomDistribution
  std::optional<std::unique_ptr<AtomDistribution>> createAtomDistribution(Log& log, unsigned nat) {
    auto toret= parseAtomDistribution(log);
    if (!toret.has_value()) {
      return std::nullopt;
    }
    ///@todo: add a the possibility to add the pcbbox via CLI
    /// this is necessary for some molfile plugins
    int repeatX=0;
    int repeatY=0;
    int repeatZ=0;
    parse("--repeatX",repeatX);
    log << "Using --repeatX=" << repeatX << "\n";
    parse("--repeatY",repeatY);
    log << "Using --repeatY=" << repeatY << "\n";
    parse("--repeatZ",repeatZ);
    log << "Using --repeatZ=" << repeatZ << "\n";
    if (repeatX<1 || repeatY<1 || repeatZ<1) {
      log << "ERROR: repetitions of the trajectory must be >=1\n";
      return std::nullopt;
    }
    if (repeatX*repeatY*repeatZ >1) {
      //In case it is needed
      (*toret)->overrideNat(nat);
      return std::make_unique<repliedTrajectory>(std::move(*toret),
             repeatX,
             repeatY,
             repeatZ,
             nat
                                                );
    }
    return toret;
  }
};

PLUMED_REGISTER_CLTOOL(Benchmark,"benchmark")

void Benchmark::registerKeywords( Keywords& keys ) {
  CLTool::registerKeywords( keys );
  keys.add("compulsory","--plumed","plumed.dat","colon separated path(s) to the input file(s)");
  keys.add("compulsory","--kernel","this","colon separated path(s) to kernel(s)");
  keys.add("compulsory","--natoms","100000","the number of atoms to use for the simulation");
  // Maybe "--natoms" can be more clear when calling --help if we use reset_style to "atoms"
  keys.add("compulsory","--nsteps","2000","number of steps of MD to perform (-1 means forever)");
  keys.add("compulsory","--maxtime","-1","maximum number of seconds (-1 means forever)");
  keys.add("compulsory","--sleep","0","number of seconds of sleep, mimicking MD calculation");
  keys.add("compulsory","--atom-distribution","line","the kind of possible atomic displacement at each step");
  // Maybe "--atom-distribution" can be more clear when calling --help if we use reset_style to "atoms"
  keys.add("optional","--dump-trajectory","dump the trajectory to this file");
  keys.addFlag("--domain-decomposition",false,"simulate domain decomposition, implies --shuffle");
  keys.addFlag("--shuffled",false,"reshuffle atoms");
  TrajectoryParser::registerKeywords(keys);
  keys.add("compulsory","--repeatX","1","number of time to align the read trajectory along the fist box component,"
           " ingnored with a atomic distribution");
  keys.add("compulsory","--repeatY","1","number of time to align the read trajectory along the second box component,"
           " ingnored with a atomic distribution");
  keys.add("compulsory","--repeatZ","1","number of time to align the read trajectory along the third box component,"
           " ingnored with a atomic distribution");
}

Benchmark::Benchmark(const CLToolOptions& co ):
  CLTool(co) {
  inputdata=inputType::commandline;
}


int Benchmark::main(FILE* in, FILE*out,Communicator& pc) {
  // deterministic initializations to avoid issues with MPI
  generator rng;
  PLMD::Random atomicGenerator;

  struct FileDeleter {
    void operator()(FILE*f) const noexcept {
      if(f) {
        std::fclose(f);
      }
    }
  };

  std::unique_ptr<FILE,FileDeleter> log_dev_null{std::fopen("/dev/null","w")};

  Log log;
  if(pc.Get_rank()==0) {
    log.link(out);
  } else {
    log.link(log_dev_null.get());
  }
  log.setLinePrefix("BENCH:  ");
  log <<"Welcome to PLUMED benchmark\n";
  std::vector<Kernel> kernels;

  // perform comparative analysis
  // ensure that kernels vector is destroyed from last to first element upon exit
  auto kernels_deleter=[&log](auto f) {
    if(!f) {
      return;
    }
    if(f->empty()) {
      return;
    }
    generator bootstrapRng;

    const auto size=f->back().timings.size();
    //B are the bootstrap iterations
    constexpr int B=200;
    const size_t numblocks=size;
    // for some reasons, blocked bootstrap underestimates error
    // For now I keep it off. If I remove it, the code below could be simplified
    // if(numblocks>20) numblocks=20;
    const auto blocksize=size/numblocks;

    if(f->size()<2) {
      log<<"Single run, skipping comparative analysis\n";
    } else if(size<10) {
      log<<"Too small sample, skipping comparative analysis\n";
    } else
      try {

        log<<"Running comparative analysis, "<<numblocks<<" blocks with size "<<blocksize<<"\n";

        std::vector<std::size_t> choice(size);
        std::uniform_int_distribution<> distrib(0, numblocks-1);
        std::vector<std::vector<long long int>> blocks(f->size());

        {
          int i=0;
          for(auto it = f->rbegin(); it != f->rend(); ++it,++i) {
            size_t l=0;
            blocks[i].assign(numblocks,0);
            for(auto j=0ULL; j<numblocks; j++) {
              for(auto k=0ULL; k<blocksize; k++) {
                plumed_assert(l<it->timings.size());
                blocks[i][j]+=it->timings[l];
                l++;
              }
            }
          }
        }

        std::vector<std::vector<double>> ratios(f->size());
        for(auto & r : ratios) {
          //B are the bootstrap iterations
          r.resize(B);
        }

        //B are the bootstrap iterations
        for(unsigned b=0; b<B; b++) {
          for(auto & c : choice) {
            c=distrib(bootstrapRng);
          }
          long long int reference=0;
          for(const auto & c : choice) {
            reference+=blocks[0][c];
          }
          for(auto i=0ULL; i<blocks.size(); i++) {
            long long int estimate=0;
            // this would lead to separate bootstrap samples for each estimate:
            // for(auto & c : choice){c=distrib(bootstrapRng);}
            for(const auto & c : choice) {
              estimate+=blocks[i][c];
            }
            ratios[i][b]=double(estimate)/double(reference);
          }
        }

        {
          int i=0;
          for(auto it = f->rbegin(); it != f->rend(); ++it,++i) {
            double sum=0.0;
            double sum2=0.0;
            for(auto r : ratios[i]) {
              sum+=r;
              sum2+=r*r;
            }
            //B are the bootstrap iterations
            it->comparative_timing=sum/B;
            it->comparative_timing_error=std::sqrt(sum2/B-sum*sum/(B*B));
          }
        }

      } catch(std::exception & e) {
        log<<"Unexpected error during comparative analysis\n";
        log<<e.what()<<"\n";
      }
    while(!f->empty()) {
      f->pop_back();
    }

  };
  std::unique_ptr<decltype(kernels),decltype(kernels_deleter)> kernels_deleter_obj(&kernels,kernels_deleter);


  // construct the kernels vector:
  {
    std::vector<std::string> allpaths;

    {
      std::string paths;
      parse("--kernel",paths);
      log <<"Using --kernel=" << paths << "\n";
      allpaths=Tools::getWords(paths,":");
    }

    std::vector<std::string> allplumed;
    {
      std::string paths;
      parse("--plumed",paths);
      log <<"Using --plumed=" << paths << "\n";
      allplumed=Tools::getWords(paths,":");
    }

    plumed_assert(allplumed.size()>0 && allpaths.size()>0);

    // this check only works on MacOS
#if defined(__APPLE__)
    // if any of the paths if different from "this", we check if libplumed was loaded locally to avoid conflicts.
    if(std::any_of(allpaths.begin(),allpaths.end(),[](auto value) {
    return value != "this";
  })) {
      if(DLLoader::isPlumedGlobal()) {
        plumed_error()<<"It looks like libplumed is loaded in the global namespace, you cannot load a different version of the kernel\n"
                      <<"Please make sure you use the plumed-runtime executable and that the env var PLUMED_LOAD_NAMESPACE is not set to GLOBAL";
      }
    }
#endif

    if(allplumed.size()>1 && allpaths.size()>1 && allplumed.size() != allpaths.size()) {
      plumed_error() << "--kernel and --plumed should have either one element or the same number of elements";
    }

    if(allplumed.size()>1 && allpaths.size()==1)
      for(unsigned i=1; i<allplumed.size(); i++) {
        allpaths.push_back(allpaths[0]);
      }
    if(allplumed.size()==1 && allpaths.size()>1)
      for(unsigned i=1; i<allpaths.size(); i++) {
        allplumed.push_back(allplumed[0]);
      }

    for(unsigned i=0; i<allpaths.size(); i++) {
      kernels.emplace_back(allpaths[i],allplumed[i],&log);
    }
  }

  // reverse order so that log happens in the forward order:
  std::reverse(kernels.begin(),kernels.end());

  // read other flags:
  bool shuffled=false;
  parseFlag("--shuffled",shuffled);

  int nf;
  parse("--nsteps",nf);
  log << "Using --nsteps=" << nf << "\n";
  unsigned natoms;
  parse("--natoms",natoms);
  log << "Using --natoms=" << natoms << "\n";
  double maxtime;
  parse("--maxtime",maxtime);
  log << "Using --maxtime=" << maxtime << "\n";

  bool domain_decomposition=false;
  parseFlag("--domain-decomposition",domain_decomposition);

  if(pc.Get_size()>1) {
    domain_decomposition=true;
  }
  if(domain_decomposition) {
    shuffled=true;
  }

  if (shuffled) {
    log << "Using --shuffled\n";
  }
  if (domain_decomposition) {
    log << "Using --domain-decomposition\n";
  }

  double timeToSleep;
  parse("--sleep",timeToSleep);
  log << "Using --sleep=" << timeToSleep << "\n";

  std::vector<int> shuffled_indexes;
  std::unique_ptr<AtomDistribution> distribution;
  if(auto checkDistr = createAtomDistribution(log,natoms);
      checkDistr.has_value()) {
    distribution = std::move (*checkDistr);
    if (distribution->overrideNat(natoms)) {
      log << "Distribution overrode --natoms, Using --natoms=" << natoms << "\n";
    }
  } else {
    std::fprintf(stderr,"ERROR: problem with setting up the trajectory for the benchmark\n");
    return 1;
  }

  {
    std::string fileToDump;
    if(parse("--dump-trajectory",fileToDump)) {
      log << "Saving the trajectory to \"" << fileToDump << "\" and exiting\n";
      std::vector<double> cell(9);
      std::vector<Vector> pos(natoms);
      std::ofstream ofile(fileToDump);
      if (nf<0) {
        //if the user accidentally sets infinite steps, we set it to print only one
        nf=1;
      }
      for(int step=0; step<nf; ++step) {
        auto sw=kernels[0].stopwatch.startStop("TrajectoryGeneration");
        distribution->frame(pos,cell,step,atomicGenerator);

        ofile << natoms << "\n"
              << cell[0] << " " << cell[1] << " " << cell[2] << " "
              << cell[3] << " " << cell[4] << " " << cell[5] << " "
              << cell[6] << " " << cell[7] << " " << cell[8] << "\n";
        for(unsigned i=0; i<natoms; ++i) {
          ofile << "X\t" << pos[i]<< "\n";
        }
      }
      ofile.close();
      return 0;
    }
  }

  log <<"Initializing the setup of the kernel(s)\n";
  const auto initial_time=std::chrono::high_resolution_clock::now();

  for(auto & k : kernels) {
    auto & p(k.handle);
    auto sw=k.stopwatch.startStop("A Initialization");
    if(Communicator::plumedHasMPI() && domain_decomposition) {
      p.cmd("setMPIComm",&pc.Get_comm());
    }
    p.cmd("setRealPrecision",(int)sizeof(double));
    p.cmd("setMDLengthUnits",1.0);
    p.cmd("setMDChargeUnits",1.0);
    p.cmd("setMDMassUnits",1.0);
    p.cmd("setMDEngine","benchmarks");
    p.cmd("setTimestep",1.0);
    p.cmd("setPlumedDat",k.plumed_dat.c_str());
    p.cmd("setLog",out);
    p.cmd("setNatoms",natoms);
    p.cmd("init");
  }

  std::vector<double> cell( 9 ), virial( 9 );
  std::vector<Vector> pos( natoms ), forces( natoms );
  std::vector<double> masses( natoms, 1 ), charges( natoms, 0 );

  if(shuffled) {
    shuffled_indexes.resize(natoms);
    for(unsigned i=0; i<natoms; i++) {
      shuffled_indexes[i]=i;
    }
    std::shuffle(shuffled_indexes.begin(),shuffled_indexes.end(),rng);
  }

  // non owning pointers, used for shuffling the execution order
  std::vector<Kernel*> kernels_ptr;
  for(unsigned i=0; i<kernels.size(); i++) {
    kernels_ptr.push_back(&kernels[i]);
  }

  int plumedStopCondition=0;
  bool fast_finish=false;
  int part=0;

  log<<"Starting MD loop\n";
  log<<"Use CTRL+C to stop at any time and collect timers (not working in MPI runs)\n";
  // trap signals:
  SignalHandlerGuard sigIntGuard(SIGINT, signalHandler);
  SignalHandlerGuard sigTermGuard(SIGTERM, signalHandler);

  for(int step=0; nf<0 || step<nf; ++step) {
    std::shuffle(kernels_ptr.begin(),kernels_ptr.end(),rng);
    distribution->frame(pos,cell,step,atomicGenerator);

    double* pos_ptr;
    double* for_ptr;
    double* charges_ptr;
    double* masses_ptr;
    int* indexes_ptr=nullptr;
    int n_local_atoms;

    if(domain_decomposition) {
      const auto nproc=pc.Get_size();
      const auto nn=natoms/nproc;
      //using int to remove warning, MPI don't work with unsigned
      int excess=natoms%nproc;
      const auto myrank=pc.Get_rank();
      auto shift=0;
      n_local_atoms=nn;
      if(myrank<excess) {
        n_local_atoms+=1;
      }
      for(int i=0; i<myrank; i++) {
        shift+=nn;
        if(i<excess) {
          shift+=1;
        }
      }
      pos_ptr=&pos[shift][0];
      for_ptr=&forces[shift][0];
      charges_ptr=&charges[shift];
      masses_ptr=&masses[shift];
      indexes_ptr=shuffled_indexes.data()+shift;
    } else {
      pos_ptr=&pos[0][0];
      for_ptr=&forces[0][0];
      charges_ptr=&charges[0];
      masses_ptr=&masses[0];
      n_local_atoms=natoms;
      indexes_ptr=shuffled_indexes.data();
    }

    const char* sw_name;
    if(part==0) {
      sw_name="B0 First step";
    } else if(part==1) {
      sw_name="B1 Warm-up";
    } else if(part==2) {
      sw_name="B2 Calculation part 1";
    } else {
      sw_name="B3 Calculation part 2";
    }


    for(unsigned i=0; i<kernels_ptr.size(); i++) {
      auto & p(kernels_ptr[i]->handle);

      {
        auto sw=kernels_ptr[i]->stopwatch.startPause(sw_name);
        p.cmd("setStep",step);
        p.cmd("setStopFlag",&plumedStopCondition);
        p.cmd("setForces",for_ptr, {n_local_atoms,3});
        p.cmd("setBox",&cell[0], {3,3});
        p.cmd("setVirial",&virial[0], {3,3});
        p.cmd("setPositions",pos_ptr, {n_local_atoms,3});
        p.cmd("setMasses",masses_ptr, {n_local_atoms});
        p.cmd("setCharges",charges_ptr, {n_local_atoms});
        if(shuffled) {
          p.cmd("setAtomsNlocal",n_local_atoms);
          p.cmd("setAtomsGatindex",indexes_ptr, {n_local_atoms});
        }
        p.cmd("prepareCalc");
      }

      // mimick MD calculation here
      {
        unsigned k=0;
        auto start=std::chrono::high_resolution_clock::now();
        while(std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-start).count()<(long long int)1e9*timeToSleep) {
          k+=i*i;
        }
        std::fprintf(log_dev_null.get(),"%u",k);
      }

      {
        auto sw=kernels_ptr[i]->stopwatch.startStop(sw_name);
        p.cmd("performCalc");
      }

      if(kernels_ptr.size()>1 && part>1) {
        kernels_ptr[i]->timings.push_back(kernels_ptr[i]->stopwatch.getLastCycle(sw_name));
      }
      if(plumedStopCondition || signalReceived.load()) {
        fast_finish=true;
      }
    }
    auto elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-initial_time).count();
    if(part==0) {
      part=1;
    }
    if(part<2) {
      if((maxtime>0 && elapsed>(long long int)(0.2*1e9*maxtime)) || (nf>0 && step+1>=nf/5) || (maxtime<0 && nf<0 && step+1>=100)) {
        part=2;
        log<<"Warm-up completed\n";
      }
    }
    if(part<3) {
      if((maxtime>0 && elapsed>(long long int)(0.6*1e9*maxtime)) || (nf>0 && step+1>=3*nf/5)) {
        part=3;
        log<<"60% completed\n";
      }
    }

    if(maxtime>0 && elapsed>(long long int)(1e9*maxtime)) {
      fast_finish=true;
    }

    {
      unsigned tmp=fast_finish;
      pc.Bcast(tmp,0);
      fast_finish=tmp;
    }
    if(fast_finish) {
      break;
    }
  }

  return 0;
}

} // namespace unnamed
} // namespace cltools
} // namespace PLMD
