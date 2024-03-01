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
#include <cstdio>
#include <string>
#include <vector>
#include <fstream>
#include <atomic>
#include <csignal>
#include <cstdio>
#include <random>
#include <algorithm>
#if defined(__PLUMED_HAS_DLOPEN)
#include <dlfcn.h>
#endif

namespace PLMD {
namespace cltools {

//+PLUMEDOC TOOLS benchmark
/*
benchmark is a lightweight reimplementation of driver focused on running benchmarks

The main difference wrt driver is that it generates a trajectory in memory rather than reading it
from a file. This allows to better time the overhead of the plumed library, without including
the time needed to read the trajectory.

As of now it does not test cases where atoms are scattered over processors (TODO).

It is also possible to load a separate version of the plumed kernel. This enables running
benchmarks agaist previous plumed versions

\par Examples

First, you should create a sample `plumed.dat` file for testing. For instance:
```
WHOLEMOLECULES ENTITY0=1-10000
p: POSITION ATOM=1
RESTRAINT ARG=p.x KAPPA=1 AT=0

```

Then you can test the performance of this input with the following command:
\verbatim
plumed benchmark
\endverbatim

You can also test a different (older) version of PLUMED with the same input. To do so,
you should run
\verbatim
plumed-runtime benchmark --kernel /path/to/lib/libplumedKernel.so
\endverbatim

\warning It is necessary to use the `plumed-runtime` executable here to avoid conflicts between different
plumed versions. You will find it in `src/lib` if you are using the non installed version of plumed,
and in `$prefix/lib/plumed` if you installed plumed in $prefix,.

The best way to compare two versions of plumed on the same input is to pass multiple colon-separated kernels:

\verbatim
plumed-runtime benchmark --kernel /path/to/lib/libplumedKernel.so:/path2/to/lib/libplumedKernel.so:this
\endverbatim

Here `this` means the kernel of the version with which you are running the benchmark. This comparison runs the three
instances simultaneously (alternating them) so that systematic differences in the load of your machine will affect them
to the same extent.

In case the different versions require modified plumed.dat files, or if you simply want to compare
two different plumed input files that compute the same thing, you can also use multiple plumed input files:

\verbatim
plumed-runtime benchmark --kernel /path/to/lib/libplumedKernel.so:this --plumed plumed1.dat:plumed2.dat
\endverbatim

\par Profiling

If you want to attach a profiler on the fly to the process, you might find it convenient to use `--nsteps -1`.
The simulation will run forever and can be interrupted with CTRL-C. When interrupted, the result of the timers
should be displayed anyway.

\par Output

In the output you will see the usual reports about timing produced by the internal
timers of the tested plumed instances.
In addition, this tool will monitor the timing externally, with some slightly different criterion:
- First, the initialization (construction of the input) will be shown with a separate timer
- Second, the timer corresponding to the calculation will be split in two parts, reporting
  execution of the first half of the calculation and the second half.
- Finally, you might notice some discrepancy because some of the actions that are usually
  not expensive are not included in the internal timers. The external timer will
  thus provide a better estimate of the total elapsed time, including everything.

The internal timers are still useful to monitor what happens at the different stages
and, with \ref DEBUG `DETAILED_TIMERS`, what happens in each action.

*/
//+ENDPLUMEDOC

// We use an anonymous namespace here to avoid clashes with variables
// declared in other parts of the code
namespace {

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
    fprintf(stderr, "Signal handler called\n");
  }
}

/// Local structure handling a kernel and the related timers
struct Kernel {
  std::string path;
  std::string plumed_dat;
  PlumedHandle handle;
  Stopwatch stopwatch;
  Log* log=nullptr;
  Kernel(const std::string & path_,const std::string & plumed_dat, Log* log_):
    path(path_),
    plumed_dat(plumed_dat),
    stopwatch(*log_),
    log(log_)
  {
    if(path_!="this") handle=PlumedHandle::dlopen(path_.c_str());
  }
  ~Kernel() {
    if(log) {
      (*log)<<"\n";
      (*log)<<"Kernel: "<<path<<"\n";
      (*log)<<"Input:  "<<plumed_dat<<"\n";
    }
  }
  Kernel(Kernel && other) noexcept:
    path(std::move(other.path)),
    plumed_dat(std::move(other.plumed_dat)),
    handle(std::move(other.handle)),
    stopwatch(std::move(other.stopwatch)),
    log(other.log)
  {
    other.log=nullptr;
  }
  Kernel & operator=(Kernel && other) noexcept
  {
    if(this != &other) {
      path=std::move(other.path);
      plumed_dat=std::move(other.plumed_dat);
      handle=std::move(other.handle);
      stopwatch=std::move(other.stopwatch);
      log=other.log;
      other.log=nullptr;
    }
    return *this;
  }
};

class Benchmark:
  public CLTool
{
public:
  static void registerKeywords( Keywords& keys );
  explicit Benchmark(const CLToolOptions& co );
  int main(FILE* in, FILE*out,Communicator& pc) override;
  std::string description()const override {
    return "run a calculation with a fixed trajectory to find bottlenecks in PLUMED";
  }
};

PLUMED_REGISTER_CLTOOL(Benchmark,"benchmark")

void Benchmark::registerKeywords( Keywords& keys ) {
  CLTool::registerKeywords( keys );
  keys.add("compulsory","--plumed","plumed.dat","convert the input in this file to the html manual");
  keys.add("compulsory","--kernel","this","colon separated path(s) to kernel(s)");
  keys.add("compulsory","--natoms","100000","the number of atoms to use for the simulation");
  keys.add("compulsory","--nsteps","2000","number of steps of MD to perform (-1 means forever)");
  keys.addFlag("--shuffled",false,"reshuffle atoms");
}

Benchmark::Benchmark(const CLToolOptions& co ):
  CLTool(co)
{
  inputdata=commandline;
}


int Benchmark::main(FILE* in, FILE*out,Communicator& pc) {

  Log log;
  log.link(out);
  log.setLinePrefix("BENCH:  ");

  std::vector<Kernel> kernels;

  // ensure that kernels vector is destroyed from last to first element upon exit
  auto kernels_deleter=[](auto f) { while(!f->empty()) f->pop_back();};
  std::unique_ptr<decltype(kernels),decltype(kernels_deleter)> kernels_deleter_obj(&kernels,kernels_deleter);


  std::random_device rd;
  std::mt19937 g(rd());

  // construct the kernels vector:
  {
    std::vector<std::string> allpaths;

    {
      std::string paths;
      parse("--kernel",paths);
      allpaths=Tools::getWords(paths,":");
    }

    std::vector<std::string> allplumed;
    {
      std::string paths;
      parse("--plumed",paths);
      allplumed=Tools::getWords(paths,":");
    }

    plumed_assert(allplumed.size()>0 && allpaths.size()>0);

#if defined(__PLUMED_HAS_DLOPEN)
    if(std::any_of(allpaths.begin(),allpaths.end(),[](auto value) {return value != "this";})) {
#if defined(__PLUMED_HAS_RTLD_DEFAULT)
      void* handle=RTLD_DEFAULT;
#else
      handle=dlopen(NULL,RTLD_LOCAL);
#endif
      if(dlsym(handle,"plumed_plumedmain_create") || dlsym(handle,"plumedmain_create")) {
        plumed_error()<<"It looks like libplumed is loaded in the global namespace, you cannot load a different version of the kernel\n"
                      <<"Please make sure you use the plumed-runtime executable and that the env var PLUMED_LOAD_NAMESPACE is not set to GLOBAL";
      }
      if(handle) dlclose(handle);
    }
#endif

    if(allplumed.size()>1 && allpaths.size()>1 && allplumed.size() != allpaths.size()) {
      plumed_error() << "--kernel and --plumed should have either one element or the same number of elements";
    }

    if(allplumed.size()>1 && allpaths.size()==1) for(unsigned i=1; i<allplumed.size(); i++) allpaths.push_back(allpaths[0]);
    if(allplumed.size()==1 && allpaths.size()>1) for(unsigned i=1; i<allpaths.size(); i++) allplumed.push_back(allplumed[0]);

    for(unsigned i=0; i<allpaths.size(); i++) kernels.emplace_back(allpaths[i],allplumed[i],&log);
  }

  // reverse order so that log happens in the forward order:
  std::reverse(kernels.begin(),kernels.end());

  // read other flags:
  bool shuffled=false;
  parseFlag("--shuffled",shuffled);
  int nf; parse("--nsteps",nf);
  unsigned natoms; parse("--natoms",natoms);

  std::vector<int> shuffled_indexes;

  // trap signals:
  SignalHandlerGuard sigIntGuard(SIGINT, signalHandler);

  for(auto & k : kernels) {
    auto & p(k.handle);
    auto sw=k.stopwatch.startStop("A Initialization");
    if(Communicator::plumedHasMPI()) p.cmd("setMPIComm",&pc.Get_comm());
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
    for(unsigned i=0; i<natoms; i++) shuffled_indexes[i]=i;
    std::shuffle(shuffled_indexes.begin(),shuffled_indexes.end(),g);
  }

  // non owning pointers, used for shuffling the execution order
  std::vector<Kernel*> kernels_ptr;
  for(unsigned i=0; i<kernels.size(); i++) kernels_ptr.push_back(&kernels[i]);

  int plumedStopCondition=0;
  bool fast_finish=false;
  for(int step=0; nf<0 || step<nf; ++step) {
    std::shuffle(kernels_ptr.begin(),kernels_ptr.end(),g);
    for(unsigned j=0; j<natoms; ++j) pos[j] = Vector(step*j, step*j+1, step*j+2);
    for(unsigned i=0; i<kernels_ptr.size(); i++) {
      auto & p(kernels_ptr[i]->handle);

      const char* sw_name;
      if(nf<0) sw_name="B Calculation";
      else if(step<nf/2) sw_name="B1 Calculation part 1";
      else sw_name="B2 Calculation part 2";
      auto sw=kernels_ptr[i]->stopwatch.startStop(sw_name);
      p.cmd("setStep",step);
      p.cmd("setStopFlag",&plumedStopCondition);
      p.cmd("setForces",&forces[0][0],3*natoms);
      p.cmd("setBox",&cell[0],9);
      p.cmd("setVirial",&virial[0],9);
      p.cmd("setPositions",&pos[0][0],3*natoms);
      p.cmd("setMasses",&masses[0],natoms);
      p.cmd("setCharges",&charges[0],natoms);
      if(shuffled) {
        p.cmd("setAtomsNlocal",natoms);
        p.cmd("setAtomsGatindex",&shuffled_indexes[0],shuffled_indexes.size());
      }
      p.cmd("calc");
      if(plumedStopCondition || signalReceived.load()) fast_finish=true;
    }
    if(fast_finish) break;
  }

  return 0;
}

}
} // End of namespace
}
