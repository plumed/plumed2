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
#include <cstdio>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <atomic>
#include <csignal>
#include <cstdio>

namespace PLMD {
namespace cltools {

//+PLUMEDOC TOOLS driver
/*
benchmark is a lightweight reimplementation of driver focused on running benchmarks

The main difference wrt driver is that it generates a trajectory in memory rather than reading it
from a file. This allows to better time the overhead of the plumed library, without including
the time needed to read the trajectory.

As of now it does not test cases where atoms are scattered over processors (TODO).

It is also possible to load a separate version of the plumed kernel. This enables running
benchmarks agaist previous plumed versions

\par Examples

\verbatim
plumed benchmark --plumed plumed.dat
\endverbatim

\verbatim
plumed benchmark --plumed plumed.dat --kernel /path/to/libplumedKernel.so
\endverbatim

*/
//+ENDPLUMEDOC

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

}

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
  keys.add("compulsory","--natoms","100000","the number of atoms to use for the simulation");
  keys.add("compulsory","--nsteps","2000","number of steps of MD to perform (-1 means forever)");
  keys.add("optional","--kernel","path to kernel (default=current kernel)");
}

Benchmark::Benchmark(const CLToolOptions& co ):
  CLTool(co)
{
  inputdata=commandline;
}

int Benchmark::main(FILE* in, FILE*out,Communicator& pc) {
  PlumedHandle p([&]() {
    std::string kernel;
    parse("--kernel",kernel);
    if(kernel.length()>0) return PlumedHandle::dlopen(kernel.c_str());
    else return PlumedHandle();
  }());

  if(Communicator::plumedHasMPI()) p.cmd("setMPIComm",&pc.Get_comm());
  p.cmd("setRealPrecision",(int)sizeof(double));
  p.cmd("setMDLengthUnits",1.0);
  p.cmd("setMDChargeUnits",1.0);
  p.cmd("setMDMassUnits",1.0);
  p.cmd("setMDEngine","benchmarks");
  p.cmd("setTimestep",1.0);
  std::string plumedFile; parse("--plumed",plumedFile);
  p.cmd("setPlumedDat",plumedFile.c_str());
  p.cmd("setLog",out);

  int nf; parse("--nsteps",nf);
  unsigned natoms; parse("--natoms",natoms);
  p.cmd("setNatoms",natoms); p.cmd("init");
  std::vector<double> cell( 9 ), virial( 9 );
  std::vector<Vector> pos( natoms ), forces( natoms );
  std::vector<double> masses( natoms, 1 ), charges( natoms, 0 );


  SignalHandlerGuard sigIntGuard(SIGINT, signalHandler);

  int plumedStopCondition=0;
  for(int step=0; nf<0 || step<nf; ++step) {
    for(unsigned j=0; j<natoms; ++j) pos[j] = Vector(step*j, step*j+1, step*j+2);
    p.cmd("setStep",step);
    p.cmd("setStopFlag",&plumedStopCondition);
    p.cmd("setForces",&forces[0][0],3*natoms);
    p.cmd("setBox",&cell[0],9);
    p.cmd("setVirial",&virial[0],9);
    p.cmd("setPositions",&pos[0][0],3*natoms);
    p.cmd("setMasses",&masses[0],natoms);
    p.cmd("setCharges",&charges[0],natoms);
    p.cmd("calc");
    if(plumedStopCondition || signalReceived.load()) break;
  }
  return 0;
}

} // End of namespace
}
