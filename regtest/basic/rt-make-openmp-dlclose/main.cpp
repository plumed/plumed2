/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2023 The plumed team
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

#include "plumed/wrapper/Plumed.h"
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <vector>
#if defined(__linux__) && defined(__PLUMED_HAS_DLOPEN) && defined(__PLUMED_HAS_DLADDR)
#include <dlfcn.h>
#include <unistd.h>
#include <sys/wait.h>
#endif

int main(int argc, char** argv) {
#if defined(__linux__) && defined(__PLUMED_HAS_DLOPEN) && defined(__PLUMED_HAS_DLADDR) && defined(RTLD_NODELETE)
  const bool expectUnload=argc==2 && std::strcmp(argv[1], "--expect-unload")==0;
  bool kernelCanUnload=true;
  if(!expectUnload) {
    // The opt-out must be tested in a fresh process, before NODELETE is set.
    const pid_t child=fork();
    if(child<0) {
      throw std::runtime_error("Cannot fork opt-out test");
    }
    if(child==0) {
      setenv("PLUMED_OPENMP_RETAIN_RUNTIME", "no", 1);
      execl(argv[0], argv[0], "--expect-unload", static_cast<char*>(nullptr));
      _exit(127);
    }
    int status=0;
    if(waitpid(child, &status, 0)!=child || !WIFEXITED(status) ||
        (WEXITSTATUS(status)!=0 && WEXITSTATUS(status)!=77)) {
      throw std::runtime_error("Runtime-retention opt-out test failed");
    }
    // ELF loaders may retain a C++ DSO independently of our runtime policy
    // (for example because of GNU-unique symbols). The opt-out child measures
    // whether this particular kernel can unload before enforcing that property.
    kernelCanUnload=WEXITSTATUS(status)==0;
  }
  const char* path=std::getenv("PLUMED_KERNEL");
  if(!path) {
    throw std::runtime_error("PLUMED_KERNEL is required");
  }
  // Check retention before starting any worker threads: without the fix this
  // fails deterministically, rather than relying on the timing of a crash.
  plumed p=plumed_create_dlopen(path);
  if(!plumed_valid(p)) {
    throw std::runtime_error("Cannot create PLUMED");
  }
  void* handle=dlopen(path, RTLD_NOW | RTLD_LOCAL);
  void* runtime=dlsym(handle, "omp_get_thread_num");
  void* kernel=dlsym(handle, "plumed_plumedmain_create");
  dlclose(handle);
  plumed_finalize(p);
  Dl_info info{};
  const bool kernelResident=dladdr(kernel, &info)!=0;
  if(kernelResident && expectUnload) {
    // Report a baseline loader restriction, not a runtime-retention failure.
    return 77;
  }
  if(kernelResident && kernelCanUnload) {
    throw std::runtime_error("Runtime retention unexpectedly kept the PLUMED kernel loaded");
  }
  if(!kernelResident && runtime && static_cast<bool>(dladdr(runtime, &info))==expectUnload) {
    throw std::runtime_error("Unexpected OpenMP runtime retention");
  }
  if(expectUnload) {
    return 0;
  }

  for(int cycle=0; cycle<100; ++cycle) {
    p=plumed_create_dlopen(path);
    if(!plumed_valid(p)) {
      throw std::runtime_error("Cannot reload PLUMED");
    }
    int n=108, step=0;
    double dt=0.001;
    std::vector<double> positions(3*n), masses(n,1.0), forces(3*n,0.0);
    double box[9]= {6,0,0,0,6,0,0,0,6}, virial[9]= {};
    for(int i=0; i<3*n; ++i) {
      positions[i]=0.01*i;
    }
    plumed_cmd(p,"setLogFile","/dev/null");
    plumed_cmd(p,"setNatoms",&n);
    plumed_cmd(p,"setTimestep",&dt);
    plumed_cmd(p,"init",nullptr);
    plumed_cmd(p,"readInputLine","c: COORDINATION GROUPA=1-108 GROUPB=1-108 R_0=1");
    plumed_cmd(p,"readInputLine","RESTRAINT ARG=c KAPPA=0 SLOPE=1 AT=0");
    plumed_cmd(p,"setStep",&step);
    plumed_cmd(p,"setPositions",positions.data());
    plumed_cmd(p,"setMasses",masses.data());
    plumed_cmd(p,"setForces",forces.data());
    plumed_cmd(p,"setBox",box);
    plumed_cmd(p,"setVirial",virial);
    plumed_cmd(p,"calc",nullptr);
    plumed_finalize(p);
    if(kernelCanUnload && dladdr(kernel, &info)) {
      throw std::runtime_error("The PLUMED kernel did not unload after calculation");
    }
    usleep(1000);
  }
#endif
  (void)argc;
  (void)argv;
  std::ofstream("result") << "OK\n";
}
