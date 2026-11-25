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

#include "OpenMP.h"
#include "Tools.h"
#include <cstdlib>
#include <cstring>
#if defined(_OPENMP)
#include <omp.h>
#endif
#if defined(__linux__) && defined(_OPENMP) && defined(__PLUMED_HAS_DLOPEN) && defined(__PLUMED_HAS_DLADDR)
#include <dlfcn.h>
#endif

namespace PLMD {

namespace OpenMP {

void retainRuntime() {
#if defined(__linux__) && defined(_OPENMP) && defined(__PLUMED_HAS_DLOPEN) && defined(__PLUMED_HAS_DLADDR) && defined(RTLD_NODELETE)
  const char* retain=std::getenv("PLUMED_OPENMP_RETAIN_RUNTIME");
  if(retain && std::strcmp(retain, "no")==0) {
    return;
  }
  // A kernel can be unloaded while idle OpenMP workers are still executing in
  // the runtime. Retain that runtime, rather than the kernel or its other
  // dependencies. The function-local static also serializes concurrent callers.
  static const bool retained=[] {
    Dl_info runtime{}, kernel{};
    plumed_massert(dladdr(reinterpret_cast<void*>(&omp_get_thread_num), &runtime),
                   "Cannot locate the OpenMP runtime");
    plumed_massert(dladdr(reinterpret_cast<void*>(&retainRuntime), &kernel),
                   "Cannot locate the PLUMED kernel");
    // With a statically linked OpenMP runtime, retaining this object would also
    // retain the whole kernel (or executable). Leave that case unchanged.
    if(runtime.dli_fbase==kernel.dli_fbase) {
      return true;
    }
    void* handle=dlopen(runtime.dli_fname, RTLD_NOW | RTLD_LOCAL | RTLD_NODELETE);
    plumed_massert(handle, "Cannot retain the OpenMP runtime");
    // NODELETE persists after closing this extra reference. No per-load handle
    // is leaked, and normal dlclose() of individual kernels remains possible.
    dlclose(handle);
    return true;
  }();
  (void)retained;
#endif
}

///singleton struct to treat the openMP setting as a global variables, but with a layer of encapsulation
struct OpenMPVars {
  unsigned cacheline_size=512;
  bool cache_set=false;
  unsigned num_threads=1;
  bool nt_env_set=false;
  static OpenMPVars & get() {
    static OpenMPVars vars;
    return vars;
  }
private:
  OpenMPVars()=default;
  OpenMPVars(OpenMPVars&)=delete;
  OpenMPVars& operator=(OpenMPVars&&)=delete;
};

void setNumThreads(const unsigned nt) {
  OpenMPVars::get().num_threads=nt;
}

unsigned getCachelineSize() {
  if(!OpenMPVars::get().cache_set) {
    if(std::getenv("PLUMED_CACHELINE_SIZE")) {
      Tools::convert(std::getenv("PLUMED_CACHELINE_SIZE"),OpenMPVars::get().cacheline_size);
    }
    OpenMPVars::get().cache_set = true;
  }
  return OpenMPVars::get().cacheline_size;
}

unsigned getNumThreads() {
  if(!OpenMPVars::get().nt_env_set) {
    if(std::getenv("PLUMED_NUM_THREADS")) {
      Tools::convert(std::getenv("PLUMED_NUM_THREADS"),OpenMPVars::get().num_threads);
    }
    OpenMPVars::get().nt_env_set = true;
  }
  return OpenMPVars::get().num_threads;
}

unsigned getThreadNum() {
#if defined(_OPENMP)
  return omp_get_thread_num();
#else
  return 0;
#endif
}

}//namespace OpenMP
}//namespace PLMD
