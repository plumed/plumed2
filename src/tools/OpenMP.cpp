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
#if defined(_OPENMP)
#include <omp.h>
#endif

namespace PLMD {

namespace OpenMP {
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
