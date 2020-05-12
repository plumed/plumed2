/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2020 The plumed team
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


struct OpenMPVars {
  unsigned cacheline_size=512;
  bool cache_set=false;
  unsigned num_threads=1;
  bool nt_env_set=false;
};

static OpenMPVars & getOpenMPVars() {
  static OpenMPVars vars;
  return vars;
}

void OpenMP::setNumThreads(const unsigned nt) {
  getOpenMPVars().num_threads=nt;
}

unsigned OpenMP::getCachelineSize() {
  if(!getOpenMPVars().cache_set) {
    if(std::getenv("PLUMED_CACHELINE_SIZE")) Tools::convert(std::getenv("PLUMED_CACHELINE_SIZE"),getOpenMPVars().cacheline_size);
    getOpenMPVars().cache_set = true;
  }
  return getOpenMPVars().cacheline_size;
}

unsigned OpenMP::getNumThreads() {
  if(!getOpenMPVars().nt_env_set) {
    if(std::getenv("PLUMED_NUM_THREADS")) Tools::convert(std::getenv("PLUMED_NUM_THREADS"),getOpenMPVars().num_threads);
    getOpenMPVars().nt_env_set = true;
  }
  return getOpenMPVars().num_threads;
}

unsigned OpenMP::getThreadNum() {
#if defined(_OPENMP)
  return omp_get_thread_num();
#else
  return 0;
#endif
}



}



