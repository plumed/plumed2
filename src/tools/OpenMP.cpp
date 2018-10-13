/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2018 The plumed team
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

static unsigned plmd_numThreads = 1;
static bool plmd_numThreads_set = false;
static unsigned plmd_cachelineSize = 512;
static bool plmd_cachelineSize_set = false;

void OpenMP::setNumThreads(const unsigned nt) {
  plmd_numThreads = nt;
  plmd_numThreads_set = true;
}

unsigned OpenMP::getCachelineSize() {
  if(!plmd_cachelineSize_set) {
    if(std::getenv("PLUMED_CACHELINE_SIZE")) Tools::convert(std::getenv("PLUMED_CACHELINE_SIZE"),plmd_cachelineSize);
    plmd_cachelineSize_set = true;
  }
  return plmd_cachelineSize;
}

unsigned OpenMP::getNumThreads() {
  if(!plmd_numThreads_set) {
    if(std::getenv("PLUMED_NUM_THREADS")) Tools::convert(std::getenv("PLUMED_NUM_THREADS"),plmd_numThreads);
    plmd_numThreads_set = true;
  }
  return plmd_numThreads;
}

unsigned OpenMP::getThreadNum() {
#if defined(_OPENMP)
  return omp_get_thread_num();
#else
  return 0;
#endif
}



}



