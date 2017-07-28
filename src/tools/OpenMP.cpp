/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2017 The plumed team
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

unsigned OpenMP::getCachelineSize() {
  static unsigned cachelineSize=512;
  if(std::getenv("PLUMED_CACHELINE_SIZE")) Tools::convert(std::getenv("PLUMED_CACHELINE_SIZE"),cachelineSize);
  return cachelineSize;
}

unsigned OpenMP::getNumThreads() {
  static unsigned numThreads=1;
  if(std::getenv("PLUMED_NUM_THREADS")) Tools::convert(std::getenv("PLUMED_NUM_THREADS"),numThreads);
  return numThreads;
}

unsigned OpenMP::getThreadNum() {
#if defined(_OPENMP)
  return omp_get_thread_num();
#else
  return 0;
#endif
}



}



