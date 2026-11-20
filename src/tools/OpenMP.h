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
#ifndef __PLUMED_tools_OpenMP_h
#define __PLUMED_tools_OpenMP_h

#include <vector>

#include <cstddef>

namespace PLMD {

namespace OpenMP {

/// Keep a dynamically linked OpenMP runtime resident when the kernel is unloaded.
/// Currently needed on Linux, where runtime worker threads can outlive dlclose().
/// Set PLUMED_OPENMP_RETAIN_RUNTIME=no before constructing any PLUMED instance
/// to disable retention. Once retained, a runtime cannot be made unloadable again.
void retainRuntime();

/// Set number of threads that can be used by openMP.
/// While an initialized PLUMED instance exists (see lockNumThreads()) the number can only be decreased:
/// a larger value is ignored and false is returned.
bool setNumThreads(const unsigned nt);

/// Called when a PLUMED instance is initialized: until the matching unlockNumThreads() the number of threads
/// cannot grow. The actions size their per-thread storage when they are created, so more threads later
/// would access it out of bounds. The limit is process-wide, as the number of threads is shared by all
/// the PLUMED instances in the process.
void lockNumThreads();

/// Called when a PLUMED instance that called lockNumThreads() is destroyed
void unlockNumThreads();

/// Get number of threads that can be used by openMP
unsigned getNumThreads();

/// Returns a unique thread identification number within the current team
unsigned getThreadNum();

/// get cacheline size
unsigned getCachelineSize();

/*
How many threads a parallel loop uses.

All the threads (getNumThreads()) or one, never a number in between. With the GNU
OpenMP runtime a region with fewer threads than the previous one (but more than one)
destroys the threads that are not needed, and the next larger region creates new ones
(issue #415). When PLUMED runs inside an MD code these are the threads of the MD code:
the new ones lose the core the MD code had assigned to each of them and all run on the
core of the calling thread. Both functions below return getNumThreads() or 1.

A loop is shared only when each thread gets enough work to pay for the region, which
costs about a microsecond when the threads are awake. The thresholds assume awake threads:
when PLUMED is called by an MD code the threads may have gone to sleep, and waking them
costs tens of microseconds, once per step if all the regions of the step are close. Whether
that pays is a property of the whole input, and is chosen with the number of threads.

- Loops that do a few operations for each element of an array (copy, zero, sum, scale)
  use getGoodNumThreads(array): the work is measured by the bytes of the array.
  The default of PLUMED_CACHELINE_SIZE (512: 1 kB per thread) is lower than what these
  loops need (16 to 32 kB per thread in a first measurement): until the default is
  calibrated, some of them (data passing, atom copies, clearing derivatives) also
  require a minimum number of elements.
- Loops whose elements cost much more than reading them (a pair with a switching
  function, the search of the neighbors of an atom, a task of an action) use
  getGoodNumThreads(n,minPerThread) with the number of elements each thread needs.
*/

/// All the threads if each of them gets at least 2*getCachelineSize() bytes of an array of
/// n elements of type T (the first argument only gives the type), else one.
/// The factor two is because the array is not necessarily aligned to the cache lines.
template<typename T>
unsigned getGoodNumThreads(const T* /*getTheType*/,unsigned n) {
  const unsigned numThreads=getNumThreads();
  return n*sizeof(T)/(2*getCachelineSize())>=numThreads ? numThreads : 1;
}

/// All the threads if each of them gets at least minPerThread of the n elements of a loop, else one.
/// For loops whose elements cost much more than reading them, see above.
inline unsigned getGoodNumThreads(const std::size_t n,const std::size_t minPerThread) {
  const unsigned numThreads=getNumThreads();
  return n/(minPerThread>0 ? minPerThread : 1)>=numThreads ? numThreads : 1;
}

/// Get a reasonable number of threads so as to access to vector v
template<typename T>
unsigned getGoodNumThreads(const std::vector<T> & v) {
  if(v.size()==0) {
    return 1;
  } else {
    return getGoodNumThreads(&v[0],v.size());
  }
}

}//namespace OpenMP
}//namespace PLMD

#endif
