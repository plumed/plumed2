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

#include <chrono>
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
core of the calling thread. The functions below return getNumThreads() or 1.

A parallel region costs 1 to 6 microseconds (4 to 32 threads), so a loop uses the threads only
if each thread gets about 10 microseconds of serial work, or more. Measured with GROMACS
(GPU-resident, 4 threads, Leonardo): copy loops with 1 kB per thread made an input with many
small actions 6% slower, 8 kB per thread made an input copying 6829 atoms 16% faster and left
the first unchanged.

These thresholds assume that the threads are awake. Inside an MD code that calls PLUMED at
every step they are: the OpenMP runtime keeps them spinning for about a millisecond after a
region (OMP_WAIT_POLICY=active, which never lets them sleep, changed nothing in the runs above).
If PLUMED is called after a longer idle time, waking them costs 20 to 45 microseconds per
thread, and an input with small actions is faster with one thread (PLUMED_NUM_THREADS=1).

- Loops that do a few operations for each element of an array (copy, zero, sum, scale)
  use getGoodNumThreads(array): at least 2*PLUMED_CACHELINE_SIZE bytes per thread, 8 kB
  with the default.
- Loops whose elements cost much more than reading them (a pair with a switching
  function, a link of a molecule, the search of the neighbors of an atom) use
  getGoodNumThreads(n,minPerThread) with the number of elements each thread needs.
- Loops whose elements have a cost that is not known in advance (the tasks of an action,
  from a fraction of a microsecond to milliseconds) use a TimedLoop.
*/

/// All the threads for an array of n elements of type T (the first argument only gives the
/// type) with at least 2*getCachelineSize() bytes per thread, else one.
template<typename T>
unsigned getGoodNumThreads(const T* /*getTheType*/,unsigned n) {
  const unsigned numThreads=getNumThreads();
  return std::size_t(n)*sizeof(T)/(2*getCachelineSize())>=numThreads ? numThreads : 1;
}

/// All the threads if each of them gets at least minPerThread of the n elements of a loop, else one.
/// For loops whose elements cost much more than reading them, see above.
inline unsigned getGoodNumThreads(const std::size_t n,const std::size_t minPerThread) {
  const unsigned numThreads=getNumThreads();
  return n/(minPerThread>0 ? minPerThread : 1)>=numThreads ? numThreads : 1;
}

/// All the threads or one for a loop whose cost per element is not known in advance: it is
/// measured. The first calibrationRuns runs are serial and timed. Then the loop uses all the
/// threads if its elements give each thread at least minMicrosecondsPerThread of serial work,
/// estimated with the fastest serial run per element. The serial runs are always timed.
/// \code
/// const unsigned nt=timedLoop.getNumThreads(n);
/// const auto start=OpenMP::TimedLoop::now();
/// #pragma omp parallel for num_threads(nt)
/// ...
/// timedLoop.record(start,n,nt);
/// \endcode
class TimedLoop {
public:
  using Clock=std::chrono::steady_clock;
/// Serial work each thread needs (see above)
  static constexpr double minMicrosecondsPerThread=10.0;
/// Serial runs before the threads can be used
  static constexpr unsigned calibrationRuns=2;
/// Number of threads for a run of n elements
  unsigned getNumThreads(const std::size_t n) const {
    const unsigned all=OpenMP::getNumThreads();
    if( all==1 || n==0 || serialRuns<calibrationRuns ) {
      return 1;
    }
    return microsecondsPerElement*n>=minMicrosecondsPerThread*all ? all : 1;
  }
/// The start of a run
  static Clock::time_point now() {
    return Clock::now();
  }
/// Record a run of n elements with nt threads that began at start: only the serial runs are used
  void record(const Clock::time_point start,const std::size_t n,const unsigned nt) {
    if( nt!=1 || n==0 ) {
      return;
    }
    const double us=std::chrono::duration<double,std::micro>(Clock::now()-start).count()/n;
    if( serialRuns==0 || us<microsecondsPerElement ) {
      microsecondsPerElement=us;
    }
    if( serialRuns<calibrationRuns ) {
      serialRuns++;
    }
  }
private:
  unsigned serialRuns=0;
  double microsecondsPerElement=0.0;
};

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
