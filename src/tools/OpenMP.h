/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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

namespace PLMD{

class OpenMP{

/// Number of threads that can be used by openmp
static unsigned numThreads;
/// Empirical cacheline size (to be used to estimate correct number of threads)
static unsigned cachelineSize;

static bool initialized;

static void init();

public:

/// Get number of threads that can be used by openMP
static unsigned getNumThreads();

/// get cacheline size
static unsigned getCachelineSize();

/// Get a reasonable number of threads so as to access to an array of size s located at x
template<typename T>
static unsigned getGoodNumThreads(const T*x,unsigned s);

/// Get a reasonable number of threads so as to access to vector v;
template<typename T>
static unsigned getGoodNumThreads(const std::vector<T> & v);

};

inline
unsigned OpenMP::getCachelineSize(){
  if(!initialized) init();
  return cachelineSize;
}

inline
unsigned OpenMP::getNumThreads(){
  if(!initialized) init();
  return numThreads;
}

template<typename T>
unsigned OpenMP::getGoodNumThreads(const T*x,unsigned n){
  if(!initialized) init();
  unsigned long p=(unsigned long) x;
// a factor two is necessary since there is no guarantee that x is aligned
// to cache line boundary
  unsigned m=n/(2*cachelineSize*sizeof(T));
  if(m>numThreads) m=numThreads;
  if(m==0) m=1;
  return m;
}


template<typename T>
unsigned OpenMP::getGoodNumThreads(const std::vector<T> & v){
  if(!initialized) init();
  if(v.size()==0) return 0;
  else return getGoodNumThreads(&v[0],v.size());
}


}

#endif
