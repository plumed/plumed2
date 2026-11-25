/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2024 Daniele Rapetti, The plumed team

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   cudaOnPlumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   cudaOnPlumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "cudaHelpers.cuh"
namespace CUDAHELPERS {
/// finds the nearest upper multiple of the given reference
template <typename T,
          typename std::enable_if_t<std::is_integral_v<T>, bool> = true>
inline T nearestUpperMultipleTo (T number, T reference) {
  return ((number - 1) | (reference - 1)) + 1;
}

size_t idealGroups (const size_t numberOfElements,
                    const size_t runningThreads) {
  // nearest upper multiple to the numberof threads
  const size_t nnToGPU =
    nearestUpperMultipleTo (numberOfElements, runningThreads);
  /// Brent’s theorem says each thread should sum O(log n) elements
  // https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf
  // const size_t elementsPerThread=log2(runningThreads);
  const size_t expectedTotalThreads = ceil (nnToGPU / log2 (runningThreads));
  // hence the blocks should have this size:
  const unsigned ngroups =
    nearestUpperMultipleTo (expectedTotalThreads, runningThreads) /
    runningThreads;
  return ngroups;
}

size_t threadsPerBlock (const unsigned N, const unsigned maxNumThreads) {
  // this seeks the minimum number of threads to use a single block (and ends
  // the recursion)
  size_t dim = 32;
  for (dim = 32; dim < 1024; dim <<= 1) {
    if (maxNumThreads < dim) {
      dim >>= 1;
      break;
    }
    if (N < dim) {
      break;
    }
  }
  return dim;
}
} // namespace CUDAHELPERS
