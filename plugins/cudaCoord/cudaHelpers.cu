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
