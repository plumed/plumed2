#include <cstdlib>
#include <fstream>
#include <vector>

#include "plumed/tools/OpenMP.h"
#include "plumed/tools/Tools.h"

using namespace PLMD;

int main(int, char **) {
  unsigned cachelineSize;
  Tools::convert(std::getenv("PLUMED_CACHELINE_SIZE"), cachelineSize);

  unsigned requestedThreads;
  Tools::convert(std::getenv("PLUMED_NUM_THREADS"), requestedThreads);

  std::ofstream out("output");
  out << "Plumed reads \"PLUMED_CACHELINE_SIZE\" correctly: "
      << (OpenMP::getCachelineSize() == cachelineSize ? "true" : "false")
      << "\n";
  out << "Plumed reads \"PLUMED_NUM_THREADS\": correctly: "
      << (OpenMP::getNumThreads() == requestedThreads ? "true" : "false")
      << "\n";
  // constexpr auto size_ofDouble= sizeof(double);
  // out << size_ofDouble <<'\n';
  out << "Checking getGoodNumThreads\n";
  out << "std::vector version:\n";
  std::vector<double> smallv(1);
  out << "getGoodNumThreads with a small vector returns 1: "
      << (OpenMP::getGoodNumThreads(smallv) == 1 ? "true" : "false") << '\n';
  // 2*cachelineSize bytes per thread
  const std::size_t nbig=2*std::size_t(cachelineSize)*requestedThreads/sizeof(double);
  std::vector<double> bigv(nbig);
  out << "getGoodNumThreads with returns \"PLUMED_NUM_THREADS\": "
      << (OpenMP::getGoodNumThreads(bigv) == requestedThreads ? "true"
                                                              : "false")
      << '\n';
  out << "Pointer version:\n";
  // getGoodNumThreads does not access to the values pointed, so passing a
  // nullptr should not give errors
  double *ptd = nullptr;

  out << "getGoodNumThreads(pointer) with a small vector returns 1: "
      << (OpenMP::getGoodNumThreads(ptd, 1) == 1 ? "true" : "false") << '\n';
  out << "getGoodNumThreads(pointer) with returns \"PLUMED_NUM_THREADS\": "
      << (OpenMP::getGoodNumThreads(ptd, nbig) == requestedThreads
              ? "true"
              : "false")
      << '\n';
  out << "getGoodNumThreads(pointer) with one element less returns 1: "
      << (OpenMP::getGoodNumThreads(ptd, nbig-1) == 1 ? "true" : "false") << '\n';
  out << "Version with a minimum number of elements per thread:\n";
  out << "getGoodNumThreads(n,10) with 10 elements per thread returns \"PLUMED_NUM_THREADS\": "
      << (OpenMP::getGoodNumThreads(10*requestedThreads, 10) == requestedThreads ? "true" : "false") << '\n';
  out << "getGoodNumThreads(n,10) with one element less returns 1: "
      << (OpenMP::getGoodNumThreads(10*requestedThreads-1, 10) == 1 ? "true" : "false") << '\n';
  out << "getGoodNumThreads(n,0) is getGoodNumThreads(n,1): "
      << (OpenMP::getGoodNumThreads(requestedThreads, 0) == requestedThreads
          && OpenMP::getGoodNumThreads(requestedThreads-1, 0) == 1 ? "true" : "false") << '\n';
  out << "TimedLoop:\n";
  OpenMP::TimedLoop loop;
  out << "before the calibration runs it returns 1: "
      << (loop.getNumThreads(1000000) == 1 ? "true" : "false") << '\n';
  // two serial runs of 100 elements that took (slightly more than) 1000 microseconds:
  // 10 microseconds per element
  for(unsigned i=0; i<OpenMP::TimedLoop::calibrationRuns; ++i) {
    loop.record(OpenMP::TimedLoop::now()-std::chrono::microseconds(1000),100,1);
  }
  // each thread needs minMicrosecondsPerThread of work
  const std::size_t enough=static_cast<std::size_t>(OpenMP::TimedLoop::minMicrosecondsPerThread/10.0*requestedThreads+0.5);
  out << "with enough work for all the threads it returns \"PLUMED_NUM_THREADS\": "
      << (loop.getNumThreads(enough) == requestedThreads ? "true" : "false") << '\n';
  out << "with less work it returns 1: "
      << (loop.getNumThreads(enough-1) == 1 ? "true" : "false") << '\n';
  out << "After setting the number of threads to 3:\n";
  OpenMP::setNumThreads(3);
  out << "getNumThreads reaturns 3: "
      << (OpenMP::getNumThreads() == 3 ? "true" : "false")
      << "\n";
}