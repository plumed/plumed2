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
  std::vector<double> bigv(cachelineSize);
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
      << (OpenMP::getGoodNumThreads(ptd, cachelineSize) == requestedThreads
              ? "true"
              : "false")
      << '\n';
  out << "After setting the number of threads to 3:\n";
  OpenMP::setNumThreads(3);
  out << "getNumThreads reaturns 3: "
      << (OpenMP::getNumThreads() == 3 ? "true" : "false")
      << "\n";
}