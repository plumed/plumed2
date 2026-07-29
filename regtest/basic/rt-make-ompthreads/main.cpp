/*
 * Checks the mechanism the MD-code patches rely on to tell PLUMED how many OpenMP threads the engine is
 * using: cmd("setNumOMPthreads") must be reflected by PLMD::OpenMP::getNumThreads().
 *
 * The GROMACS 2025 patch issues this command from the PlumedForceProvider constructor. Without it
 * getNumThreads() reports 1 and every OpenMP-parallel action silently runs serially, which is invisible
 * in results and only shows up as lost performance -- exactly the kind of regression a test should catch.
 *
 * Setting the value twice also pins down that the command is re-issuable rather than latch-once.
 */
#include "plumed/wrapper/Plumed.h"
#include "plumed/tools/OpenMP.h"
#include <fstream>

int main() {
  std::ofstream ofs("output");

  // Query first so any PLUMED_NUM_THREADS in the environment has already been latched: the point is that
  // an explicit setNumOMPthreads overrides whatever the default resolved to.
  PLMD::OpenMP::getNumThreads();

  PLMD::Plumed p;

  unsigned nt = 3;
  p.cmd("setNumOMPthreads", &nt);
  ofs << "after setNumOMPthreads(3): " << PLMD::OpenMP::getNumThreads() << "\n";

  nt = 5;
  p.cmd("setNumOMPthreads", &nt);
  ofs << "after setNumOMPthreads(5): " << PLMD::OpenMP::getNumThreads() << "\n";

  // 0 threads is meaningless; PlumedMain maps it to 1 rather than passing it through.
  nt = 0;
  p.cmd("setNumOMPthreads", &nt);
  ofs << "after setNumOMPthreads(0): " << PLMD::OpenMP::getNumThreads() << "\n";

  return 0;
}
