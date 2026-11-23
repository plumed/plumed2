/*
 * PLUMED shares the OpenMP threads of the MD code. With the GNU OpenMP runtime, regions with a number of
 * threads different from that of the MD code (other than 1) replace threads of the MD code at every step
 * (see PLMD::OpenMP::getGoodNumThreads). PLUMED_NUM_THREADS takes precedence over setNumOMPthreads, so
 * when the two differ the log has a warning.
 * Here the MD code sets 4 threads and PLUMED_NUM_THREADS=2 (in config). Without OpenMP there is no warning.
 */
#include "plumed/wrapper/Plumed.h"
#include <fstream>
#include <sstream>
#include <string>

int main() {
  {
    PLMD::Plumed p;
    unsigned nt=4;
    p.cmd("setNumOMPthreads",&nt);
    p.cmd("setLogFile","plumed.log");
    p.cmd("init");
  }
  std::ifstream log("plumed.log");
  std::stringstream ss;
  ss<<log.rdbuf();
  const bool warned=ss.str().find("WARNING: PLUMED uses 2 OpenMP threads, but the MD code set 4")!=std::string::npos;
#ifdef _OPENMP
  const bool expected=true;
#else
  const bool expected=false;
#endif
  std::ofstream ofs("output");
  ofs<<"warning as expected: "<<(warned==expected ? "yes" : "no")<<"\n";
  return 0;
}
