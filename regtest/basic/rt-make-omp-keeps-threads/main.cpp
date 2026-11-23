/*
 * PLUMED must not destroy the OpenMP threads of the code that calls it.
 *
 * With the GNU OpenMP runtime, a parallel region with fewer threads than the previous one
 * (but more than one) destroys the threads that are not needed, and a later, larger region
 * creates new ones. An MD code assigns a core to each of its threads once, at the start:
 * the new threads do not have it, they all run on the core of the calling thread, and the
 * MD code becomes many times slower (seen with GROMACS, 4 threads, and the cells of a list
 * of 2200 atoms built by 2 threads).
 * So PLUMED uses all the threads or one, see PLMD::OpenMP::getGoodNumThreads.
 *
 * Here the caller has 4 threads, PLUMED computes a COORDINATION with a neighbor list of
 * 3000 atoms (cells, lists, copies of the atoms), and the threads of the caller must be the
 * same after it. With other runtimes, or without OpenMP, the test passes in any case.
 */
#include "plumed/wrapper/Plumed.h"
#include <fstream>
#include <set>
#include <vector>
#ifdef __linux__
#include <unistd.h>
#include <sys/syscall.h>
#endif

// The identifiers of the threads of a region with 4 threads. They are the ones of the
// operating system: those of the C++ library are reused by the threads created later.
// On other systems the test only checks that PLUMED runs.
static std::set<long> team() {
  std::set<long> ids;
  #pragma omp parallel num_threads(4)
  {
#ifdef __linux__
    const long id=static_cast<long>(syscall(SYS_gettid));
    #pragma omp critical
    ids.insert(id);
#endif
  }
  return ids;
}

int main() {
  std::ofstream ofs("output");
  const std::set<long> before=team();

  const int natoms=3000;
  const double side=4.0;
  std::vector<double> positions(3*natoms), forces(3*natoms,0.0), masses(natoms,1.0), charges(natoms,0.0);
  // a deterministic, roughly uniform set of positions
  unsigned long long seed=12345;
  for(auto & x : positions) {
    seed=seed*6364136223846793005ULL+1442695040888963407ULL;
    x=side*static_cast<double>(seed>>11)/9007199254740992.0;
  }
  double box[9]= {side,0,0,0,side,0,0,0,side};
  double virial[9]= {0,0,0,0,0,0,0,0,0};
  {
    PLMD::Plumed p;
    int nt=4;
    p.cmd("setNumOMPthreads",&nt);
    p.cmd("setNatoms",&natoms);
    p.cmd("setLogFile","plumed.log");
    p.cmd("init");
    p.cmd("readInputLine","c: COORDINATION GROUPA=1-500 GROUPB=501-3000 R_0=0.3 NLISTCELLS NL_CUTOFF=0.6 NL_STRIDE=1");
    p.cmd("readInputLine","RESTRAINT ARG=c AT=0 SLOPE=1");
    for(int step=0; step<3; ++step) {
      p.cmd("setStep",&step);
      p.cmd("setPositions",positions.data());
      p.cmd("setMasses",masses.data());
      p.cmd("setCharges",charges.data());
      p.cmd("setBox",&box[0]);
      p.cmd("setForces",forces.data());
      p.cmd("setVirial",&virial[0]);
      p.cmd("calc");
    }
  }
  const std::set<long> after=team();
  ofs<<"the threads of the caller are the same after PLUMED: "<<(before==after ? "yes" : "no")<<"\n";
  return 0;
}
