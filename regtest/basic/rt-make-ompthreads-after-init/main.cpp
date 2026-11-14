/*
 * The MD code tells PLUMED the number of OpenMP threads with cmd("setNumOMPthreads"). The actions size their
 * per-thread storage when they are created, so once a PLUMED instance has been initialized the number of
 * threads can only decrease: a larger value is ignored (with a warning) instead of making the actions write
 * outside their buffers. The number of threads is shared by all the PLUMED instances in the process, so the
 * limit holds as long as an initialized instance exists.
 *
 * Each run computes 100 DISTANCEs -> CUSTOM -> SUM -> RESTRAINT on 200 atoms (enough tasks to use 4 threads
 * in the task manager) and the forces are compared with those of a run with 4 threads.
 */
#include "plumed/wrapper/Plumed.h"
#include "plumed/tools/OpenMP.h"
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

static const int natoms = 200;

static void setup(PLMD::Plumed& p) {
  p.cmd("setNatoms", natoms);
  p.cmd("setMDEngine", "test");
  p.cmd("setTimestep", 1.0);
  p.cmd("setLogFile", "/dev/null");
  p.cmd("init");
  std::string line = "d: DISTANCE";
  for (int i = 0; i < natoms / 2; ++i) {
    line += " ATOMS" + std::to_string(i + 1) + "=" + std::to_string(2 * i + 1) + "," + std::to_string(2 * i + 2);
  }
  p.cmd("readInputLine", line.c_str());
  p.cmd("readInputLine", "c: CUSTOM ARG=d FUNC=x*x PERIODIC=NO");
  p.cmd("readInputLine", "s: SUM ARG=c PERIODIC=NO");
  p.cmd("readInputLine", "RESTRAINT ARG=s AT=0 KAPPA=1");
}

static std::vector<double> calc(PLMD::Plumed& p) {
  std::vector<double> pos(3 * natoms), forces(3 * natoms, 0.0), masses(natoms, 1.0), box(9, 0.0), virial(9, 0.0);
  for (int i = 0; i < natoms; ++i) {
    pos[3 * i + 0] = 0.37 * i;
    pos[3 * i + 1] = std::sin(0.1 * i);
    pos[3 * i + 2] = std::cos(0.07 * i);
  }
  for (int step = 0; step < 2; ++step) {
    std::fill(forces.begin(), forces.end(), 0.0);
    std::fill(virial.begin(), virial.end(), 0.0);
    p.cmd("setStep", step);
    p.cmd("setBox", box.data(), {3, 3});
    p.cmd("setMasses", masses.data(), {natoms});
    p.cmd("setPositions", pos.data(), {natoms, 3});
    p.cmd("setForces", forces.data(), {natoms, 3});
    p.cmd("setVirial", virial.data(), {3, 3});
    p.cmd("calc");
  }
  return forces;
}

static bool same(const std::vector<double>& a, const std::vector<double>& b) {
  double maxdiff = 0.0, maxf = 0.0;
  for (unsigned i = 0; i < a.size(); ++i) {
    maxdiff = std::max(maxdiff, std::fabs(a[i] - b[i]));
    maxf = std::max(maxf, std::fabs(a[i]));
  }
  return maxf > 0 && maxdiff <= 1e-10 * (1 + maxf);
}

static void setThreads(PLMD::Plumed& p, unsigned nt) {
  p.cmd("setNumOMPthreads", &nt);
}

int main() {
  // read PLUMED_NUM_THREADS (if set) now, so that it does not override the values set below
  PLMD::OpenMP::getNumThreads();
  std::ofstream ofs("output");

  std::vector<double> ref;
  {
    PLMD::Plumed p;
    setThreads(p, 4);
    setup(p);
    ref = calc(p);
  }

  // more threads after init: ignored
  {
    PLMD::Plumed p;
    setThreads(p, 1);
    setup(p);
    setThreads(p, 4);
    ofs << "increase after init, threads: " << PLMD::OpenMP::getNumThreads() << "\n";
    ofs << "increase after init, same forces: " << (same(ref, calc(p)) ? "yes" : "no") << "\n";
  }

  // fewer threads after init: accepted
  {
    PLMD::Plumed p;
    setThreads(p, 4);
    setup(p);
    setThreads(p, 2);
    ofs << "decrease after init, threads: " << PLMD::OpenMP::getNumThreads() << "\n";
    // an action created now sizes its storage for 2 threads: going back to more threads is not allowed
    p.cmd("readInputLine", "c2: CUSTOM ARG=d FUNC=2*x PERIODIC=NO");
    setThreads(p, 3);
    ofs << "increase after a decrease, threads: " << PLMD::OpenMP::getNumThreads() << "\n";
    ofs << "decrease after init, same forces: " << (same(ref, calc(p)) ? "yes" : "no") << "\n";
  }

  // another instance cannot raise the number of threads while an initialized instance exists
  {
    PLMD::Plumed b;
    {
      PLMD::Plumed a;
      setThreads(a, 2);
      setup(a);
      setThreads(b, 4);
      ofs << "other instance, threads: " << PLMD::OpenMP::getNumThreads() << "\n";
      ofs << "other instance, same forces: " << (same(ref, calc(a)) ? "yes" : "no") << "\n";
    }
    // once the initialized instance is gone the limit is released
    setThreads(b, 4);
    ofs << "after destruction, threads: " << PLMD::OpenMP::getNumThreads() << "\n";
  }
  return 0;
}
