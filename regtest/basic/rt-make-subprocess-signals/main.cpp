#include "plumed/tools/Subprocess.h"
#include <cerrno>
#include <chrono>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <thread>
#include <unistd.h>

using namespace PLMD;

// Returns true if the process does not exist anymore (or is a zombie waiting
// to be reaped by init, which is equivalent for our purposes).
static bool gone(pid_t p) {
  for(unsigned i=0; i<500; i++) {
    if(kill(p,0)<0 && errno==ESRCH) {
      return true;
    }
    std::string cmd="ps -o stat= -p "+std::to_string(p)+" 2>/dev/null";
    FILE* f=popen(cmd.c_str(),"r");
    char buf[16]= {};
    bool zombie=false;
    if(f) {
      if(std::fgets(buf,sizeof(buf),f)) {
        zombie=(buf[0]=='Z' || buf[1]=='Z');
      }
      pclose(f);
    }
    if(zombie) {
      return true;
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
  }
  return false;
}

static bool exists(const std::string & path) {
  std::ifstream ifs(path);
  return ifs.good();
}

int main() {
  // enable SIGSTOP/SIGCONT in Subprocess::stop()/cont()
  setenv("PLUMED_ENABLE_SIGNALS","1",1);

  std::ofstream ofs("output");

  // 1. A bash script that spawns a grandchild and has an EXIT trap,
  //    similar to what scripts/selector.sh does.
  //    The subprocess is stopped before being destroyed.
  //    Signals must reach the grandchild, the destructor should not hang
  //    on the stopped process group, and the trap should be executed.
  {
    std::remove("trap_executed");
    std::string line;
    {
      Subprocess sp("bash -c 'trap \"echo yes > trap_executed\" EXIT; sleep 1000 & echo $!; wait'");
      sp.getline(line);
      sp.stop();
    }
    pid_t grandchild=std::stoi(line);
    ofs<<"test1 grandchild terminated: "<<(gone(grandchild)?"yes":"no")<<"\n";
    ofs<<"test1 trap executed: "<<(exists("trap_executed")?"yes":"no")<<"\n";
  }

  // 2. Same, but stop() is called immediately after construction,
  //    possibly before the child process is running.
  {
    std::remove("trap_executed");
    {
      Subprocess sp("bash -c 'trap \"echo yes > trap_executed\" EXIT; echo ready; while true; do sleep 0.1; done'");
      sp.stop();
      // the process is stopped, so we resume it to read the first line
      auto h=sp.contStop();
      std::string line;
      sp.getline(line);
      ofs<<"test2 read: "<<line<<"\n";
    }
    ofs<<"test2 trap executed: "<<(exists("trap_executed")?"yes":"no")<<"\n";
  }

  // 3. A process group ignoring SIGTERM must be killed anyway (with SIGKILL)
  {
    std::string line;
    {
      Subprocess sp("trap '' TERM; echo $$; while true; do sleep 0.1; done");
      sp.getline(line);
    }
    pid_t shell=std::stoi(line);
    ofs<<"test3 shell terminated: "<<(gone(shell)?"yes":"no")<<"\n";
  }

  return 0;
}
