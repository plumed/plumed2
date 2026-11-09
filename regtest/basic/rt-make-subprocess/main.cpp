#include "plumed/tools/Subprocess.h"
#include <cerrno>
#include <cstring>
#include <exception>
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <sys/wait.h>
#include <csignal>
#include <cstdlib>

using namespace PLMD;

int main() {
  // if processes are not properly terminated and cleaned
  // only a finite number of them will be allowed

  std::ofstream ofs("should_be_empty");

  // Literal arguments, inherited variables, empty overrides and exit status.
  const std::string literal="space ' $(touch injected_subprocess)";
  const int status=Subprocess::run("/bin/sh", {
    "-c", "test \"$1\" = \"$SUBPROCESS_RUN_VALUE\" && "
    "test \"$SUBPROCESS_RUN_INHERITED\" = inherited && "
    "test \"$SUBPROCESS_RUN_EMPTY\" = '' && exit 7", "test",literal
  }, {{"SUBPROCESS_RUN_VALUE",literal},{"SUBPROCESS_RUN_EMPTY",""}});
  if(status<0 || !WIFEXITED(status) || WEXITSTATUS(status)!=7 ||
      std::string(std::getenv("SUBPROCESS_RUN_VALUE"))!="parent" ||
      std::string(std::getenv("SUBPROCESS_RUN_EMPTY"))!="parent" ||
      std::ifstream("injected_subprocess")) {
    ofs<<"literal argument, environment or exit status failure\n";
  }
  const int signaled=Subprocess::run("/bin/sh", {"-c","kill -TERM $$"});
  if(signaled<0 || !WIFSIGNALED(signaled) || WTERMSIG(signaled)!=SIGTERM) {
    ofs<<"signal status failure\n";
  }
  if(Subprocess::run("/plumed-nonexistent-executable")!=-1 || errno!=ENOENT) {
    ofs<<"missing executable failure\n";
  }
  if(Subprocess::run("/bin/sh", {}, {{"invalid=key","value"}})!=-1 || errno!=EINVAL) {
    ofs<<"invalid environment key failure\n";
  }

  // the number of iteration required to trigger the bugs fixed in
  // 2d45954bbb26471eaf48c72d9a9f0a2dcd8536cc
  // are, on MacOS:
  // ~ 250 for the missing close call
  // ~ 2500 for the missing waitpid call

  // notice that we are not spawning simultaneous processes.
  // they are called in a sequential way and killed (by destructor)
  for(unsigned i=0; i<5000; i++) {
    try {
      Subprocess sp("yes");
    } catch(const std::exception& e) {
      const int saved_errno=errno;
      ofs<<"failed after "<<i<<"\n";
      // Preserve the exception that the reference-file comparison cannot show.
      // errno is recorded at the catch, not assumed to identify the cause.
      std::cerr<<"Subprocess failed after "<<i<<": "<<e.what()
               <<"\nerrno at catch: "<<saved_errno<<" ("<<std::strerror(saved_errno)<<")\n";
      break;
    } catch(...) {
      ofs<<"failed after "<<i<<"\n";
      std::cerr<<"Subprocess failed after "<<i<<": unknown exception\n";
      break;
    }
  }
  return 0;
}

