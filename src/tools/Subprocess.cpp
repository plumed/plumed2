/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2018-2023 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Subprocess.h"
#include "Exception.h"
#include "Tools.h"
#ifdef __PLUMED_HAS_SUBPROCESS
#include <unistd.h>
#include <csignal>
#include <cerrno>
#include <chrono>
#include <thread>
#include <sys/wait.h>
#endif

namespace PLMD {

/// Retrieve PLUMED_ENABLE_SIGNALS.
/// Inline static so that it can store a static variable (for quicker access)
/// without adding a unique global symbol to a library including this header file.
inline static bool SubprocessPidGetenvSignals() noexcept {
  static const bool res=std::getenv("PLUMED_ENABLE_SIGNALS");
  return res;
}

/// Small utility class, used to avoid inclusion of unistd.h> in a header file.
///
/// The subprocess is launched as `/bin/sh -c cmd` and is made leader of a new
/// process group (with pgid equal to its pid). All signals are sent to the whole
/// process group (`kill(-pid,...)`), so that they also reach the processes
/// spawned by the shell (e.g. the python interpreter launched by a script).
/// Sending them only to the shell would leave the actual worker untouched.
class SubprocessPid {
#ifdef __PLUMED_HAS_SUBPROCESS
  /// Time given to the process group to exit after SIGTERM before sending SIGKILL.
  static constexpr std::chrono::milliseconds gracePeriod{1000};
  /// Returns true if the shell process has exited.
  /// The process is NOT reaped (WNOWAIT), so that its pid (and thus the process group id)
  /// cannot be recycled by the system before we are done signaling the group.
  bool exited() const noexcept {
    while(true) {
      siginfo_t info;
      info.si_pid=0;
      if(waitid(P_PID,pid,&info,WEXITED|WNOHANG|WNOWAIT)==0) {
        return info.si_pid==pid;
      }
      if(errno!=EINTR) {
        // ECHILD or other errors: nothing to wait for
        return true;
      }
    }
  }
public:
  const pid_t pid;
  explicit SubprocessPid(pid_t pid_):
    pid(pid_) {
    plumed_assert(pid!=0 && pid!=-1);
  }
  void stop() noexcept {
    // Signals give problems with MPI on Travis.
    // I disable them for now.
    if(SubprocessPidGetenvSignals()) {
      kill(-pid,SIGSTOP);
    }
  }
  void cont() noexcept {
    // Signals give problems with MPI on Travis.
    // I disable them for now.
    if(SubprocessPidGetenvSignals()) {
      kill(-pid,SIGCONT);
    }
  }
  ~SubprocessPid() {
    // The destructor implies we do not need the subprocess anymore.
    // At this point its stdin has already been closed.
    // We first ask the whole process group to terminate with SIGTERM, which
    // allows a clean exit (e.g., shell EXIT traps removing temporary files).
    // SIGCONT is needed because the group might have been stopped with SIGSTOP
    // (see stop()), in which case SIGTERM would remain pending forever.
    kill(-pid,SIGTERM);
    kill(-pid,SIGCONT);
    // Poll with an increasing interval, so that the common case
    // (immediate exit) is fast.
    const auto start=std::chrono::steady_clock::now();
    std::chrono::microseconds interval{50};
    while(!exited() && std::chrono::steady_clock::now()-start<gracePeriod) {
      std::this_thread::sleep_for(interval);
      if(interval<std::chrono::milliseconds(20)) {
        interval*=2;
      }
    }
    // Kill whatever is left in the process group: the shell itself if it ignored
    // SIGTERM, or processes that it spawned and that are still alive.
    // This is safe since the shell has not been reaped yet, so the pgid cannot
    // have been recycled. If nothing is left, kill() just fails with ESRCH.
    kill(-pid,SIGKILL);
    // Reap the shell. This is anyway required to avoid leaks
    int status;
    while(waitpid(pid, &status, 0)<0 && errno==EINTR) {}
  }
#endif
};

Subprocess::Subprocess(const std::string & cmd) {
#ifdef __PLUMED_HAS_SUBPROCESS
  char* arr [] = {
    // const_cast are necessary here due to the declaration of execv
    const_cast<char*>("/bin/sh"),
    const_cast<char*>("-c"),
    const_cast<char*>(cmd.c_str()),
    nullptr
  };
  int cp[2];
  int pc[2];
  if(pipe(pc)<0) {
    plumed_error()<<"error creating parent to child pipe";
  }
  if(pipe(cp)<0) {
    plumed_error()<<"error creating child to parent pipe";
  }
  pid_t forkedPid=fork();
  switch(forkedPid) {
  case -1:
    plumed_error()<<"error forking";
    break;
// CHILD:
  case 0: {
    if(close(1)<0) {
      plumed_error()<<"error closing file";
    }
    if(dup(cp[1])<0) {
      plumed_error()<<"error duplicating file";
    }
    if(close(0)<0) {
      plumed_error()<<"error closing file";
    }
    if(dup(pc[0])<0) {
      plumed_error()<<"error duplicating file";
    }
    if(close(pc[1])<0) {
      plumed_error()<<"error closing file";
    }
    if(close(cp[0])<0) {
      plumed_error()<<"error closing file";
    }
    // Make the child the leader of a new process group, so that signals
    // can be sent to it and to all the processes it spawns.
    if(setpgid(0,0)<0) {
      plumed_error()<<"error setting process group";
    }
    auto err=execv(arr[0],arr);
    plumed_error()<<"error in script file " << cmd << ", execv returned "<<err;
  }
// PARENT::
  default:
    // setpgid is also called here to avoid a race condition:
    // the parent might send a signal to the group before the child has created it.
    // Errors are ignored: the call fails with EACCES if the child has already called execv,
    // but in that case the child has already set its process group.
    setpgid(forkedPid,forkedPid);
    pid=Tools::make_unique<SubprocessPid>(forkedPid);
    if(close(pc[0])<0) {
      plumed_error()<<"error closing file";
    }
    if(close(cp[1])<0) {
      plumed_error()<<"error closing file";
    }
    fpc=pc[1];
    fcp=cp[0];
    fppc=fdopen(fpc,"w");
    parent_to_child.link(fppc);
    fpcp=fdopen(fcp,"r");
    child_to_parent.link(fpcp);
  }
#else
  plumed_error()<<"Subprocess not supported";
#endif
}

Subprocess::~Subprocess() {
#ifdef __PLUMED_HAS_SUBPROCESS
// close files:
  fclose(fppc);
  fclose(fpcp);
// fclose also closes the underlying descriptors,
// so this is not needed:
  close(fpc);
  close(fcp);
// after closing the communication, the subprocess is killed
// in pid's destructor
#endif
}

bool Subprocess::available() noexcept {
#ifdef __PLUMED_HAS_SUBPROCESS
  return true;
#else
  return false;
#endif
}

void Subprocess::stop() noexcept {
#ifdef __PLUMED_HAS_SUBPROCESS
  pid->stop();
#endif
}

void Subprocess::cont() noexcept {
#ifdef __PLUMED_HAS_SUBPROCESS
  pid->cont();
#endif
}

void Subprocess::flush() {
  parent_to_child.flush();
}

Subprocess & Subprocess::getline(std::string & line) {
  child_to_parent.getline(line);
  if(!child_to_parent) {
    plumed_error() <<"error reading subprocess";
  }
  return (*this);
}

Subprocess::Handler::Handler(Subprocess *subp) noexcept:
  sp(subp) {
  subp->cont();
}

Subprocess::Handler::~Handler() {
  if(sp) {
    sp->stop();
  }
}

Subprocess::Handler::Handler(Handler && handler) noexcept :
  sp(handler.sp) {
  handler.sp=nullptr;
}

Subprocess::Handler & Subprocess::Handler::operator=(Handler && handler) noexcept {
  if(this!=&handler) {
    if(sp) {
      sp->stop();
    }
    sp=handler.sp;
    handler.sp=nullptr;
  }
  return *this;
}


}
