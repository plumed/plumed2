/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2018-2020 The plumed team
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
class SubprocessPid {
#ifdef __PLUMED_HAS_SUBPROCESS
public:
  pid_t pid;
  explicit SubprocessPid(pid_t pid):
    pid(pid)
  {
    plumed_assert(pid!=0 && pid!=-1);
  }
  void stop() noexcept {
    // Signals give problems with MPI on Travis.
    // I disable them for now.
    if(SubprocessPidGetenvSignals()) if(pid!=0 && pid!=-1) kill(pid,SIGSTOP);
  }
  void cont() noexcept {
    // Signals give problems with MPI on Travis.
    // I disable them for now.
    if(SubprocessPidGetenvSignals()) if(pid!=0 && pid!=-1) kill(pid,SIGCONT);
  }
  ~SubprocessPid() {
    // this is apparently working also with MPI on Travis.
    if(pid!=0 && pid!=-1) kill(pid,SIGINT);
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
  if(pipe(pc)<0) plumed_error()<<"error creating parent to child pipe";
  if(pipe(cp)<0) plumed_error()<<"error creating child to parent pipe";
  pid_t pid=fork();
  switch(pid) {
  case -1:
    plumed_error()<<"error forking";
    break;
// CHILD:
  case 0:
  {
    if(close(1)<0) plumed_error()<<"error closing file";
    if(dup(cp[1])<0) plumed_error()<<"error duplicating file";
    if(close(0)<0) plumed_error()<<"error closing file";
    if(dup(pc[0])<0) plumed_error()<<"error duplicating file";
    if(close(pc[1])<0) plumed_error()<<"error closing file";
    if(close(cp[0])<0) plumed_error()<<"error closing file";
    execv(arr[0],arr);
    plumed_error()<<"error in script file";
  }
// PARENT::
  default:
    this->pid=Tools::make_unique<SubprocessPid>(pid);
    if(close(pc[0])<0) plumed_error()<<"error closing file";
    if(close(cp[1])<0) plumed_error()<<"error closing file";
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
// fpc should be closed to terminate the child executable
  fclose(fppc);
  close(fpc);
// fcp should not be closed because it could make the child executable fail
/// TODO: check if this is necessary and make this class exception safe!
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
  if(!child_to_parent) plumed_error() <<"error reading subprocess";
  return (*this);
}

Subprocess::Handler::Handler(Subprocess *sp) noexcept:
  sp(sp)
{
  sp->cont();
}

Subprocess::Handler::~Handler() {
  if(sp) sp->stop();
}

Subprocess::Handler::Handler(Handler && handler) noexcept :
  sp(handler.sp)
{
  handler.sp=nullptr;
}

Subprocess::Handler & Subprocess::Handler::operator=(Handler && handler) noexcept {
  if(this!=&handler) {
    if(sp) sp->stop();
    sp=handler.sp;
    handler.sp=nullptr;
  }
  return *this;
}


}
