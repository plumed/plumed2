/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2019 The plumed team
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
#endif

using namespace std;
namespace PLMD {

Subprocess::Subprocess(const std::string & cmd) {
#ifdef __PLUMED_HAS_SUBPROCESS
  char* arr[4];
  arr[0]=const_cast<char*>("/bin/sh");
  arr[1]=const_cast<char*>("-c");
  arr[2]=const_cast<char*>(cmd.c_str());
  arr[3]=nullptr;
  int cp[2];
  int pc[2];
  if(pipe(pc)<0) plumed_error()<<"error creating parent to child pipe";
  if(pipe(cp)<0) plumed_error()<<"error creating child to parent pipe";
  int pid=fork();
  switch(pid) {
  case -1:
    plumed_error()<<"error forking";
    break;
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
  default:
    if(close(pc[0])<0) plumed_error()<<"error closing file";
    if(close(cp[1])<0) plumed_error()<<"error closing file";
    fpc=pc[1];
    fcp=cp[0];
    parent_to_child.link(fdopen(fpc,"w"));
    child_to_parent.link(fdopen(fcp,"r"));
  }
#else
  plumed_error()<<"Subprocess not supported";
#endif
}

Subprocess::~Subprocess() {
#ifdef __PLUMED_HAS_SUBPROCESS
// fpc should be closed to terminate the child executable
  if(close(fpc)<0) plumed_error()<<"error closing file";
// fcp should not be closed because it could make the child executable fail
#endif
}

bool Subprocess::available() {
#ifdef __PLUMED_HAS_SUBPROCESS
  return true;
#else
  return false;
#endif
}

}
