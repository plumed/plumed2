/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2021 The plumed team
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
#include "Exception.h"

#if defined(__PLUMED_HAS_EXECINFO)
#include <execinfo.h>
#endif

#include <cstdio>
#include <cstring>
#include <cstdlib>

namespace PLMD {

Exception::Exception()
{
  callstack.fill(nullptr);
#ifdef __PLUMED_HAS_EXECINFO
  callstack_n = backtrace(&callstack[0], callstack.size()-1);
  const char* env=std::getenv("PLUMED_STACK_TRACE");
  if(env && !std::strcmp(env,"yes")) {
    msg+="\n\n********** STACK DUMP **********\n";
    msg+=stack();
    msg+="\n********** END STACK DUMP **********\n";
  }
#endif
  msg+="\n+++ PLUMED error";
}

Exception& Exception::operator<<(const std::string&msg)
{
  if(msg.length()>0) {
    if(note) this->msg +="\n+++ message follows +++\n";
    this->msg +=msg;
    note=false;
  }
  return *this;
}

Exception& Exception::operator<<(const Location&loc)
{
  if(loc.file) {
    char cline[1000];
    std::sprintf(cline,"%u",loc.line);
    this->msg += "\n+++ at ";
    this->msg += loc.file;
    this->msg += ":";
    this->msg += cline;
    if(loc.pretty && loc.pretty[0]) {
      this->msg += ", function ";
      this->msg += loc.pretty;
    }
  }
  note=true;
  return *this;
}

Exception& Exception::operator<<(const Assertion&as)
{
  if(as.assertion) {
    this->msg += "\n+++ assertion failed: ";
    this->msg += as.assertion;
  }
  note=true;
  return *this;
}

const char* Exception::stack() const {
#ifdef __PLUMED_HAS_EXECINFO
  if(stackTrace.length()==0) {
    char** strs = backtrace_symbols(&callstack[0], callstack_n);
    for (int i = 0; i < callstack_n; ++i) {stackTrace+=strs[i]; stackTrace+="\n";}
    free(strs);
  }
#endif
  return stackTrace.c_str();
}

}


