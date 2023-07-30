/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include <vector>

namespace PLMD {

namespace {
// see https://www.geeksforgeeks.org/simplify-directory-path-unix-like/

// function to simplify a Unix - styled
// absolute path
std::string simplify(const std::string & path)
{
  // using vector in place of stack
  std::vector<std::string> v;
  int n = path.length();
  std::string ans;
  for (int i = 0; i < n; i++) {
    std::string dir = "";
    // forming the current directory.
    while (i < n && path[i] != '/') {
      dir += path[i];
      i++;
    }

    // if ".." , we pop.
    if (dir == "..") {
      if (!v.empty())
        v.pop_back();
    }
    else if (dir == "." || dir == "") {
      // do nothing (added for better understanding.)
    }
    else {
      // push the current directory into the vector.
      v.push_back(dir);
    }
  }

  // forming the ans
  bool first=true;
  for (auto i : v) {
    if(!first) ans += "/";
    first=false;
    ans += i;
  }

  // vector is empty
  if (ans == "")
    return "/";

  return ans;
}

}

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
}

Exception& Exception::operator<<(const std::string&msg)
{
  if(msg.length()>0) {
    if(note) this->msg +="\n";
    this->msg +=msg;
    note=false;
  }
  return *this;
}

Exception& Exception::operator<<(const Location&loc)
{
  if(loc.file) {
    const std::size_t clinelen=1000;
    char cline[clinelen];
    std::snprintf(cline,clinelen,"%u",loc.line);
    this->msg += "\n(";
    try {
      this->msg += simplify(loc.file);
    } catch(...) {
      // ignore
    }
    this->msg += ":";
    this->msg += cline;
    this->msg += ")";
    if(loc.pretty && loc.pretty[0]) {
      this->msg += " ";
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


