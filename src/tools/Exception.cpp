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

// Stack traces: execinfo (found by ./configure) or, with MSVC on Windows,
// CaptureStackBackTrace and DbgHelp for the names
#if !defined(__PLUMED_HAS_EXECINFO) && defined(_MSC_VER)
#define PLUMED_EXCEPTION_DBGHELP 1
#endif

#if defined(__PLUMED_HAS_EXECINFO)
#include <execinfo.h>
#elif defined(PLUMED_EXCEPTION_DBGHELP)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#include <dbghelp.h>
#include <mutex>
// also for the programs linking the static kernel
#pragma comment(lib, "dbghelp.lib")
#endif

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <new>
#include <vector>

namespace PLMD {

namespace {
// see https://www.geeksforgeeks.org/simplify-directory-path-unix-like/

// function to simplify a Unix - styled
// absolute path
std::string simplify(const std::string & path) {
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
      if (!v.empty()) {
        v.pop_back();
      }
    } else if (dir == "." || dir == "") {
      // do nothing (added for better understanding.)
    } else {
      // push the current directory into the vector.
      v.push_back(dir);
    }
  }

  // forming the ans
  bool first=true;
  for (const auto & i : v) {
    if(!first) {
      ans += "/";
    }
    first=false;
    ans += i;
  }

  // vector is empty
  if (ans == "") {
    return "/";
  }

  return ans;
}

}

namespace {
/// A copy of the stack trace, or nullptr (also if it cannot be allocated)
std::array<void*,128>* copyCallstack(const std::unique_ptr<std::array<void*,128>> & callstack) {
  return callstack ? new(std::nothrow) std::array<void*,128>(*callstack) : nullptr;
}

#ifdef PLUMED_EXCEPTION_DBGHELP
/// Load the symbols of the process.
/// Names and lines are found only if the debug information (pdb) is available.
bool initializeSymbols(HANDLE process) {
  SymSetOptions(SymGetOptions() | SYMOPT_UNDNAME | SYMOPT_DEFERRED_LOADS | SYMOPT_LOAD_LINES);
  return SymInitialize(process,nullptr,TRUE)!=FALSE;
}
#endif
}

Exception::Exception() {
#if defined(__PLUMED_HAS_EXECINFO) || defined(PLUMED_EXCEPTION_DBGHELP)
  // no exception here if there is no memory: no stack trace in that case
  callstack.reset(new(std::nothrow) std::array<void*,128>);
  if(callstack) {
    callstack->fill(nullptr);
#ifdef __PLUMED_HAS_EXECINFO
    callstack_n = backtrace(callstack->data(), callstack->size()-1);
#else
    callstack_n = CaptureStackBackTrace(0,static_cast<DWORD>(callstack->size()-1),callstack->data(),nullptr);
#endif
  }
  const char* env=std::getenv("PLUMED_STACK_TRACE");
  if(env && !std::strcmp(env,"yes")) {
    msg+="\n\n********** STACK DUMP **********\n";
    msg+=stack();
    msg+="\n********** END STACK DUMP **********\n";
  }
#endif
}

Exception::Exception(const Exception & e):
  msg(e.msg),
  note(e.note),
  callstack(copyCallstack(e.callstack)),
  callstack_n(callstack ? e.callstack_n : 0),
  stackTrace(e.stackTrace) {
}

Exception & Exception::operator=(const Exception & e) {
  if(this!=&e) {
    msg=e.msg;
    note=e.note;
    callstack.reset(copyCallstack(e.callstack));
    callstack_n=callstack ? e.callstack_n : 0;
    stackTrace=e.stackTrace;
  }
  if(stream) {
    stream->str("");
  }
  return *this;
}

std::stringstream & Exception::getStream() {
  if(!stream) {
    stream=std::make_unique<std::stringstream>();
  }
  return *stream;
}

const std::array<void*,128> & Exception::trace() const noexcept {
  static const std::array<void*,128> empty{};
  return callstack ? *callstack : empty;
}

Exception& Exception::operator<<(const std::string&addmsg) {
  if(addmsg.length()>0) {
    if(note) {
      msg +="\n";
    }
    msg +=addmsg;
    note=false;
  }
  return *this;
}

Exception& Exception::operator<<(const Location&loc) {
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

Exception& Exception::operator<<(const Assertion&as) {
  if(as.assertion) {
    this->msg += "\n+++ assertion failed: ";
    this->msg += as.assertion;
  }
  note=true;
  return *this;
}

const char* Exception::stack() const {
#ifdef __PLUMED_HAS_EXECINFO
  if(stackTrace.length()==0 && callstack) {
    char** strs = backtrace_symbols(callstack->data(), callstack_n);
    for (int i = 0; i < callstack_n; ++i) {
      stackTrace+=strs[i];
      stackTrace+="\n";
    }
    free(strs);
  }
#elif defined(PLUMED_EXCEPTION_DBGHELP)
  if(stackTrace.length()==0 && callstack) {
    // the DbgHelp functions are not thread safe
    static std::mutex mutex;
    std::lock_guard<std::mutex> lock(mutex);
    const HANDLE process=GetCurrentProcess();
    static const bool symbols=initializeSymbols(process);
    for(int i=0; i<callstack_n; ++i) {
      const DWORD64 address=reinterpret_cast<DWORD64>((*callstack)[i]);
      if(!address) {
        continue;
      }
      // module(function+offset) [address] file:line, as backtrace_symbols
      char module[MAX_PATH]="?";
      HMODULE handle=nullptr;
      if(GetModuleHandleExA(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS | GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
                            reinterpret_cast<LPCSTR>((*callstack)[i]),&handle)) {
        GetModuleFileNameA(handle,module,MAX_PATH);
      }
      const char* base=std::strrchr(module,'\\');
      stackTrace+=base ? base+1 : module;
      char line[64];
      if(symbols) {
        alignas(SYMBOL_INFO) char buffer[sizeof(SYMBOL_INFO)+256];
        auto info=reinterpret_cast<SYMBOL_INFO*>(buffer);
        info->SizeOfStruct=sizeof(SYMBOL_INFO);
        info->MaxNameLen=256;
        DWORD64 displacement=0;
        if(SymFromAddr(process,address,&displacement,info)) {
          std::snprintf(line,sizeof(line),"+0x%llx",static_cast<unsigned long long>(displacement));
          stackTrace+="(";
          stackTrace+=info->Name;
          stackTrace+=line;
          stackTrace+=")";
        }
      }
      std::snprintf(line,sizeof(line)," [0x%llx]",static_cast<unsigned long long>(address));
      stackTrace+=line;
      if(symbols) {
        IMAGEHLP_LINE64 source{};
        source.SizeOfStruct=sizeof(IMAGEHLP_LINE64);
        DWORD column=0;
        if(SymGetLineFromAddr64(process,address,&column,&source) && source.FileName) {
          std::snprintf(line,sizeof(line),":%lu",static_cast<unsigned long>(source.LineNumber));
          stackTrace+=" ";
          stackTrace+=source.FileName;
          stackTrace+=line;
        }
      }
      stackTrace+="\n";
    }
  }
#endif
  return stackTrace.c_str();
}

}


