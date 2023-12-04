/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "DLLoader.h"
#include "Exception.h"

#include <cstdlib>

#ifdef __PLUMED_HAS_DLOPEN
#include <dlfcn.h>
#endif

namespace PLMD {

bool DLLoader::installed() {
#ifdef __PLUMED_HAS_DLOPEN
  return true;
#else
  return false;
#endif
}


void* DLLoader::load(const std::string&s) {
#ifdef __PLUMED_HAS_DLOPEN
  void* p=dlopen(s.c_str(),RTLD_NOW|RTLD_LOCAL);
  if(!p) {
    lastError=dlerror();
  } else {
    lastError="";
    handles.push(p);
  }
  return p;
#else
  return NULL;
#endif
}

const std::string & DLLoader::error() {
  return lastError;
}

DLLoader::~DLLoader() {
  auto debug=std::getenv("PLUMED_LOAD_DEBUG");
#ifdef __PLUMED_HAS_DLOPEN
  if(debug) std::fprintf(stderr,"delete dlloader\n");
  while(!handles.empty()) {
    int ret=dlclose(handles.top());
    if(ret) {
      std::fprintf(stderr,"+++ error reported by dlclose: %s\n",dlerror());
    }
    handles.pop();
  }
  if(debug) std::fprintf(stderr,"end delete dlloader\n");
#endif
}

DLLoader::DLLoader() {
  // do nothing
}

DLLoader::EnsureGlobalDLOpen::EnsureGlobalDLOpen(const void *symbol) noexcept {

  // dladddr might be not available
  // this part should likely be protected in a __PLUMED_HAS_DLADDR ifdef
#ifdef __PLUMED_HAS_DLOPEN
  Dl_info info;
  int zeroIsError=dladdr(symbol, &info);
  // from the manual:
  // If the address specified in addr could not be matched to a shared
  //      object, then these functions return 0.  In this case, an error
  //      message is not available via dlerror(3).
  if(zeroIsError==0) {
    plumed_error() << "Failure in finding any object that contains the symbol "<< symbol;
  }
  // std::cerr << "Path: " << info.dli_fname << "\n";
  //This "promotes" to GLOBAL the object with the symbol pointed by ptr
  handle_ = dlopen(info.dli_fname, RTLD_GLOBAL | RTLD_NOW);
  // std::cerr << "Handle: " << handle << "\n";
#else
  //errormessage
#endif

}

DLLoader::EnsureGlobalDLOpen::~EnsureGlobalDLOpen() {
#ifdef __PLUMED_HAS_DLOPEN
  if (handle_) {
    dlclose(handle_);
  }
#endif //__PLUMED_HAS_DLOPEN
}


}
