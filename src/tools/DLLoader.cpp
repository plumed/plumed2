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
#ifdef __PLUMED_HAS_DLOPEN
#ifdef __PLUMED_HAS_DLADDR
  Dl_info info;
  // from the manual:
  // If the address specified in addr could not be matched to a shared
  //      object, then these functions return 0.  In this case, an error
  //      message is not available via dlerror(3).
  if(dladdr(symbol, &info)!=0) {
    //This "promotes" to GLOBAL the object with the symbol pointed by ptr
    handle_ = dlopen(info.dli_fname, RTLD_GLOBAL | RTLD_NOW);
  } else {
    std::fprintf(stderr,
                 "+++WARNING+++"
                 "Failure in finding any object that contains the symbol %p.\n",
                 symbol);
    
  }
#else
  std::fprintf(stderr,
               "+++WARNING+++"
               "I can't use dladdr for promoting the library containing the symbol %p.\n"
               "This system seems not to support dladdr",
               symbol);
#endif //__PLUMED_HAS_DLADDR
#endif //__PLUMED_HAS_DLOPEN
}

DLLoader::EnsureGlobalDLOpen::~EnsureGlobalDLOpen() {
#ifdef __PLUMED_HAS_DLOPEN
  if (handle_) {
    dlclose(handle_);
  }
#endif //__PLUMED_HAS_DLOPEN
}

} // namespace PLMD
