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
#include <iostream>
#include "core/ActionRegister.h"
#include "core/CLToolRegister.h"

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
  auto lockerAction=Register::registrationLock(s);
  void* p=dlopen(s.c_str(),RTLD_NOW|RTLD_LOCAL);
  if(!p) {
    plumed_error()<<"Could not load library "<<s<<"\n"<<dlerror();
  }
  handles.push_back(p);
  Register::completeAllRegistrations(p);
  return p;
#else
  plumed_error()<<"you are trying to use dlopen but it's not configured on your system";
#endif
}

DLLoader::~DLLoader() {
  auto debug=std::getenv("PLUMED_LOAD_DEBUG");
#ifdef __PLUMED_HAS_DLOPEN
  if(debug) {
    std::fprintf(stderr,"delete dlloader\n");
  }
  while(!handles.empty()) {
    int ret=dlclose(handles.back());
    if(ret) {
      std::fprintf(stderr,"+++ error reported by dlclose: %s\n",dlerror());
    }
    handles.pop_back();
  }
  if(debug) {
    std::fprintf(stderr,"end delete dlloader\n");
  }
#endif
}

DLLoader::DLLoader() {
  // do nothing
}

const std::vector<void*> & DLLoader::getHandles() const noexcept {
  return handles;
}

DLLoader::EnsureGlobalDLOpen::EnsureGlobalDLOpen(const void *symbol) noexcept {
#ifdef __PLUMED_HAS_DLOPEN
#ifdef __PLUMED_HAS_DLADDR
  Dl_info info;
  // from the manual:
  // If the address specified in addr could not be matched to a shared
  //      object, then these functions return 0.  In this case, an error
  //      message is not available via dlerror(3).
  int zeroIsError=dladdr(symbol, &info);
  if(zeroIsError!=0) {
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

bool DLLoader::isPlumedGlobal() {
#if defined(__APPLE__)
  bool result;
#ifdef __PLUMED_HAS_DLOPEN
  void* handle;
#if defined(__PLUMED_HAS_RTLD_DEFAULT)
  handle=RTLD_DEFAULT;
#else
  handle=dlopen(NULL,RTLD_LOCAL);
#endif
  // we check for two variants, see wrapper/Plumed.h for an explanation:
  result=dlsym(handle,"plumed_plumedmain_create") || dlsym(handle,"plumedmain_create");
  if(handle) {
    dlclose(handle);
  }
#else
  // if a system cannot use dlopen, we assume plumed is globally available
  result=true;
#endif //__PLUMED_HAS_DLOPEN
  return result;
#else
  plumed_error()<<"DLLoader::isPlumedGlobal() is only functional with APPLE dlsym";
#endif
}

} // namespace PLMD
