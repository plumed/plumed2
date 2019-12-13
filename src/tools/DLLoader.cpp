/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2019 The plumed team
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
  if(debug) std::fprintf(stderr,"delete dlloader\n");
#ifdef __PLUMED_HAS_DLOPEN
  while(!handles.empty()) {
    int ret=dlclose(handles.top());
    if(ret) {
      std::fprintf(stderr,"+++ error reported by dlclose: %s\n",dlerror());
    }
    handles.pop();
  }
#endif
  if(debug) std::fprintf(stderr,"end delete dlloader\n");
}

DLLoader::DLLoader() {
  // do nothing
}


}
