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
#include "PlumedHandle.h"
#include "core/PlumedMain.h"
#include "Tools.h"
#include "lepton/Exception.h"
#include <cstring>
#ifdef __PLUMED_HAS_DLOPEN
#include <dlfcn.h>
#endif

// Including Plumed.h in this manner allows to create a local
// implementation of the wrapper in an anonymous namespace.
// This allows to avoid recoding all the Plumed.h stuff here
// and at the same time avoids possible conflicts.
#define __PLUMED_WRAPPER_IMPLEMENTATION 1
#define __PLUMED_WRAPPER_EXTERN 0
#define __PLUMED_WRAPPER_CXX_ANONYMOUS_NAMESPACE 1
#define __PLUMED_WRAPPER_CXX_ANONYMOUS_NAMESPACE_PLMD_EXCEPTIONS 1
#include "../wrapper/Plumed.h"

namespace PLMD
{


PlumedHandle::PlumedHandle():
  local(Tools::make_unique<PlumedMain>())
{
}

PlumedHandle::PlumedHandle(const char* kernel)
#ifdef __PLUMED_HAS_DLOPEN
  :
  loaded(plumed_c2v(plumed_create_dlopen(kernel)))
{
  if(!plumed_valid(plumed_v2c(loaded))) {
    // this is necessary to make sure loaded is properly destroyed
    plumed_finalize(plumed_v2c(loaded));
    plumed_error() << "You are trying to dynamically load a kernel, but the path " << kernel <<" could not be opened";
  }
}
#else
{
  plumed_error() << "You are trying to dynamically load a kernel, but PLUMED was compiled without dlopen";
}
#endif

PlumedHandle::~PlumedHandle() {
  if(loaded) plumed_finalize(plumed_v2c(loaded));
}

PlumedHandle PlumedHandle::dlopen(const char* path) {
  return PlumedHandle(path);
}

void PlumedHandle::cmd(const std::string & key,const TypesafePtr & ptr) {
  if(local) {
    local->cmd(key.c_str(),ptr);
  } else if(loaded) {
    plumed_safeptr safe;
    safe.ptr=ptr.getRaw();
    safe.nelem=ptr.getNelem();
    safe.shape=const_cast<std::size_t*>(ptr.getShape());
    safe.flags=ptr.getFlags();
    safe.opt=nullptr;
    plumed_cmd(plumed_v2c(loaded),key.c_str(),safe);
  } else plumed_error() << "should never arrive here (either one or the other should work)";
}

}
