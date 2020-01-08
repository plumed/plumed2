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
#include "PlumedHandle.h"
#include "core/PlumedMain.h"
#include "Tools.h"
#include <cstring>
#ifdef __PLUMED_HAS_DLOPEN
#include <dlfcn.h>
#endif

namespace PLMD
{

PlumedHandle::DlHandle::~DlHandle() {
#ifdef __PLUMED_HAS_DLOPEN
  if(handle) dlclose(handle);
#endif
}

PlumedHandle::PlumedHandle():
  local(new PlumedMain)
{
}

PlumedHandle::PlumedHandle(const char* kernel)
#ifdef __PLUMED_HAS_DLOPEN
  :
  handle([&]() {
  dlerror();
  int mode = RTLD_LOCAL | RTLD_NOW;
#ifdef RTLD_DEEPBIND
// Needed on Linux to avoid namespace clashes
  mode |= RTLD_DEEPBIND;
#endif
  void* h=::dlopen(kernel,mode);
// try to remove the "Kernel" string.
// needed to load old versions
  if(!h) {
    std::string k(kernel);
    auto i=k.rfind("Kernel");
    if(i!=std::string::npos) {
      k=k.substr(0,i) + k.substr(i+6);
      h=::dlopen(k.c_str(),mode);
    }
  }
  plumed_assert(h) << "there was a problem loading kernel "<<kernel <<"\n"<<dlerror();
  return DlHandle(h);
// once the DlHandle has been constructed we know that later exceptions will also call dlclose().
}()),
symbol_((plumed_symbol_table_type*) dlsym(handle,"plumed_symbol_table")),
create_([&]() {
  if(symbol_) {
    plumed_assert(symbol_->functions.create);
    return symbol_->functions.create;
  }
  void* c=nullptr;
  if(!c) c=dlsym(handle,"plumedmain_create");
  if(!c) c=dlsym(handle,"plumed_plumedmain_create");
  plumed_assert(c) << "in kernel "<<kernel<<" I could not find (plumed_)plumedmain_create";
  plumed_create_pointer cc;
  *(void **)(&cc)=c;
  return cc;
}()),
cmd_([&]() {
  if(symbol_) {
    plumed_assert(symbol_->functions.cmd);
    return symbol_->functions.cmd;
  }
  void* c=nullptr;
  if(!c) c=dlsym(handle,"plumedmain_cmd");
  if(!c) c=dlsym(handle,"plumed_plumedmain_cmd");
  plumed_assert(c) << "in kernel "<<kernel<<" I could not find (plumed_)plumedmain_cmd";
  plumed_cmd_pointer cc;
  *(void **)(&cc)=c;
  return cc;
}()),
finalize_([&]() {
  if(symbol_) {
    plumed_assert(symbol_->functions.finalize);
    return symbol_->functions.finalize;
  }
  void* f=nullptr;
  if(!f) f=dlsym(handle,"plumedmain_finalize");
  if(!f) f=dlsym(handle,"plumed_plumedmain_finalize");
  plumed_assert(f) << "in kernel "<<kernel<<" I could not find (plumed_)plumedmain_finalize";
  plumed_finalize_pointer ff;
  *(void **)(&ff)=f;
  return ff;
}()),
p(create_())
// No exceptions thrown past this point.
// Thus, destructor PlumedHandle::~PlumedHandle() will always be called and p will always be finalized.
{}
#else
{
  plumed_error() << "You are trying to dynamically load a kernel, but PLUMED was compiled without dlopen";
}
#endif

PlumedHandle::~PlumedHandle() {
#ifdef __PLUMED_HAS_DLOPEN
  if(p) finalize_(p);
#endif
}

PlumedHandle PlumedHandle::dlopen(const char* path) {
  return PlumedHandle(path);
}

void PlumedHandle::cmd(const char*key,const void*ptr) {
  if(local) local->cmd(key,ptr);
  else if(p && cmd_) cmd_(p,key,ptr);
  else plumed_error() << "should never arrive here (either one or the other should work)";
}

}
