/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
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
#include "PlumedMainInitializer.h"
#include "PlumedMain.h"
#include "tools/Exception.h"
#include <cstdlib>
#include <cstring>
#if defined __PLUMED_HAS_DLOPEN
#include <dlfcn.h>
#endif

using namespace std;

extern "C" void*plumed_plumedmain_create() {
  return new PLMD::PlumedMain;
}

extern "C" void plumed_plumedmain_cmd(void*plumed,const char*key,const void*val) {
  plumed_massert(plumed,"trying to use a plumed object which is not initialized");
  auto p=static_cast<PLMD::PlumedMain*>(plumed);
  try {
    p->cmd(key,val);
  } catch(std::exception & e) {
// we are at library boundaries.
// if a error_handler was provided, we use it to manage this exception.
// this allows an exception to be catched also if the MD code
// was linked against a different C++ library
    if(!p->callErrorHandler(2,e.what())) throw;
  }
}

extern "C" void plumed_plumedmain_finalize(void*plumed) {
  plumed_massert(plumed,"trying to deallocate a plumed object which is not initialized");
// I think it is not possible to replace this delete with a smart pointer
// since the ownership of this pointer is in a C structure. GB
  delete static_cast<PLMD::PlumedMain*>(plumed);
}

// values here should be consistent with those in plumed_symbol_table_init !!!!
plumed_symbol_table_type plumed_symbol_table=
{1,{plumed_plumedmain_create,plumed_plumedmain_cmd,plumed_plumedmain_finalize}};

// values here should be consistent with those above !!!!
extern "C" void plumed_symbol_table_init() {
  plumed_symbol_table.version=1;
  plumed_symbol_table.functions.create=plumed_plumedmain_create;
  plumed_symbol_table.functions.cmd=plumed_plumedmain_cmd;
  plumed_symbol_table.functions.finalize=plumed_plumedmain_finalize;
}

namespace PLMD {

#define plumed_convert_fptr(ptr,fptr) { ptr=NULL; std::memcpy(&ptr,&fptr,(sizeof(fptr)>sizeof(ptr)?sizeof(ptr):sizeof(fptr))); }

/// Static object which registers Plumed.
/// This is a static object which, during its construction at startup,
/// registers the pointers to plumed_plumedmain_create, plumed_plumedmain_cmd and plumed_plumedmain_finalize
/// to the plumed_kernel_register function.
/// Registration is only required with plumed loader <=2.4, but we do it anyway in order to maintain
/// backward compatibility. Notice that as of plumed 2.5 the plumed_kernel_register is found
/// using dlsym, in order to allow the libplumedKernel library to be loadable also when
/// the plumed_kernel_register symbol is not available.
static class PlumedMainInitializer {
  const bool debug;
public:
  PlumedMainInitializer():
    debug(std::getenv("PLUMED_LOAD_DEBUG"))
  {
// make sure static plumed_function_pointers is initialized here
    plumed_symbol_table_init();
    if(debug) fprintf(stderr,"+++ Initializing PLUMED with plumed_symbol_table version %i at %p\n",plumed_symbol_table.version,(void*)&plumed_symbol_table);
#if defined(__PLUMED_HAS_DLOPEN)
    if(std::getenv("PLUMED_LOAD_SKIP_REGISTRATION")) {
      if(debug) fprintf(stderr,"+++ Skipping registration +++\n");
      return;
    }
    typedef plumed_plumedmain_function_holder* (*plumed_kernel_register_type)(const plumed_plumedmain_function_holder*);
    plumed_kernel_register_type plumed_kernel_register=nullptr;
    void* handle=nullptr;
#if defined(__PLUMED_HAS_RTLD_DEFAULT)
    if(debug) fprintf(stderr,"+++ Registering functions. Looking in RTLD_DEFAULT +++\n");
    void* dls=dlsym(RTLD_DEFAULT,"plumed_kernel_register");
#else
    handle=dlopen(NULL,RTLD_LOCAL);
    if(debug) fprintf(stderr,"+++ Registering functions. dlopen handle at %p +++\n",handle);
    void* dls=dlsym(handle,"plumed_kernel_register");
#endif
    *(void **)(&plumed_kernel_register)=dls;
    if(debug) {
      if(plumed_kernel_register) {
        fprintf(stderr,"+++ plumed_kernel_register found at %p +++\n",dls);
      }
      else fprintf(stderr,"+++ plumed_kernel_register not found +++\n");
    }
    void*createp;
    void*cmdp;
    void*finalizep;
    plumed_convert_fptr(createp,plumed_symbol_table.functions.create);
    plumed_convert_fptr(cmdp,plumed_symbol_table.functions.cmd);
    plumed_convert_fptr(finalizep,plumed_symbol_table.functions.finalize);
    if(plumed_kernel_register && debug) fprintf(stderr,"+++ Registering functions at %p (%p,%p,%p) +++\n",
          (void*)&plumed_symbol_table.functions,createp,cmdp,finalizep);
    if(plumed_kernel_register) (*plumed_kernel_register)(&plumed_symbol_table.functions);
// Notice that handle could be null in the following cases:
// - if we use RTLD_DEFAULT
// - on Linux if we don't use RTLD_DEFAULT, since dlopen(NULL,RTLD_LOCAL) returns a null pointer.
    if(handle) dlclose(handle);
#endif
  }
  ~PlumedMainInitializer() {
    if(debug) fprintf(stderr,"+++ Finalizing PLUMED with plumed_symbol_table at %p\n",(void*)&plumed_symbol_table);
  }
} PlumedMainInitializerRegisterMe;

}


