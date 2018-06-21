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
#include "PlumedMain.h"
#include "tools/Exception.h"
#include <cstdlib>
#if defined __PLUMED_HAS_DLOPEN
#include <dlfcn.h>
#endif

using namespace std;

// !!!!!!!!!!!!!!!!!!!!!!    DANGER   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
// THE FOLLOWING ARE UTILITIES WHICH ARE NECESSARY FOR DYNAMIC LOADING OF THE PLUMED KERNEL:
// This section should be consistent with the Plumed.h file.
// Since the Plumed.h file may be included in host MD codes, **NEVER** MODIFY THE CODE IN THIS FILE
// unless you know exactly what you are doing

/**
  Container for plumedmain function pointers (create, cmd and finalize).
*/
typedef struct {
  void*(*create)();
  void(*cmd)(void*,const char*,const void*);
  void(*finalize)(void*);
} plumed_plumedmain_function_holder;

/**
  Container for symbol table. Presently only contains a version number and a plumed_plumedmain_function_holder object.
*/
typedef struct {
  int version;
  plumed_plumedmain_function_holder functions;
} plumed_symbol_table_type;

/* These functions should be accessible from C, since they might be statically
   used from Plumed.c (for static binding) */

extern "C" void*plumed_plumedmain_create() {
  return new PLMD::PlumedMain;
}

extern "C" void plumed_plumedmain_cmd(void*plumed,const char*key,const void*val) {
  plumed_massert(plumed,"trying to use a plumed object which is not initialized");
  static_cast<PLMD::PlumedMain*>(plumed)->cmd(key,val);
}

extern "C" void plumed_plumedmain_finalize(void*plumed) {
  plumed_massert(plumed,"trying to deallocate a plumed object which is not initialized");
// I think it is not possible to replace this delete with a smart pointer
// since the ownership of this pointer is in a C structure. GB
  delete static_cast<PLMD::PlumedMain*>(plumed);
}

/// This is a static structure that is searched by the plumed loader, so it should be here.
extern "C" plumed_symbol_table_type plumed_symbol_table;


/// Notice that this object is only searched with dlsym after a dlopen on libplumedKernel
/// has completed. Thus, it should not suffer for any static initialization problem.
/// It is initialized using the plumed_symbol_table_init
plumed_symbol_table_type plumed_symbol_table;

extern "C" void plumed_symbol_table_init() {
  plumed_symbol_table.version=1;
  plumed_symbol_table.functions.create=plumed_plumedmain_create;
  plumed_symbol_table.functions.cmd=plumed_plumedmain_cmd;
  plumed_symbol_table.functions.finalize=plumed_plumedmain_finalize;
}

namespace PLMD {

/// Static object which registers Plumed.
/// This is a static object which, during its construction at startup,
/// registers the pointers to plumed_plumedmain_create, plumed_plumedmain_cmd and plumed_plumedmain_finalize
/// to the plumed_kernel_register function.
/// Registration is only required with plumed loader <=2.4, but we do it anyway in order to maintain
/// backward compatibility. Notice that as of plumed 2.5 the plumed_kernel_register is found
/// using dlsym, in order to allow the libplumedKernel library to be loadable also when
/// the plumed_kernel_register symbol is not available.
static class PlumedMainInitializer {
public:
  PlumedMainInitializer() {
#if defined(__PLUMED_HAS_DLOPEN)
    bool debug=std::getenv("PLUMED_LOAD_DEBUG");
    if(std::getenv("PLUMED_LOAD_SKIP_REGISTRATION")) {
      if(debug) fprintf(stderr,"+++ Skipping registration +++\n");
      return;
    }
    typedef plumed_plumedmain_function_holder* (*plumed_kernel_register_type)(const plumed_plumedmain_function_holder*);
    plumed_kernel_register_type plumed_kernel_register=nullptr;
    void* handle=nullptr;
#if defined(__PLUMED_HAS_RTLD_DEFAULT)
    if(debug) fprintf(stderr,"+++ Registering functions. Looking in RTLD_DEFAULT +++\n");
    plumed_kernel_register=(plumed_kernel_register_type) dlsym(RTLD_DEFAULT,"plumed_kernel_register");
#else
    handle=dlopen(NULL,RTLD_LOCAL);
    if(debug) fprintf(stderr,"+++ Registering functions. dlopen handle at %p +++\n",handle);
    plumed_kernel_register=(plumed_kernel_register_type) dlsym(handle,"plumed_kernel_register");
#endif
    if(debug) {
      if(plumed_kernel_register) {
        fprintf(stderr,"+++ plumed_kernel_register found at %p +++\n",(void*)plumed_kernel_register);
      }
      else fprintf(stderr,"+++ plumed_kernel_register not found +++\n");
    }
// make sure static plumed_function_pointers is initialized here
    plumed_symbol_table_init();
    if(plumed_kernel_register && debug) fprintf(stderr,"+++ Registering functions at %p (%p,%p,%p) +++\n",
          (void*)&plumed_symbol_table.functions,(void*)plumed_symbol_table.functions.create,(void*)plumed_symbol_table.functions.cmd,(void*)plumed_symbol_table.functions.finalize);
    if(plumed_kernel_register) (*plumed_kernel_register)(&plumed_symbol_table.functions);
// Notice that handle could be null in the following cases:
// - if we use RTLD_DEFAULT
// - on Linux if we don't use RTLD_DEFAULT, since dlopen(NULL,RTLD_LOCAL) returns a null pointer.
    if(handle) dlclose(handle);
#endif
  }
} RegisterMe;

}


