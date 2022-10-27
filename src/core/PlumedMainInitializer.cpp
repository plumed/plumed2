/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2022 The plumed team
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
#include "lepton/Exception.h"
#include <cstdlib>
#include <cstring>
#include <iostream>
#if defined __PLUMED_HAS_DLOPEN
#include <dlfcn.h>
#endif
#include <exception>
#include <stdexcept>
#include <ios>
#include <new>
#include <typeinfo>
#ifdef __PLUMED_LIBCXX11
#include <system_error>
#include <future>
#include <memory>
#include <functional>
#endif
#include "tools/TypesafePtr.h"


static bool getenvTypesafeDebug() noexcept {
  static const auto* res=std::getenv("PLUMED_TYPESAFE_DEBUG");
  return res;
}

static void typesafeDebug(const char*key,plumed_safeptr_x safe) noexcept {
  std::fprintf(stderr,"+++ PLUMED_TYPESAFE_DEBUG %s %p %zu",key,safe.ptr,safe.nelem);
  const size_t* shape=safe.shape;
  if(shape) {
    std::fprintf(stderr," (");
    while(*shape!=0) {
      std::fprintf(stderr," %zu",*shape);
      shape++;
    }
    std::fprintf(stderr," )");
  }
  std::fprintf(stderr," %zx %p\n",safe.flags,safe.opt);
}

// create should never throw
// in case of a problem, it logs the error and return a null pointer
// when loaded by an interface >=2.5, this will result in a non valid plumed object.
// earlier interfaces will just give a segfault or a failed assertion.
extern "C" void*plumed_plumedmain_create() {
  try {
    return new PLMD::PlumedMain;
  } catch(const std::exception & e) {
    std::cerr<<"+++ an error happened while creating a plumed object\n";
    std::cerr<<e.what()<<std::endl;
    return nullptr;
  } catch(...) {
    std::cerr<<"+++ an unknown error happened while creating a plumed object"<<std::endl;
    return nullptr;
  }
}

extern "C" unsigned plumed_plumedmain_create_reference(void*plumed) {
  plumed_massert(plumed,"trying to create a reference to a plumed object which is not initialized");
  auto p=static_cast<PLMD::PlumedMain*>(plumed);
  return p->increaseReferenceCounter();
}

extern "C" unsigned plumed_plumedmain_delete_reference(void*plumed) {
  plumed_massert(plumed,"trying to delete a reference to a plumed object which is not initialized");
  auto p=static_cast<PLMD::PlumedMain*>(plumed);
  return p->decreaseReferenceCounter();
}

extern "C" unsigned plumed_plumedmain_use_count(void*plumed) {
  plumed_massert(plumed,"trying to delete a reference to a plumed object which is not initialized");
  auto p=static_cast<PLMD::PlumedMain*>(plumed);
  return p->useCountReferenceCounter();
}

extern "C" void plumed_plumedmain_cmd(void*plumed,const char*key,const void*val) {
  plumed_massert(plumed,"trying to use a plumed object which is not initialized");
  auto p=static_cast<PLMD::PlumedMain*>(plumed);
  p->cmd(key,PLMD::TypesafePtr::unchecked(val));
}

extern "C" {
  static void plumed_plumedmain_cmd_safe(void*plumed,const char*key,plumed_safeptr_x safe) {
    plumed_massert(plumed,"trying to use a plumed object which is not initialized");
    auto p=static_cast<PLMD::PlumedMain*>(plumed);
    if(getenvTypesafeDebug()) typesafeDebug(key,safe);
    p->cmd(key,PLMD::TypesafePtr::fromSafePtr(&safe));
  }
}

extern "C" {
  static void plumed_plumedmain_cmd_safe_nothrow(void*plumed,const char*key,plumed_safeptr_x safe,plumed_nothrow_handler_x nothrow) {
// This is a workaround for a suboptimal choice in PLUMED <2.8
// In particular, the only way to bypass the exception handling process was to call the plumed_plumedmain_cmd_safe
// function directly.
// With this modification, it is possible to just call the plumed_plumedmain_cmd_safe_nothrow function
// passing a null error handler.
    if(!nothrow.handler) {
      plumed_plumedmain_cmd_safe(plumed,key,safe);
      return;
    }
// At library boundaries we translate exceptions to error codes.
// This allows an exception to be catched also if the MD code
// was linked against a different C++ library
    try {
      plumed_massert(plumed,"trying to use a plumed object which is not initialized");
      auto p=static_cast<PLMD::PlumedMain*>(plumed);
      if(getenvTypesafeDebug()) typesafeDebug(key,safe);
      p->cmd(key,PLMD::TypesafePtr::fromSafePtr(&safe));
    } catch(const PLMD::ExceptionTypeError & e) {
      nothrow.handler(nothrow.ptr,20300,e.what(),nullptr);
    } catch(const PLMD::ExceptionError & e) {
      nothrow.handler(nothrow.ptr,20200,e.what(),nullptr);
    } catch(const PLMD::ExceptionDebug & e) {
      nothrow.handler(nothrow.ptr,20100,e.what(),nullptr);
    } catch(const PLMD::Exception & e) {
      nothrow.handler(nothrow.ptr,20000,e.what(),nullptr);
    } catch(const PLMD::lepton::Exception & e) {
      nothrow.handler(nothrow.ptr,19900,e.what(),nullptr);
      // 11000 to 12000 are "bad exceptions". message will be copied without new allocations
    } catch(const std::bad_exception & e) {
      nothrow.handler(nothrow.ptr,11500,e.what(),nullptr);
#ifdef __PLUMED_LIBCXX11
    } catch(const std::bad_array_new_length & e) {
      nothrow.handler(nothrow.ptr,11410,e.what(),nullptr);
#endif
    } catch(const std::bad_alloc & e) {
      nothrow.handler(nothrow.ptr,11400,e.what(),nullptr);
#ifdef __PLUMED_LIBCXX11
    } catch(const std::bad_function_call & e) {
      nothrow.handler(nothrow.ptr,11300,e.what(),nullptr);
    } catch(const std::bad_weak_ptr & e) {
      nothrow.handler(nothrow.ptr,11200,e.what(),nullptr);
#endif
    } catch(const std::bad_cast & e) {
      nothrow.handler(nothrow.ptr,11100,e.what(),nullptr);
    } catch(const std::bad_typeid & e) {
      nothrow.handler(nothrow.ptr,11000,e.what(),nullptr);
      // not implemented yet: std::regex_error
      // we do not allow regex yet due to portability problems with gcc 4.8
      // as soon as we transition to using <regex> it should be straightforward to add
    } catch(const std::ios_base::failure & e) {
#ifdef __PLUMED_LIBCXX11
      int value=e.code().value();
      const void* opt[3]= {"c",&value,nullptr}; // "c" passes the error code. nullptr terminates the optional part.
      if(e.code().category()==std::generic_category()) nothrow.handler(nothrow.ptr,10230,e.what(),opt);
      else if(e.code().category()==std::system_category()) nothrow.handler(nothrow.ptr,10231,e.what(),opt);
      else if(e.code().category()==std::iostream_category()) nothrow.handler(nothrow.ptr,10232,e.what(),opt);
      else if(e.code().category()==std::future_category()) nothrow.handler(nothrow.ptr,10233,e.what(),opt);
      else
#endif
        // 10239 represents std::ios_base::failure with default constructur
        nothrow.handler(nothrow.ptr,10239,e.what(),nullptr);
#ifdef __PLUMED_LIBCXX11
    } catch(const std::system_error & e) {
      int value=e.code().value();
      const void* opt[3]= {"c",&value,nullptr}; // "c" passes the error code. nullptr terminates the optional part.
      if(e.code().category()==std::generic_category()) nothrow.handler(nothrow.ptr,10220,e.what(),opt);
      else if(e.code().category()==std::system_category()) nothrow.handler(nothrow.ptr,10221,e.what(),opt);
      else if(e.code().category()==std::iostream_category()) nothrow.handler(nothrow.ptr,10222,e.what(),opt);
      else if(e.code().category()==std::future_category()) nothrow.handler(nothrow.ptr,10223,e.what(),opt);
      // fallback to generic runtime_error
      else nothrow.handler(nothrow.ptr,10200,e.what(),nullptr);
#endif
    } catch(const std::underflow_error &e) {
      nothrow.handler(nothrow.ptr,10215,e.what(),nullptr);
    } catch(const std::overflow_error &e) {
      nothrow.handler(nothrow.ptr,10210,e.what(),nullptr);
    } catch(const std::range_error &e) {
      nothrow.handler(nothrow.ptr,10205,e.what(),nullptr);
    } catch(const std::runtime_error & e) {
      nothrow.handler(nothrow.ptr,10200,e.what(),nullptr);
      // not implemented yet: std::future_error
      // not clear how useful it would be.
    } catch(const std::out_of_range & e) {
      nothrow.handler(nothrow.ptr,10120,e.what(),nullptr);
    } catch(const std::length_error & e) {
      nothrow.handler(nothrow.ptr,10115,e.what(),nullptr);
    } catch(const std::domain_error & e) {
      nothrow.handler(nothrow.ptr,10110,e.what(),nullptr);
    } catch(const std::invalid_argument & e) {
      nothrow.handler(nothrow.ptr,10105,e.what(),nullptr);
    } catch(const std::logic_error & e) {
      nothrow.handler(nothrow.ptr,10100,e.what(),nullptr);
      // generic exception. message will be copied without new allocations
      // reports all non caught exceptions that are derived from std::exception
      // for instance, boost exceptions would end up here
    } catch(const std::exception & e) {
      nothrow.handler(nothrow.ptr,10000,e.what(),nullptr);
    } catch(...) {
      // if exception cannot be translated, we throw a bad_exception
      nothrow.handler(nothrow.ptr,11500,"plumed could not translate exception",nullptr);
      throw;
    }
  }
}

extern "C" {
  static void plumed_plumedmain_cmd_nothrow(void*plumed,const char*key,const void*val,plumed_nothrow_handler_x nothrow) {
    plumed_safeptr_x safe;
    plumed_assert(nothrow.handler) << "Accepting a null pointer here would make the calling code non compatible with plumed 2.5 to 2.7";
    safe.ptr=val;
    safe.nelem=0;
    safe.shape=NULL;
    safe.flags=0;
    safe.opt=NULL;
    plumed_plumedmain_cmd_safe_nothrow(plumed,key,safe,nothrow);
  }
}

extern "C" void plumed_plumedmain_finalize(void*plumed) {
  plumed_massert(plumed,"trying to deallocate a plumed object which is not initialized");
// I think it is not possible to replace this delete with a smart pointer
// since the ownership of this pointer is in a C structure. GB
  delete static_cast<PLMD::PlumedMain*>(plumed);
}

// values here should be consistent with those in plumed_symbol_table_init !!!!
plumed_symbol_table_type_x plumed_symbol_table= {
  4,
  {plumed_plumedmain_create,plumed_plumedmain_cmd,plumed_plumedmain_finalize},
  plumed_plumedmain_cmd_nothrow,
  plumed_plumedmain_cmd_safe,
  plumed_plumedmain_cmd_safe_nothrow,
  plumed_plumedmain_create_reference,
  plumed_plumedmain_delete_reference,
  plumed_plumedmain_use_count
};

// values here should be consistent with those above !!!!
extern "C" void plumed_symbol_table_init() {
  plumed_symbol_table.version=4;
  plumed_symbol_table.functions.create=plumed_plumedmain_create;
  plumed_symbol_table.functions.cmd=plumed_plumedmain_cmd;
  plumed_symbol_table.functions.finalize=plumed_plumedmain_finalize;
  plumed_symbol_table.cmd_nothrow=plumed_plumedmain_cmd_nothrow;
  plumed_symbol_table.cmd_safe=plumed_plumedmain_cmd_safe;
  plumed_symbol_table.cmd_safe_nothrow=plumed_plumedmain_cmd_safe_nothrow;
  plumed_symbol_table.create_reference=plumed_plumedmain_create_reference;
  plumed_symbol_table.delete_reference=plumed_plumedmain_delete_reference;
  plumed_symbol_table.use_count=plumed_plumedmain_use_count;
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
namespace {
class PlumedMainInitializer {
  const bool debug;
public:
  PlumedMainInitializer():
    debug(std::getenv("PLUMED_LOAD_DEBUG"))
  {
// make sure static plumed_function_pointers is initialized here
    plumed_symbol_table_init();
    if(debug) std::fprintf(stderr,"+++ Initializing PLUMED with plumed_symbol_table version %i at %p\n",plumed_symbol_table.version,(void*)&plumed_symbol_table);
#if defined(__PLUMED_HAS_DLOPEN)
    if(std::getenv("PLUMED_LOAD_SKIP_REGISTRATION")) {
      if(debug) std::fprintf(stderr,"+++ Skipping registration +++\n");
      return;
    }
    typedef plumed_plumedmain_function_holder_x* (*plumed_kernel_register_type_x)(const plumed_plumedmain_function_holder_x*);
    plumed_kernel_register_type_x plumed_kernel_register=nullptr;
    void* handle=nullptr;
#if defined(__PLUMED_HAS_RTLD_DEFAULT)
    if(debug) std::fprintf(stderr,"+++ Registering functions. Looking in RTLD_DEFAULT +++\n");
    void* dls=dlsym(RTLD_DEFAULT,"plumed_kernel_register");
#else
    handle=dlopen(NULL,RTLD_LOCAL);
    if(debug) std::fprintf(stderr,"+++ Registering functions. dlopen handle at %p +++\n",handle);
    void* dls=dlsym(handle,"plumed_kernel_register");
#endif
    *(void **)(&plumed_kernel_register)=dls;
    if(debug) {
      if(plumed_kernel_register) {
        std::fprintf(stderr,"+++ plumed_kernel_register found at %p +++\n",dls);
      }
      else std::fprintf(stderr,"+++ plumed_kernel_register not found +++\n");
    }
    void*createp;
    void*cmdp;
    void*finalizep;
    plumed_convert_fptr(createp,plumed_symbol_table.functions.create);
    plumed_convert_fptr(cmdp,plumed_symbol_table.functions.cmd);
    plumed_convert_fptr(finalizep,plumed_symbol_table.functions.finalize);
    if(plumed_kernel_register && debug) std::fprintf(stderr,"+++ Registering functions at %p (%p,%p,%p) +++\n",
          (void*)&plumed_symbol_table.functions,createp,cmdp,finalizep);
    if(plumed_kernel_register) (*plumed_kernel_register)(&plumed_symbol_table.functions);
// Notice that handle could be null in the following cases:
// - if we use RTLD_DEFAULT
// - on Linux if we don't use RTLD_DEFAULT, since dlopen(NULL,RTLD_LOCAL) returns a null pointer.
    if(handle) dlclose(handle);
#endif
  }
  ~PlumedMainInitializer() {
    if(debug) std::fprintf(stderr,"+++ Finalizing PLUMED with plumed_symbol_table at %p\n",(void*)&plumed_symbol_table);
  }
} PlumedMainInitializerRegisterMe;
}

}


