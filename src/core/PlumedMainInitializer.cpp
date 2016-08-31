/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
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

/* These functions should be accessible from C, since they might be statically
   used from Plumed.c (for static binding) */

extern "C" void*plumedmain_create(){
  return new PLMD::PlumedMain;
}

extern "C" void plumedmain_cmd(void*plumed,const char*key,const void*val){
  plumed_massert(plumed,"trying to use a plumed object which is not initialized");
  static_cast<PLMD::PlumedMain*>(plumed)->cmd(key,val);
}

extern "C" void plumedmain_finalize(void*plumed){
  plumed_massert(plumed,"trying to deallocate a plumed object which is not initialized");
  delete static_cast<PLMD::PlumedMain*>(plumed);
}

/* This refers to a function implemented in Plumed.c */
extern "C" plumed_plumedmain_function_holder* plumed_kernel_register(const plumed_plumedmain_function_holder*);

namespace PLMD{

/// Static object which registers Plumed.
/// This is a static object which, during its construction at startup,
/// registers the pointers to plumedmain_create, plumedmain_cmd and plumedmain_finalize
/// to the plumed_kernel_register function
static class PlumedMainInitializer{
  public:
  PlumedMainInitializer(){
    plumed_plumedmain_function_holder fh={plumedmain_create,plumedmain_cmd,plumedmain_finalize};
    plumed_kernel_register(&fh);
  }
} RegisterMe;

}


