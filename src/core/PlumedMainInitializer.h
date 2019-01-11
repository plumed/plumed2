/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2018,2019 The plumed team
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
#ifndef __PLUMED_core_PlumedMainInitializer_h
#define __PLUMED_core_PlumedMainInitializer_h
// !!!!!!!!!!!!!!!!!!!!!!    DANGER   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
// THE FOLLOWING ARE DEFINITIONS WHICH ARE NECESSARY FOR DYNAMIC LOADING OF THE PLUMED KERNEL:
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

typedef struct {
  void* ptr;
  void (*handler)(void*,int,const char*,const void*);
} plumed_nothrow_handler;

/**
  Container for symbol table. Presently only contains a version number and a plumed_plumedmain_function_holder object.
  The definition of this structure might change in the future. In particular, the structure might grow adding
  new fields. However, in that case the version number should be updated as well.
*/
typedef struct {
  int version;
  plumed_plumedmain_function_holder functions;
  void (*cmd_nothrow)(void*plumed,const char*key,const void*val,plumed_nothrow_handler nothrow);
} plumed_symbol_table_type;


// additional definitions
typedef void*(*plumed_create_pointer)(void);
typedef void(*plumed_cmd_pointer)(void*,const char*,const void*);
typedef void(*plumed_finalize_pointer)(void*);

/* These functions should be accessible from C, since they might be statically
   used from Plumed.c (for static binding) */

/**
  Constructs a plumed object.
  This function returns a void pointer that can be used in \ref plumed_plumedmain_cmd and \ref plumed_plumedmain_finalize.
*/
extern "C" void*plumed_plumedmain_create();

/**
  Send a command `key` and a pointer `val` to a void pointer returned by \ref plumed_plumedmain_create.
*/
extern "C" void plumed_plumedmain_cmd(void*plumed,const char*key,const void*val);

/**
  Finalize a void pointer returned by \ref plumed_plumedmain_create
*/
extern "C" void plumed_plumedmain_finalize(void*plumed);

/**
  Static symbol table that is accessed by the plumed loader.
  Notice that this table is initialized with a static object construction.
  In principle, it should be accessed by other programs dlopening the plumed kernel.
  In that case, it is guaranteed to be already initialized.
  However, when accessed directly it might be safer to first call \ref plumed_symbol_table_init.
*/

extern "C" plumed_symbol_table_type plumed_symbol_table;

/**
  Function that makes sure that \ref plumed_symbol_table is initialized.
  Can be called multiple times.
*/
extern "C" void plumed_symbol_table_init();

namespace PLMD {
// This is just to avoid plumedcheck warnings.
// Notice that the only define C-style objects here, so namespace is not needed
}


#endif
