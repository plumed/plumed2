/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
// PLUMED_C_COMPATIBLE
#ifdef __PLUMED_HAS_MPI
#define __PLUMED_WRAPPER_MPI_HOST_ADAPTER 1
#include <mpi.h>
#endif
#define __PLUMED_WRAPPER_CXX 0
#include "wrapper/Plumed.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Use the C error interface so this same source can be compiled with mpicc.
// Keep each command's failure path identical, including object cleanup.
#define PLUMED_MAIN_CMD(key,value) do { \
    plumed_cmd_nothrow(p,key,value,nothrow); \
    if(error.code) { \
      goto failure; \
    } \
  } while(0)

/**
  This main uses only the interface published in
  Plumed.h. The object file generated from this .cpp
  is the only part of the plumed library that should
  not be linked with external MD codes, so as
  to avoid linker error.
*/
int main(int argc,char**argv) {
  int ret=0;
  plumed p;
  plumed_error error;
  plumed_nothrow_handler nothrow= {&error,plumed_error_set};
#ifdef __PLUMED_HAS_MPI
  int nompi=0;
  for(int iarg=1; iarg<argc; iarg++) {
    if(!strcmp(argv[iarg],"--no-mpi")) {
      nompi=1;
    }
    if(!strcmp(argv[iarg],"--mpi")) {
      nompi=0;
    }
// stop at first non-option
    if(argv[iarg] && argv[iarg][0]!='-') {
      break;
    }
  }
  if(!nompi) {
    MPI_Init(&argc,&argv);
  }
#endif
  p=plumed_create();
  plumed_error_init(&error);
#ifdef __PLUMED_HAS_MPI
  {
    // Older kernels predate injection; they retain the native-MPI path.
    int api=0;
    PLUMED_MAIN_CMD("getApiVersion",&api);
    if(api>=11) {
      PLUMED_MAIN_CMD("injectMPIImplementation",plumed_mpi_implementation());
    }
  }
#endif
  PLUMED_MAIN_CMD("CLTool setArgc",&argc);
  PLUMED_MAIN_CMD("CLTool setArgv",argv);
#ifdef __PLUMED_HAS_MPI
  if(!nompi) {
    MPI_Comm comm;
    MPI_Comm_dup(MPI_COMM_WORLD,&comm);
    PLUMED_MAIN_CMD("CLTool setMPIComm",&comm);
  }
#endif
  PLUMED_MAIN_CMD("CLTool run",&ret);
  plumed_finalize(p);

#ifdef __PLUMED_HAS_MPI
  if(!nompi) {
    MPI_Finalize();
  }
#endif
  return ret;

failure:
  plumed_error_merge_with_nested(&error);
  fprintf(stderr,"%s\n",plumed_error_what(error));
  plumed_error_finalize(error);
  plumed_finalize(p);
  // As with an uncaught exception, terminate instead of entering collective
  // finalization while another MPI rank might still be blocked.
  abort();
}

#undef PLUMED_MAIN_CMD
