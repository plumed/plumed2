/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
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
#include "wrapper/Plumed.h"
#include <cstring>

#ifdef __PLUMED_HAS_MPI
#include <mpi.h>
#endif

/**
  This main uses only the interface published in
  Plumed.h. The object file generated from this .cpp
  is the only part of the plumed library that should
  not be linked with external MD codes, so as
  to avoid linker error.
*/
int main(int argc,char**argv) {
#ifdef __PLUMED_HAS_MPI
  bool nompi=false;
  for(unsigned iarg=1; iarg<argc; iarg++) {
    if(!strcmp(argv[iarg],"--no-mpi")) nompi=true;
    if(!strcmp(argv[iarg],"--mpi"))    nompi=false;
// stop at first non-option
    if(argv[iarg] && argv[iarg][0]!='-') break;
  }
  if(!nompi) MPI_Init(&argc,&argv);
#endif
  int ret=0;

  try {
    PLMD::Plumed p;
    p.cmd("CLTool setArgc",&argc);
    p.cmd("CLTool setArgv",argv);
#ifdef __PLUMED_HAS_MPI
    if(!nompi) {
      MPI_Comm comm;
      MPI_Comm_dup(MPI_COMM_WORLD,&comm);
      p.cmd("CLTool setMPIComm",&comm);
    }
#endif
    p.cmd("CLTool run",&ret);
// end block deletes p also in case an exception occurs
  } catch(...) {
// exception is rethrown and results in a call to terminate
    throw;
  }

#ifdef __PLUMED_HAS_MPI
  if(!nompi) MPI_Finalize();
#endif
  return ret;
}
