#include "Plumed.h"
#include <cstring>

#ifdef __PLUMED_MPI
#include <mpi.h>
#endif

using namespace std;

/**
  This main uses only the interface published in
  Plumed.h. The object file generated from this .cpp
  is the only part of the plumed library that should
  not be linked with external MD codes, so as 
  to avoid linker error.
*/
int main(int argc,char**argv){
  bool nompi=false;
  if(argc>1 && !strcmp(argv[1],"--no-mpi")) nompi=true;
#ifdef __PLUMED_MPI
  if(!nompi) MPI_Init(&argc,&argv);
#endif
  int ret;

  PLMD::Plumed* p=new PLMD::Plumed;
  p->cmd("CLTool setArgc",&argc);
  p->cmd("CLTool setArgv",argv);
#ifdef __PLUMED_MPI
  if(!nompi){
    MPI_Comm comm;
    MPI_Comm_dup(MPI_COMM_WORLD,&comm);
    p->cmd("CLTool setMPIComm",&comm);
  }
#endif
  p->cmd("CLTool run",&ret);
  delete p;

#ifdef __PLUMED_MPI
  if(!nompi) MPI_Finalize();
#endif
  return ret;
}
