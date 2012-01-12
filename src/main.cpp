#include "CLTool.h"
#ifdef __PLUMED_MPI
#include "mpi.h"
#endif

using namespace PLMD;

/// This main just wraps the CLTool::globalMain static
/// function. The object file generated from this .cpp
/// is the only part of the plumed library that should
/// not be linked with external MD codes, so as 
/// to avoid linker error.
int main(int argc,char**argv){
#ifdef __PLUMED_MPI
  MPI_Init(&argc,&argv);
#endif
  int ret=CLTool::globalMain(argc,argv);
#ifdef __PLUMED_MPI
  MPI_Finalize();
#endif
  return ret;
}
