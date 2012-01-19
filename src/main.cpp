#include "CLTool.h"
#include <cstring>
#ifdef __PLUMED_MPI
#include "mpi.h"
#endif

using namespace PLMD;
using namespace std;

/// This main just wraps the CLTool::globalMain static
/// function. The object file generated from this .cpp
/// is the only part of the plumed library that should
/// not be linked with external MD codes, so as 
/// to avoid linker error.
int main(int argc,char**argv){
  bool nompi=false;
  if(argc>1 && !strcmp(argv[1],"--no-mpi")) nompi=true;
#ifdef __PLUMED_MPI
  if(!nompi) MPI_Init(&argc,&argv);
#endif
  int ret=CLTool::globalMain(argc,argv);
#ifdef __PLUMED_MPI
  if(!nompi) MPI_Finalize();
#endif
  return ret;
}
