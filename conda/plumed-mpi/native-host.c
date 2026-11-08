/* Packaging probe for older hosts: pass a native communicator without injection.
   The packaged provider must be discovered beside the installed kernel. */
#include <mpi.h>
#include <plumed/wrapper/Plumed.h>

int main(int argc, char **argv) {
  int result=1;
  MPI_Comm comm;
  MPI_Init(&argc,&argv);
  comm=MPI_COMM_WORLD;
  plumed p=plumed_create();
  plumed_cmd(p,"CLTool setArgc",&argc);
  plumed_cmd(p,"CLTool setArgv",argv);
  plumed_cmd(p,"CLTool setMPIComm",&comm);
  plumed_cmd(p,"CLTool run",&result);
  plumed_finalize(p);
  MPI_Finalize();
  return result;
}
