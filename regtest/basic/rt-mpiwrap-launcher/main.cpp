#include "plumed/mpiwrap/mpiwrap.h"
#include "plumed/mpiwrap/mpi.h"

#include <cstdio>
#include <cstring>

namespace PLMD {
namespace launcher_regtest {

int body(int rank, const void *native_world, std::size_t native_size, void *) {
  if(!native_world || native_size==0) {
    return 1;
  }
  // Every callback owns its communicator and communication buffers.
  const MPI_Comm world=mpiwrap::Comm_from_native(native_world);
  int actual_rank=-1, size=0;
  if(MPI_Comm_rank(world,&actual_rank)!=MPI_SUCCESS || actual_rank!=rank ||
      MPI_Comm_size(world,&size)!=MPI_SUCCESS || size!=4) {
    return 1;
  }
  const int input=rank+1;
  int sum=0;
  if(MPI_Allreduce(&input,&sum,1,MPI_INT,MPI_SUM,world)!=MPI_SUCCESS || sum!=10) {
    return 1;
  }
  int broadcast=rank==size-1 ? 73 : 0;
  if(MPI_Bcast(&broadcast,1,MPI_INT,size-1,world)!=MPI_SUCCESS || broadcast!=73) {
    return 1;
  }
  if(rank==0) {
    FILE *output=std::fopen("output","w");
    if(!output) {
      return 1;
    }
    std::fputs("four-rank standalone launcher passed\n",output);
    if(std::fclose(output)!=0) {
      return 1;
    }
  }
  return 0;
}

int run(int argc, char **argv) {
  if(argc!=3 || std::strcmp(argv[1],"-n")!=0 || std::strcmp(argv[2],"4")!=0) {
    return 2;
  }
  int result=-1;
#if defined(__PLUMED_HAS_MPI)
  const auto lifecycle=mpiwrap::MPILifecyclePolicy::standalone_owned;
#else
  const auto lifecycle=mpiwrap::MPILifecyclePolicy::internal_thread_mpi;
#endif
  const int status=mpiwrap::launchRanks(4,body,nullptr,&result,lifecycle);
  return status==0 ? result : 1;
}

}
}

int main(int argc,char **argv) {
  return PLMD::launcher_regtest::run(argc,argv);
}
