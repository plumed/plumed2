#include "plumed/mpiwrap/mpiwrap.h"
#include "plumed/mpiwrap/mpi.h"

#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>
#include "plumed/tools/Communicator.h"

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
  Communicator comm;
  comm.Set_comm(world);
  // Exercise typed reductions as well as copying: a datatype alias can hide
  // buffer/type mismatches until the provider interprets an arithmetic value.
  const unsigned long base=std::numeric_limits<unsigned long>::digits>32
                           ? static_cast<unsigned long>(1ULL<<40) : 1000000UL;
  unsigned long sum=base+rank;
  comm.Sum(sum);
  if(sum!=4*base+6) {
    return 1;
  }
  unsigned long broadcast=rank==3 ? base+73 : 0;
  comm.Bcast(broadcast,3);
  if(broadcast!=base+73) {
    return 1;
  }
  // Check the two-word transport used for full-width checkpoint/OPES counts,
  // including a high word and the carry across the low-word boundary.
  const unsigned long long count=(1ULL<<48)+0xffffffffULL+rank;
  const std::vector<unsigned long> words= {static_cast<unsigned long>(count&0xffffffffULL),
                                           static_cast<unsigned long>(count>>32)
                                          };
  std::vector<unsigned long> gathered(8);
  comm.Allgather(words,gathered);
  for(unsigned r=0; r<4; ++r) {
    const auto actual=static_cast<unsigned long long>(gathered[2*r])
                      | (static_cast<unsigned long long>(gathered[2*r+1])<<32);
    if(actual!=(1ULL<<48)+0xffffffffULL+r) {
      return 1;
    }
  }
  if(rank==0) {
    FILE *output=std::fopen("output","w");
    if(!output) {
      return 1;
    }
    std::fputs("unsigned-long collectives and full-width count transport passed\n",output);
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
#if defined(__PLUMED_HAS_THREAD_MPI) && !defined(__PLUMED_HAS_MPI)
  const auto lifecycle=mpiwrap::MPILifecyclePolicy::internal_thread_mpi;
#else
  const auto lifecycle=mpiwrap::MPILifecyclePolicy::standalone_owned;
#endif
  const int status=mpiwrap::launchRanks(4,body,nullptr,&result,lifecycle);
  return status==0 ? result : 1;
}

}
}

int main(int argc,char **argv) {
  return PLMD::launcher_regtest::run(argc,argv);
}
