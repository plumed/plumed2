#ifndef PLUMED_REGTEST_MPI_TEST_H
#define PLUMED_REGTEST_MPI_TEST_H

#if defined(__PLUMED_MPIWRAP)
#include "plumed/mpiwrap/mpiwrap.h"
#include "plumed/mpiwrap/thread-mpi-launcher.h"
#else
// Keep these tests usable with older native-MPI PLUMED installations.
#include <mpi.h>
#endif
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace PLMD {
namespace mpi_test {

using RankBody = int (*)(MPI_Comm, int, char **);

#if defined(__PLUMED_MPIWRAP)
struct Context {
  RankBody body;
  int argc;
  char **argv;
};

inline int enterRank(int, const void *native_world, std::size_t, void *opaque) {
  const auto &context=*static_cast<const Context *>(opaque);
  return context.body(mpiwrap_Comm_from_native(native_world),context.argc,context.argv);
}
#endif

// All objects and mutable buffers belong inside the rank body. The argument
// strings are shared read-only by internal ranks. Older runners supply no -n;
// their existing mpiprocs setting still selects the native process count.
inline int run(int argc, char **argv, int ranks, RankBody body) {
  if(argc>1 && std::strcmp(argv[1],"-n")==0) {
    if(argc<3) {
      return 2;
    }
    char *end=nullptr;
    const long requested=std::strtol(argv[2],&end,10);
    if(!end || *end || requested<1 || requested>INT_MAX) {
      return 2;
    }
    ranks=static_cast<int>(requested);
    for(int i=1; i+2<=argc; ++i) {
      argv[i]=argv[i+2];
    }
    argc-=2;
  }
#if defined(__PLUMED_MPIWRAP)
  Context context {body,argc,argv};
  int result=1;
#if defined(__PLUMED_HAS_MPI)
  const auto lifecycle=mpiwrap_provider::MPILifecyclePolicy::standalone_owned;
#else
  const auto lifecycle=mpiwrap_provider::MPILifecyclePolicy::internal_thread_mpi;
#endif
  const int status=mpiwrap_provider::launchRanks(ranks,enterRank,&context,&result,lifecycle);
  return status==0 ? result : 1;
#else
  if(MPI_Init(&argc,&argv)!=MPI_SUCCESS) {
    return 1;
  }
  int result=1;
  try {
    result=body(MPI_COMM_WORLD,argc,argv);
  } catch(...) {
    MPI_Abort(MPI_COMM_WORLD,1);
    return 1;
  }
  return MPI_Finalize()==MPI_SUCCESS ? result : 1;
#endif
}

} // namespace mpi_test
} // namespace PLMD
#endif
