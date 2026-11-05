#if defined(__PLUMED_MPIWRAP) && defined(__PLUMED_HAS_MPI)
#include "thread-mpi-launcher.h"
#include "../tools/Exception.h"

#include <mpi.h>

#include <cstdio>
#include <cstdlib>
#include <exception>

namespace PLMD {
namespace mpiwrap_provider {

int launchRanks(int ranks,
                int (*callback)(int, const void *, std::size_t, void *),
                void *callback_context, int *result,
                MPILifecyclePolicy lifecycle) {
  if (ranks < 1 || !callback || !result) {
    throw ExceptionTypeError("rank count, callback, and result output are required");
  }
  if (lifecycle == MPILifecyclePolicy::internal_thread_mpi) {
    throw ExceptionError("internal thread-MPI launch is unavailable in a native-MPI kernel");
  }
  if (lifecycle != MPILifecyclePolicy::host_owned &&
      lifecycle != MPILifecyclePolicy::standalone_owned) {
    throw ExceptionTypeError("unknown MPI lifecycle policy");
  }
  int initialized = 0;
  int status = MPI_Initialized(&initialized);
  if (status != MPI_SUCCESS) {
    throw ExceptionError("MPI_Initialized failed before launch");
  }
  int finalized = 0;
  status = MPI_Finalized(&finalized);
  if (status != MPI_SUCCESS || finalized) {
    throw ExceptionError("MPI is finalized or MPI_Finalized failed before launch");
  }
  const bool owns_mpi = lifecycle == MPILifecyclePolicy::standalone_owned;
  if (owns_mpi && initialized) {
    throw ExceptionError("standalone MPI launcher requires MPI to be uninitialized");
  }
  if (!owns_mpi && !initialized) {
    throw ExceptionError("host-owned MPI launcher requires MPI to be initialized");
  }
  if (owns_mpi) {
    status = MPI_Init(nullptr, nullptr);
    if (status != MPI_SUCCESS) {
      throw ExceptionError("MPI_Init failed in standalone launcher");
    }
    initialized = 0;
    status = MPI_Initialized(&initialized);
    if (status != MPI_SUCCESS || !initialized) {
      throw ExceptionError("MPI_Init did not initialize MPI");
    }
  }
  int actual = 0;
  int rank = 0;
  status = MPI_Comm_size(MPI_COMM_WORLD, &actual);
  if (status == MPI_SUCCESS) {
    status = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  }
  if (status != MPI_SUCCESS || actual != ranks) {
    if (owns_mpi) {
      const int finalize_status = MPI_Finalize();
      if (finalize_status != MPI_SUCCESS) {
        throw ExceptionError("MPI_Finalize failed after rank-count validation");
      }
    }
    throw ExceptionError("MPI process count does not match requested rank count");
  }
  int local_result = 0;
  try {
    MPI_Comm world = MPI_COMM_WORLD;
    local_result = callback(rank, &world, sizeof(world), callback_context);
  } catch (...) {
    if (owns_mpi) {
      std::fprintf(stderr, "PLUMED standalone MPI rank %d callback failed; aborting MPI\n", rank);
      std::fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      // A provider returning from Abort must not leave peers blocked.
      std::_Exit(EXIT_FAILURE);
    }
    throw;
  }
  *result = local_result;
  return owns_mpi ? MPI_Finalize() : status;
}

}  // namespace mpiwrap_provider
}  // namespace PLMD

#elif defined(__PLUMED_MPIWRAP) && defined(__PLUMED_HAS_THREAD_MPI)
// Keep the internal rank launcher in a separate kernel translation unit.
#include "thread-mpi-launcher.cpp.inc"
#endif
