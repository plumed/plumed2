#include "mpiwrap.h"
#if defined(__PLUMED_MPIWRAP) && defined(__PLUMED_HAS_MPI)
#include "../tools/Exception.h"

#include <mpi.h>

#include <cstdio>
#include <cstdlib>
#include <exception>

namespace PLMD {
namespace mpiwrap {

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

}  // namespace mpiwrap
}  // namespace PLMD

#elif defined(__PLUMED_MPIWRAP) && defined(__PLUMED_HAS_THREAD_MPI)
// Explicit launcher for the embedded GROMACS thread-MPI runtime.
#include "internal.h"

#include <thread_mpi_tmpi.h>
#include "../tools/Exception.h"

#include <cstdio>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <string>

namespace PLMD {
namespace mpiwrap {
extern "C" const plumed_mpiw_provider_descriptor *
plumed_mpiw_thread_mpi_provider(void);
} // namespace mpiwrap
}

namespace PLMD {
namespace mpiwrap {
namespace {

using Callback = int (*)(int, const void *, std::size_t, void *);

struct LaunchContext {
  Callback callback;
  void *callback_context;
  int process_result;
};

[[noreturn]] void terminateRank(int rank, const char *message) {
  std::fprintf(stderr, "PLUMED thread-MPI rank %d failed: %s\n", rank,
               message ? message : "unknown exception");
  std::fflush(stderr);
  std::_Exit(EXIT_FAILURE);
}

void rankEntry(const void *opaque) {
  LaunchContext *context = const_cast<LaunchContext *>(
                             static_cast<const LaunchContext *>(opaque));
  int rank = -1;
  if (thread_mpi::tMPI_Comm_rank(thread_mpi::TMPI_COMM_WORLD, &rank) !=
      thread_mpi::TMPI_SUCCESS) {
    terminateRank(rank, "could not determine rank");
  }

  int local_result = 0;
  try {
    alignas(std::max_align_t)
    unsigned char native_world[sizeof(plumed_mpiw_handle_storage)] {};
    const MPI_Comm world = COMM_WORLD();
    const std::size_t native_size = currentRegistry().nativeCommSize();
    if (!Comm_export_native(world, native_world,
                            sizeof(native_world)) || native_size == 0) {
      terminateRank(rank, "could not export the native world communicator");
    }
    local_result = context->callback(rank, native_world, native_size,
                                     context->callback_context);
  } catch (const std::exception &error) {
    terminateRank(rank, error.what());
  } catch (...) {
    terminateRank(rank, "uncaught non-standard exception in rank callback");
  }

  int process_result = local_result;
  int size = 0;
  if (thread_mpi::tMPI_Comm_size(thread_mpi::TMPI_COMM_WORLD, &size) !=
      thread_mpi::TMPI_SUCCESS) {
    terminateRank(rank, "could not determine world size");
  }
  // The imported upstream Allreduce has a known one-rank defect; the
  // identity result is already local in that case.
  if (size > 1 && thread_mpi::tMPI_Allreduce(
        &local_result, &process_result, 1, thread_mpi::TMPI_INT,
        thread_mpi::TMPI_MAX, thread_mpi::TMPI_COMM_WORLD) !=
      thread_mpi::TMPI_SUCCESS) {
    terminateRank(rank, "could not combine rank return codes");
  }
  if (rank == 0) {
    context->process_result = process_result;
  }
}

}  // namespace

int launchInternalThreadMPI(int ranks, Callback callback,
                            void *callback_context, int *result) {
  if (ranks < 1 || !callback || !result) {
    throw ExceptionTypeError("rank count, callback, and result output are required");
  }

  forceEnvironmentProvider();
  const plumed_mpiw_provider_descriptor *descriptor =
    plumed_mpiw_thread_mpi_provider();
  Registry &registry = currentRegistry();
  std::string error;
  if (!registry.install(descriptor, error)) {
    throw ExceptionError(error);
  }
  if (!registry.beginInternalLaunch(descriptor, error)) {
    throw ExceptionError(error);
  }

  LaunchContext launch_context{callback, callback_context, 0};
  const int status = thread_mpi::tMPI_Init_fn(
                       0, ranks, thread_mpi::TMPI_AFFINITY_NONE,
                       &rankEntry, &launch_context);
  if (status != thread_mpi::TMPI_SUCCESS) {
    std::fprintf(stderr, "PLUMED thread-MPI initialization failed (%d)\n", status);
    std::fflush(stderr);
    std::_Exit(EXIT_FAILURE);
  }
  registry.finishInternalLaunch();
  *result = launch_context.process_result;
  return 0;
}

int launchRanks(int ranks, Callback callback, void *callback_context,
                int *result, MPILifecyclePolicy lifecycle) {
  if (lifecycle == MPILifecyclePolicy::internal_thread_mpi) {
    return launchInternalThreadMPI(ranks, callback, callback_context, result);
  }
  (void)ranks;
  (void)callback;
  (void)callback_context;
  (void)result;
  throw ExceptionError("real-MPI lifecycle policy is unavailable in an internal tMPI build");
}

}  // namespace mpiwrap
}  // namespace PLMD

#endif
