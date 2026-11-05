// Explicit launcher for the embedded GROMACS thread-MPI runtime.
#include "thread-mpi-launcher.h"
#include "provider-registry.hpp"

#include "provider-abi.h"
#include <thread_mpi_tmpi.h>
#include "provider-loader.cpp.inc"
#include "../tools/Exception.h"

#include <cstdio>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <string>
#include <condition_variable>
#include <mutex>
#include <thread>

namespace PLMD {
extern "C" const plumed_mpiw_provider_descriptor *
plumed_mpiw_thread_mpi_provider(void);
}

namespace PLMD {
namespace mpiwrap_provider {
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
    const MPI_Comm world = mpiwrap_COMM_WORLD();
    const std::size_t native_size = currentRegistry().nativeCommSize();
    if (!mpiwrap_Comm_export_native(world, native_world,
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

#if !defined(__PLUMED_HAS_MPI)
void ensureProviderSelected() {
  if (forceEnvironmentProvider()) {
    return;
  }
  static std::mutex selection_mutex;
  static std::condition_variable selection_changed;
  static bool selecting = false;
  static std::thread::id selecting_thread;

  Registry &registry = currentRegistry();
  if (registry.installed()) {
    return;
  }
  {
    std::unique_lock<std::mutex> lock(selection_mutex);
    // Another thread may install a provider while we wait for the lock.
    // cppcheck-suppress identicalConditionAfterEarlyExit
    if (registry.installed()) {
      return;
    }
    while (selecting) {
      if (selecting_thread == std::this_thread::get_id()) {
        throw std::runtime_error("recursive MPI provider selection");
      }
      selection_changed.wait(lock);
      if (registry.installed()) {
        return;
      }
    }
    selecting = true;
    selecting_thread = std::this_thread::get_id();
  }

  try {
    const char *selected = std::getenv("PLUMED_MPI_WRAPPER");
    if (selected && selected[0]) {
      installEnvironmentProvider(registry, selected);
    } else {
      const plumed_mpiw_provider_descriptor *descriptor =
        plumed_mpiw_thread_mpi_provider();
      std::string error;
      if (!registry.install(descriptor, error) && !registry.installed()) {
        throw std::runtime_error(error.empty() ?
                                 "could not select the default MPI provider" : error);
      }
    }
  } catch (...) {
    std::lock_guard<std::mutex> lock(selection_mutex);
    selecting = false;
    selecting_thread = std::thread::id();
    selection_changed.notify_all();
    throw;
  }
  {
    std::lock_guard<std::mutex> lock(selection_mutex);
    selecting = false;
    selecting_thread = std::thread::id();
  }
  selection_changed.notify_all();
}

#endif

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

}  // namespace mpiwrap_provider
}  // namespace PLMD
