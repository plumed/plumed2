#include <tmpi.h>

#include <atomic>

namespace PLMD {
namespace thread_mpi_test {
namespace host_copy {

std::atomic<int> visited{0};
std::atomic<int> failures{0};

void rank_work(const void *) {
  int rank = -1;
  int size = 0;
  if (tMPI_Comm_rank(TMPI_COMM_WORLD, &rank) != TMPI_SUCCESS ||
      tMPI_Comm_size(TMPI_COMM_WORLD, &size) != TMPI_SUCCESS || size != 2 ||
      rank < 0 || rank >= size) {
    failures.fetch_add(1);
  }
  visited.fetch_add(1);
}

int run() {
  if (tMPI_Init_fn(0, 2, TMPI_AFFINITY_NONE, rank_work, nullptr) != TMPI_SUCCESS) {
    return 1;
  }
  int finalized = 0;
  return failures.load() == 0 && visited.load() == 2 &&
         tMPI_Finalized(&finalized) == TMPI_SUCCESS && finalized ? 0 : 1;
}

} // namespace host_copy
} // namespace thread_mpi_test
} // namespace PLMD
