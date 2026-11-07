#include <tmpi.h>

#include <atomic>
#include <cstdlib>
#include <iostream>

namespace PLMD {
namespace thread_mpi {

std::atomic<int> failures{0};
std::atomic<int> ranks_seen{0};

void rank_work(const void *) {
  int rank = -1;
  int size = 0;
  int input = 0;
  int total = 0;
  int received = -1;
  if (tMPI_Comm_rank(TMPI_COMM_WORLD, &rank) != TMPI_SUCCESS ||
      tMPI_Comm_size(TMPI_COMM_WORLD, &size) != TMPI_SUCCESS || size < 1) {
    failures.fetch_add(1);
    return;
  }
  input = rank + 1;
  const int send_to = (rank + 1) % size;
  const int recv_from = (rank + size - 1) % size;
  if (tMPI_Sendrecv(&input, 1, TMPI_INT, send_to, 17,
                    &received, 1, TMPI_INT, recv_from, 17,
                    TMPI_COMM_WORLD, TMPI_STATUS_IGNORE) != TMPI_SUCCESS ||
      received != recv_from + 1) {
    failures.fetch_add(1);
  }
  // Upstream v2026.4 tMPI_Reduce_fast skips the copy when N == 1 and the
  // buffers differ. The wrapper regression will cover this required case after
  // its single-rank fallback is implemented; keep the pinned import pristine.
  if (size > 1 &&
      (tMPI_Allreduce(&input, &total, 1, TMPI_INT, TMPI_SUM,
                      TMPI_COMM_WORLD) != TMPI_SUCCESS ||
       total != size * (size + 1) / 2)) {
    failures.fetch_add(1);
  }
  ranks_seen.fetch_add(1);
}

int run_test(int argc, char **argv) {
  if (argc != 2) {
    return 2;
  }
  const int ranks = std::atoi(argv[1]);
  if (ranks < 1 || ranks > 4) {
    return 2;
  }
  const int rc = tMPI_Init_fn(0, ranks, TMPI_AFFINITY_NONE, rank_work, nullptr);
  if (rc != TMPI_SUCCESS) {
    std::cerr << "tMPI_Init_fn failed: " << rc << '\n';
    return 1;
  }
  int finalized = 0;
  if (tMPI_Finalized(&finalized) != TMPI_SUCCESS || !finalized ||
      ranks_seen.load() != ranks) {
    std::cerr << "launch returned before all ranks finalized: seen="
              << ranks_seen.load() << " expected=" << ranks
              << " finalized=" << finalized << '\n';
    return 1;
  }
  if (failures.load() != 0) {
    std::cerr << "rank checks failed: " << failures.load() << '\n';
    return 1;
  }
  return 0;
}

} // namespace thread_mpi
} // namespace PLMD

int main(int argc, char **argv) {
  return PLMD::thread_mpi::run_test(argc, argv);
}
