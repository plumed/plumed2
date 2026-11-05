namespace PLMD {
namespace thread_mpi_test {
namespace host_copy {
int run();
}
namespace internal_copy {
int run();
}
} // namespace thread_mpi_test
} // namespace PLMD

int main() {
  if (PLMD::thread_mpi_test::host_copy::run() != 0) {
    return 1;
  }
  return PLMD::thread_mpi_test::internal_copy::run();
}
