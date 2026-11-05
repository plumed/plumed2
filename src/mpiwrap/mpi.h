/* Minimal MPI runtime inquiries shared by PLUMED modules. */
#ifndef PLUMED_MPIWRAP_MPI_H
#define PLUMED_MPIWRAP_MPI_H

namespace PLMD {

#if defined(__PLUMED_HAS_MPI)
inline bool mpiwrap_provider_installed() {
  return true;
}
inline bool mpiwrap_provider_ranks_are_threads() {
  return false;
}
inline bool mpiwrap_is_thread_mpi_active() {
  return false;
}
#elif !defined(__PLUMED_MPIWRAP)
inline bool mpiwrap_provider_installed() {
  return false;
}
inline bool mpiwrap_provider_ranks_are_threads() {
  return false;
}
inline bool mpiwrap_is_thread_mpi_active() {
  return false;
}
#else
bool mpiwrap_provider_installed();
// Describes the selected provider without requiring MPI initialization.
bool mpiwrap_provider_ranks_are_threads();
bool mpiwrap_is_thread_mpi_active();
#endif

}  // namespace PLMD

#endif
