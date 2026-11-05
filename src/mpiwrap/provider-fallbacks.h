// Shared kernel-side fallbacks for provider MPI operations.
#ifndef PLUMED_MPIWRAP_PROVIDER_FALLBACKS_HPP
#define PLUMED_MPIWRAP_PROVIDER_FALLBACKS_HPP

#include <cstddef>

namespace PLMD {
struct MPI_Request;
struct MPI_Status;
struct MPI_Datatype;
struct MPI_Comm;

bool mpiwrap_can_fallback_Wait();
bool mpiwrap_can_fallback_Waitsome();
bool mpiwrap_can_fallback_Testsome();
bool mpiwrap_can_fallback_Allgather();
bool mpiwrap_can_fallback_Allgatherv();
int mpiwrap_fallback_Wait(MPI_Request *, MPI_Status *);
int mpiwrap_fallback_Waitsome(int, MPI_Request [], int *, int [], MPI_Status []);
int mpiwrap_fallback_Testsome(int, MPI_Request [], int *, int [], MPI_Status []);
int mpiwrap_fallback_Allgather(const void *, int, MPI_Datatype,
                               void *, int, MPI_Datatype, MPI_Comm);
int mpiwrap_fallback_Allgatherv(const void *, int, MPI_Datatype,
                                void *, const int [], const int [],
                                MPI_Datatype, MPI_Comm);
bool mpiwrap_is_empty_status(const MPI_Status *);
void mpiwrap_set_empty_status(MPI_Status *);
bool mpiwrap_requests_all_null(int, const MPI_Request []);
bool mpiwrap_try_singleton_reduce(const void *, void *, int, MPI_Datatype,
                                  int, MPI_Comm, int *);
}

#endif
