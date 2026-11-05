// Generated from the canonical native pointer-transport renderer.
#include "provider-abi.h"
#include "thread_mpi/thread_mpi_tmpi.h"
#include <cstddef>
#include <cstdint>
#include <cstring>
namespace PLMD {
namespace thread_mpi_adapter {
using namespace thread_mpi;
static void plumed_mpiw_MPI_Initialized(void *argument_0, void *result) {
  int native_result = thread_mpi::tMPI_Initialized((int *)argument_0);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Comm_dup(void *argument_0, void *argument_1, void *result) {
  thread_mpi::tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = thread_mpi::tMPI_Comm_dup(native_argument_0, (thread_mpi::tMPI_Comm *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Comm_free(void *argument_0, void *result) {
  int native_result = thread_mpi::tMPI_Comm_free((thread_mpi::tMPI_Comm *)argument_0);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Comm_rank(void *argument_0, void *argument_1, void *result) {
  thread_mpi::tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = thread_mpi::tMPI_Comm_rank(native_argument_0, (int *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Comm_size(void *argument_0, void *argument_1, void *result) {
  thread_mpi::tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = thread_mpi::tMPI_Comm_size(native_argument_0, (int *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Comm_split(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *result) {
  thread_mpi::tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_result = thread_mpi::tMPI_Comm_split(native_argument_0, native_argument_1, native_argument_2, (thread_mpi::tMPI_Comm *)argument_3);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Abort(void *argument_0, void *argument_1, void *result) {
  thread_mpi::tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = thread_mpi::tMPI_Abort(native_argument_0, native_argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Bcast(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  thread_mpi::tMPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  thread_mpi::tMPI_Comm native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  int native_result = thread_mpi::tMPI_Bcast((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Scatter(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *argument_7, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  thread_mpi::tMPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  thread_mpi::tMPI_Datatype native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_argument_6;
  memcpy(&native_argument_6, argument_6, sizeof(native_argument_6));
  thread_mpi::tMPI_Comm native_argument_7;
  memcpy(&native_argument_7, argument_7, sizeof(native_argument_7));
  int native_result = thread_mpi::tMPI_Scatter((void *)argument_0, native_argument_1, native_argument_2, (void *)argument_3, native_argument_4, native_argument_5, native_argument_6, native_argument_7);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Gather(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *argument_7, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  thread_mpi::tMPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  thread_mpi::tMPI_Datatype native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_argument_6;
  memcpy(&native_argument_6, argument_6, sizeof(native_argument_6));
  thread_mpi::tMPI_Comm native_argument_7;
  memcpy(&native_argument_7, argument_7, sizeof(native_argument_7));
  int native_result = thread_mpi::tMPI_Gather((void *)argument_0, native_argument_1, native_argument_2, (void *)argument_3, native_argument_4, native_argument_5, native_argument_6, native_argument_7);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Allreduce(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *result) {
  int native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  thread_mpi::tMPI_Datatype native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  thread_mpi::tMPI_Op native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  thread_mpi::tMPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = thread_mpi::tMPI_Allreduce((void *)argument_0, (void *)argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Isend(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  thread_mpi::tMPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  thread_mpi::tMPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = thread_mpi::tMPI_Isend((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5, (thread_mpi::tMPI_Request *)argument_6);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Recv(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  thread_mpi::tMPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  thread_mpi::tMPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = thread_mpi::tMPI_Recv((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5, (thread_mpi::tMPI_Status *)argument_6);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Barrier(void *argument_0, void *result) {
  thread_mpi::tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = thread_mpi::tMPI_Barrier(native_argument_0);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Get_count(void *argument_0, void *argument_1, void *argument_2, void *result) {
  thread_mpi::tMPI_Datatype native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = thread_mpi::tMPI_Get_count((thread_mpi::tMPI_Status *)argument_0, native_argument_1, (int *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Sendrecv(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *argument_7, void *argument_8, void *argument_9, void *argument_10, void *argument_11, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  thread_mpi::tMPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  int native_argument_6;
  memcpy(&native_argument_6, argument_6, sizeof(native_argument_6));
  thread_mpi::tMPI_Datatype native_argument_7;
  memcpy(&native_argument_7, argument_7, sizeof(native_argument_7));
  int native_argument_8;
  memcpy(&native_argument_8, argument_8, sizeof(native_argument_8));
  int native_argument_9;
  memcpy(&native_argument_9, argument_9, sizeof(native_argument_9));
  thread_mpi::tMPI_Comm native_argument_10;
  memcpy(&native_argument_10, argument_10, sizeof(native_argument_10));
  int native_result = thread_mpi::tMPI_Sendrecv((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, (void *)argument_5, native_argument_6, native_argument_7, native_argument_8, native_argument_9, native_argument_10, (thread_mpi::tMPI_Status *)argument_11);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Error_string(void *argument_0, void *argument_1, void *argument_2, void *result) {
  int native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = thread_mpi::tMPI_Error_string(native_argument_0, (char *)argument_1, (int *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Send(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  thread_mpi::tMPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  thread_mpi::tMPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = thread_mpi::tMPI_Send((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Rsend(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  thread_mpi::tMPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  thread_mpi::tMPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = thread_mpi::tMPI_Send((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Irsend(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  thread_mpi::tMPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  thread_mpi::tMPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = thread_mpi::tMPI_Isend((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5, (thread_mpi::tMPI_Request *)argument_6);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Irecv(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  thread_mpi::tMPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  thread_mpi::tMPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = thread_mpi::tMPI_Irecv((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5, (thread_mpi::tMPI_Request *)argument_6);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Waitany(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *result) {
  int native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = thread_mpi::tMPI_Waitany(native_argument_0, (thread_mpi::tMPI_Request *)argument_1, (int *)argument_2, (thread_mpi::tMPI_Status *)argument_3);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Waitall(void *argument_0, void *argument_1, void *argument_2, void *result) {
  int native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = thread_mpi::tMPI_Waitall(native_argument_0, (thread_mpi::tMPI_Request *)argument_1, (thread_mpi::tMPI_Status *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Type_contiguous(void *argument_0, void *argument_1, void *argument_2, void *result) {
  int native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  thread_mpi::tMPI_Datatype native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = thread_mpi::tMPI_Type_contiguous(native_argument_0, native_argument_1, (thread_mpi::tMPI_Datatype *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Type_commit(void *argument_0, void *result) {
  int native_result = thread_mpi::tMPI_Type_commit((thread_mpi::tMPI_Datatype *)argument_0);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Type_free(void *argument_0, void *result) {
  int native_result = thread_mpi::tMPI_Type_free((thread_mpi::tMPI_Datatype *)argument_0);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Gatherv(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *argument_7, void *argument_8, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  thread_mpi::tMPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  thread_mpi::tMPI_Datatype native_argument_6;
  memcpy(&native_argument_6, argument_6, sizeof(native_argument_6));
  int native_argument_7;
  memcpy(&native_argument_7, argument_7, sizeof(native_argument_7));
  thread_mpi::tMPI_Comm native_argument_8;
  memcpy(&native_argument_8, argument_8, sizeof(native_argument_8));
  int native_result = thread_mpi::tMPI_Gatherv((void *)argument_0, native_argument_1, native_argument_2, (void *)argument_3, (int *)argument_4, (int *)argument_5, native_argument_6, native_argument_7, native_argument_8);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Scatterv(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *argument_7, void *argument_8, void *result) {
  thread_mpi::tMPI_Datatype native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  thread_mpi::tMPI_Datatype native_argument_6;
  memcpy(&native_argument_6, argument_6, sizeof(native_argument_6));
  int native_argument_7;
  memcpy(&native_argument_7, argument_7, sizeof(native_argument_7));
  thread_mpi::tMPI_Comm native_argument_8;
  memcpy(&native_argument_8, argument_8, sizeof(native_argument_8));
  int native_result = thread_mpi::tMPI_Scatterv((void *)argument_0, (int *)argument_1, (int *)argument_2, native_argument_3, (void *)argument_4, native_argument_5, native_argument_6, native_argument_7, native_argument_8);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Alltoall(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  thread_mpi::tMPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  thread_mpi::tMPI_Datatype native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  thread_mpi::tMPI_Comm native_argument_6;
  memcpy(&native_argument_6, argument_6, sizeof(native_argument_6));
  int native_result = thread_mpi::tMPI_Alltoall((void *)argument_0, native_argument_1, native_argument_2, (void *)argument_3, native_argument_4, native_argument_5, native_argument_6);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Alltoallv(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *argument_7, void *argument_8, void *result) {
  thread_mpi::tMPI_Datatype native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  thread_mpi::tMPI_Datatype native_argument_7;
  memcpy(&native_argument_7, argument_7, sizeof(native_argument_7));
  thread_mpi::tMPI_Comm native_argument_8;
  memcpy(&native_argument_8, argument_8, sizeof(native_argument_8));
  int native_result = thread_mpi::tMPI_Alltoallv((void *)argument_0, (int *)argument_1, (int *)argument_2, native_argument_3, (void *)argument_4, (int *)argument_5, (int *)argument_6, native_argument_7, native_argument_8);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Reduce(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *result) {
  int native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  thread_mpi::tMPI_Datatype native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  thread_mpi::tMPI_Op native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  int native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  thread_mpi::tMPI_Comm native_argument_6;
  memcpy(&native_argument_6, argument_6, sizeof(native_argument_6));
  int native_result = thread_mpi::tMPI_Reduce((void *)argument_0, (void *)argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5, native_argument_6);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Scan(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *result) {
  int native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  thread_mpi::tMPI_Datatype native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  thread_mpi::tMPI_Op native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  thread_mpi::tMPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = thread_mpi::tMPI_Scan((void *)argument_0, (void *)argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Group_size(void *argument_0, void *argument_1, void *result) {
  thread_mpi::tMPI_Group native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = thread_mpi::tMPI_Group_size(native_argument_0, (int *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Group_rank(void *argument_0, void *argument_1, void *result) {
  thread_mpi::tMPI_Group native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = thread_mpi::tMPI_Group_rank(native_argument_0, (int *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Comm_group(void *argument_0, void *argument_1, void *result) {
  thread_mpi::tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = thread_mpi::tMPI_Comm_group(native_argument_0, (thread_mpi::tMPI_Group *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Group_incl(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *result) {
  thread_mpi::tMPI_Group native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = thread_mpi::tMPI_Group_incl(native_argument_0, native_argument_1, (int *)argument_2, (thread_mpi::tMPI_Group *)argument_3);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Group_free(void *argument_0, void *result) {
  int native_result = thread_mpi::tMPI_Group_free((thread_mpi::tMPI_Group *)argument_0);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Comm_compare(void *argument_0, void *argument_1, void *argument_2, void *result) {
  thread_mpi::tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  thread_mpi::tMPI_Comm native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = thread_mpi::tMPI_Comm_compare(native_argument_0, native_argument_1, (int *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Comm_create(void *argument_0, void *argument_1, void *argument_2, void *result) {
  thread_mpi::tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  thread_mpi::tMPI_Group native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = thread_mpi::tMPI_Comm_create(native_argument_0, native_argument_1, (thread_mpi::tMPI_Comm *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Cart_create(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *result) {
  thread_mpi::tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  int native_result = thread_mpi::tMPI_Cart_create(native_argument_0, native_argument_1, (int *)argument_2, (int *)argument_3, native_argument_4, (thread_mpi::tMPI_Comm *)argument_5);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Topo_test(void *argument_0, void *argument_1, void *result) {
  thread_mpi::tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = thread_mpi::tMPI_Topo_test(native_argument_0, (int *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Cartdim_get(void *argument_0, void *argument_1, void *result) {
  thread_mpi::tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = thread_mpi::tMPI_Cartdim_get(native_argument_0, (int *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Cart_get(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *result) {
  thread_mpi::tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = thread_mpi::tMPI_Cart_get(native_argument_0, native_argument_1, (int *)argument_2, (int *)argument_3, (int *)argument_4);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Cart_rank(void *argument_0, void *argument_1, void *argument_2, void *result) {
  thread_mpi::tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = thread_mpi::tMPI_Cart_rank(native_argument_0, (int *)argument_1, (int *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Cart_coords(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *result) {
  thread_mpi::tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_result = thread_mpi::tMPI_Cart_coords(native_argument_0, native_argument_1, native_argument_2, (int *)argument_3);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Cart_sub(void *argument_0, void *argument_1, void *argument_2, void *result) {
  thread_mpi::tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = thread_mpi::tMPI_Cart_sub(native_argument_0, (int *)argument_1, (thread_mpi::tMPI_Comm *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Cart_map(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *result) {
  thread_mpi::tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = thread_mpi::tMPI_Cart_map(native_argument_0, native_argument_1, (int *)argument_2, (int *)argument_3, (int *)argument_4);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Get_processor_name(void *argument_0, void *argument_1, void *result) {
  int native_result = thread_mpi::tMPI_Get_processor_name((char *)argument_0, (int *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Wtime(void *result) {
  double native_result = thread_mpi::tMPI_Wtime();
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Wtick(void *result) {
  double native_result = thread_mpi::tMPI_Wtick();
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Finalized(void *argument_0, void *result) {
  int native_result = thread_mpi::tMPI_Finalized((int *)argument_0);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Comm_get_errhandler(void *argument_0, void *argument_1, void *result) {
  thread_mpi::tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = thread_mpi::tMPI_Comm_get_errhandler(native_argument_0, (thread_mpi::tMPI_Errhandler *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
static void plumed_mpiw_MPI_Comm_set_errhandler(void *argument_0, void *argument_1, void *result) {
  thread_mpi::tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  thread_mpi::tMPI_Errhandler native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = thread_mpi::tMPI_Comm_set_errhandler(native_argument_0, native_argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#include "thread-mpi-constant-accessors.inc"
static const plumed_mpiw_function_entry functions[] = {
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Initialized", 0xd4498369u, 1u, 1u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Initialized },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Comm_dup", 0xf0dc6608u, 1u, 2u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Comm_dup },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Comm_free", 0x4da4e1b2u, 1u, 1u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Comm_free },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Comm_rank", 0x4e7672aau, 1u, 2u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Comm_rank },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Comm_size", 0xda2369ecu, 1u, 2u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Comm_size },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Comm_split", 0x97f1adf0u, 1u, 4u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Comm_split },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Abort", 0x84f458d6u, 1u, 2u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Abort },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Bcast", 0xe809e1cau, 1u, 5u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Bcast },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Scatter", 0xb2fa0f58u, 1u, 8u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Scatter },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Gather", 0xb1f903d9u, 1u, 8u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Gather },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Allreduce", 0x788900dbu, 1u, 6u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Allreduce },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Isend", 0xf5e9d915u, 1u, 7u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Isend },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Recv", 0xdd56a4bcu, 1u, 7u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Recv },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Barrier", 0xe7a052f4u, 1u, 1u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Barrier },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Get_count", 0x71967a66u, 1u, 3u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Get_count },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Sendrecv", 0xb09a6fc1u, 1u, 12u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Sendrecv },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Error_string", 0xd272abecu, 1u, 3u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Error_string },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Send", 0xbd2b6912u, 1u, 6u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Send },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Rsend", 0xbfa09b06u, 1u, 6u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Rsend },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Irsend", 0xb13d1d70u, 1u, 7u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Irsend },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Irecv", 0xb69a0083u, 1u, 7u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Irecv },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Waitany", 0x5acd4436u, 1u, 4u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Waitany },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Waitall", 0x7a57ae23u, 1u, 3u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Waitall },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Type_contiguous", 0xe41681d2u, 1u, 3u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Type_contiguous },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Type_commit", 0xf848da1fu, 1u, 1u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Type_commit },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Type_free", 0xdacde819u, 1u, 1u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Type_free },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Gatherv", 0xc90455c9u, 1u, 9u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Gatherv },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Scatterv", 0x765eb673u, 1u, 9u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Scatterv },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Alltoall", 0xcfe270aeu, 1u, 7u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Alltoall },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Alltoallv", 0xa4d99ed0u, 1u, 9u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Alltoallv },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Reduce", 0x3a815604u, 1u, 7u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Reduce },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Scan", 0x64e2ee54u, 1u, 6u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Scan },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Group_size", 0x7de536d5u, 1u, 2u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Group_size },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Group_rank", 0x668fab96u, 1u, 2u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Group_rank },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Comm_group", 0xd11f4483u, 1u, 2u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Comm_group },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Group_incl", 0xbe8bfe71u, 1u, 4u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Group_incl },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Group_free", 0xfd1c62d1u, 1u, 1u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Group_free },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Comm_compare", 0xdacd79cdu, 1u, 3u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Comm_compare },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Comm_create", 0x1de85fa0u, 1u, 3u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Comm_create },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Cart_create", 0xfd363844u, 1u, 6u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Cart_create },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Topo_test", 0x09cbe810u, 1u, 2u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Topo_test },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Cartdim_get", 0xa979e88cu, 1u, 2u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Cartdim_get },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Cart_get", 0x0a3d9750u, 1u, 5u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Cart_get },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Cart_rank", 0xf899fe26u, 1u, 3u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Cart_rank },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Cart_coords", 0xf0bdc30cu, 1u, 4u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Cart_coords },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Cart_sub", 0x1e9d3861u, 1u, 3u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Cart_sub },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Cart_map", 0xba5a7b4cu, 1u, 5u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Cart_map },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Get_processor_name", 0xf8e4ccc3u, 1u, 2u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Get_processor_name },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Wtime", 0x75704e2eu, 1u, 0u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Wtime },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Wtick", 0x84f54b80u, 1u, 0u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Wtick },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Finalized", 0x105b65c7u, 1u, 1u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Finalized },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Comm_get_errhandler", 0x74d05159u, 1u, 2u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Comm_get_errhandler },
  { (uint32_t)sizeof(plumed_mpiw_function_entry), "MPI_Comm_set_errhandler", 0xe6af3a8fu, 1u, 2u, (plumed_mpiw_generic_fn)plumed_mpiw_MPI_Comm_set_errhandler },
};
static plumed_mpiw_generic_fn runtime_identity = (plumed_mpiw_generic_fn)&thread_mpi::tMPI_Initialized;
static const plumed_mpiw_metadata_entry metadata[] = {
#include "thread-mpi-constant-metadata.inc"
  { (uint32_t)sizeof(plumed_mpiw_metadata_entry), PLUMED_MPIW_VALUE_OPAQUE, 0, "identity:runtime", sizeof(plumed_mpiw_generic_fn), 0, &runtime_identity },
  { (uint32_t)sizeof(plumed_mpiw_metadata_entry), PLUMED_MPIW_VALUE_UTF8, 0, "capability:intracommunicators-only", sizeof("true"), 0, "true" },
};
static const plumed_mpiw_provider_descriptor descriptor = {
  (uint32_t)sizeof(plumed_mpiw_provider_descriptor), PLUMED_MPIW_PROVIDER_ABI_MAJOR, 0, 1,
  "PLUMED embedded GROMACS thread-MPI",
  (uint32_t)sizeof(thread_mpi::tMPI_Comm), (uint32_t)sizeof(thread_mpi::tMPI_Group),
  (uint32_t)sizeof(thread_mpi::tMPI_Datatype), (uint32_t)sizeof(thread_mpi::tMPI_Request),
  (uint32_t)sizeof(thread_mpi::tMPI_Op), (uint32_t)sizeof(thread_mpi::tMPI_Errhandler),
  (uint32_t)sizeof(thread_mpi::tMPI_Status),
  (uint32_t)offsetof(thread_mpi::tMPI_Status, TMPI_SOURCE),
  (uint32_t)offsetof(thread_mpi::tMPI_Status, TMPI_TAG),
  (uint32_t)offsetof(thread_mpi::tMPI_Status, TMPI_ERROR), 0, 0,
  (uint64_t)(sizeof(functions) / sizeof(functions[0])), functions,
  (uint64_t)(sizeof(metadata) / sizeof(metadata[0])), metadata
};
}
extern "C" const plumed_mpiw_provider_descriptor *plumed_mpiw_thread_mpi_provider(void) {
  return &thread_mpi_adapter::descriptor;
}
}
