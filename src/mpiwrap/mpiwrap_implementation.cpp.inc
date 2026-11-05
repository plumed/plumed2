/* Minimal native-MPI provider adapter. All conversions happen in the kernel. */

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
/* For standalone host builds without gromacs/version.h, callers may define the verified numeric version explicitly. */
#ifndef PLUMED_MPIW_GROMACS_VERSION
#ifdef GMX_VERSION
#define PLUMED_MPIW_GROMACS_VERSION GMX_VERSION
#else
#error "GROMACS tMPI host adapter requires GMX_VERSION or PLUMED_MPIW_GROMACS_VERSION"
#endif
#endif
#if PLUMED_MPIW_GROMACS_VERSION < 20260000
#error "GROMACS tMPI host adapter requires GROMACS 2026.0 or newer"
#endif
#endif

#define PLUMED_MPIW_JOIN_INNER(a, b) a##b
#define PLUMED_MPIW_JOIN(a, b) PLUMED_MPIW_JOIN_INNER(a, b)
#define PLUMED_MPIW_STATIC_ASSERT(c, m) \
  typedef char PLUMED_MPIW_JOIN(plumed_mpiw_static_assert_, __LINE__)[(c) ? 1 : -1]

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Initialized(void *argument_0, void *result) {
  int native_result = tMPI_Initialized((int *)argument_0);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Initialized(void *argument_0, void *result) {
  int native_result = MPI_Initialized((int *)argument_0);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Comm_dup(void *argument_0, void *argument_1, void *result) {
  tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = tMPI_Comm_dup(native_argument_0, (tMPI_Comm *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Comm_dup(void *argument_0, void *argument_1, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = MPI_Comm_dup(native_argument_0, (MPI_Comm *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Comm_free(void *argument_0, void *result) {
  int native_result = tMPI_Comm_free((tMPI_Comm *)argument_0);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Comm_free(void *argument_0, void *result) {
  int native_result = MPI_Comm_free((MPI_Comm *)argument_0);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Comm_rank(void *argument_0, void *argument_1, void *result) {
  tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = tMPI_Comm_rank(native_argument_0, (int *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Comm_rank(void *argument_0, void *argument_1, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = MPI_Comm_rank(native_argument_0, (int *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Comm_size(void *argument_0, void *argument_1, void *result) {
  tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = tMPI_Comm_size(native_argument_0, (int *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Comm_size(void *argument_0, void *argument_1, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = MPI_Comm_size(native_argument_0, (int *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Comm_split(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *result) {
  tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_result = tMPI_Comm_split(native_argument_0, native_argument_1, native_argument_2, (tMPI_Comm *)argument_3);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Comm_split(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_result = MPI_Comm_split(native_argument_0, native_argument_1, native_argument_2, (MPI_Comm *)argument_3);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Abort(void *argument_0, void *argument_1, void *result) {
  tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = tMPI_Abort(native_argument_0, native_argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Abort(void *argument_0, void *argument_1, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = MPI_Abort(native_argument_0, native_argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Bcast(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  tMPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  tMPI_Comm native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  int native_result = tMPI_Bcast((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Bcast(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  MPI_Comm native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  int native_result = MPI_Bcast((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Scatter(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *argument_7, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  tMPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  tMPI_Datatype native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_argument_6;
  memcpy(&native_argument_6, argument_6, sizeof(native_argument_6));
  tMPI_Comm native_argument_7;
  memcpy(&native_argument_7, argument_7, sizeof(native_argument_7));
  int native_result = tMPI_Scatter((void *)argument_0, native_argument_1, native_argument_2, (void *)argument_3, native_argument_4, native_argument_5, native_argument_6, native_argument_7);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Scatter(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *argument_7, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  MPI_Datatype native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_argument_6;
  memcpy(&native_argument_6, argument_6, sizeof(native_argument_6));
  MPI_Comm native_argument_7;
  memcpy(&native_argument_7, argument_7, sizeof(native_argument_7));
  int native_result = MPI_Scatter((void *)argument_0, native_argument_1, native_argument_2, (void *)argument_3, native_argument_4, native_argument_5, native_argument_6, native_argument_7);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Gather(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *argument_7, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  tMPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  tMPI_Datatype native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_argument_6;
  memcpy(&native_argument_6, argument_6, sizeof(native_argument_6));
  tMPI_Comm native_argument_7;
  memcpy(&native_argument_7, argument_7, sizeof(native_argument_7));
  int native_result = tMPI_Gather((void *)argument_0, native_argument_1, native_argument_2, (void *)argument_3, native_argument_4, native_argument_5, native_argument_6, native_argument_7);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Gather(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *argument_7, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  MPI_Datatype native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_argument_6;
  memcpy(&native_argument_6, argument_6, sizeof(native_argument_6));
  MPI_Comm native_argument_7;
  memcpy(&native_argument_7, argument_7, sizeof(native_argument_7));
  int native_result = MPI_Gather((void *)argument_0, native_argument_1, native_argument_2, (void *)argument_3, native_argument_4, native_argument_5, native_argument_6, native_argument_7);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Allreduce(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *result) {
  int native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  tMPI_Datatype native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  tMPI_Op native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  tMPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = tMPI_Allreduce((void *)argument_0, (void *)argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Allreduce(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *result) {
  int native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  MPI_Datatype native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  MPI_Op native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  MPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = MPI_Allreduce((void *)argument_0, (void *)argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Isend(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  tMPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  tMPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = tMPI_Isend((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5, (tMPI_Request *)argument_6);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Isend(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  MPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = MPI_Isend((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5, (MPI_Request *)argument_6);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Wait(void *argument_0, void *argument_1, void *result) {
  int native_result = MPI_Wait((MPI_Request *)argument_0, (MPI_Status *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Allgather(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  MPI_Datatype native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  MPI_Comm native_argument_6;
  memcpy(&native_argument_6, argument_6, sizeof(native_argument_6));
  int native_result = MPI_Allgather((void *)argument_0, native_argument_1, native_argument_2, (void *)argument_3, native_argument_4, native_argument_5, native_argument_6);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Allgatherv(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *argument_7, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  MPI_Datatype native_argument_6;
  memcpy(&native_argument_6, argument_6, sizeof(native_argument_6));
  MPI_Comm native_argument_7;
  memcpy(&native_argument_7, argument_7, sizeof(native_argument_7));
  int native_result = MPI_Allgatherv((void *)argument_0, native_argument_1, native_argument_2, (void *)argument_3, (int *)argument_4, (int *)argument_5, native_argument_6, native_argument_7);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Recv(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  tMPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  tMPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = tMPI_Recv((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5, (tMPI_Status *)argument_6);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Recv(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  MPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = MPI_Recv((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5, (MPI_Status *)argument_6);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Barrier(void *argument_0, void *result) {
  tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = tMPI_Barrier(native_argument_0);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Barrier(void *argument_0, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = MPI_Barrier(native_argument_0);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Get_count(void *argument_0, void *argument_1, void *argument_2, void *result) {
  tMPI_Datatype native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = tMPI_Get_count((tMPI_Status *)argument_0, native_argument_1, (int *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Get_count(void *argument_0, void *argument_1, void *argument_2, void *result) {
  MPI_Datatype native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = MPI_Get_count((MPI_Status *)argument_0, native_argument_1, (int *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Sendrecv(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *argument_7, void *argument_8, void *argument_9, void *argument_10, void *argument_11, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  tMPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  int native_argument_6;
  memcpy(&native_argument_6, argument_6, sizeof(native_argument_6));
  tMPI_Datatype native_argument_7;
  memcpy(&native_argument_7, argument_7, sizeof(native_argument_7));
  int native_argument_8;
  memcpy(&native_argument_8, argument_8, sizeof(native_argument_8));
  int native_argument_9;
  memcpy(&native_argument_9, argument_9, sizeof(native_argument_9));
  tMPI_Comm native_argument_10;
  memcpy(&native_argument_10, argument_10, sizeof(native_argument_10));
  int native_result = tMPI_Sendrecv((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, (void *)argument_5, native_argument_6, native_argument_7, native_argument_8, native_argument_9, native_argument_10, (tMPI_Status *)argument_11);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Sendrecv(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *argument_7, void *argument_8, void *argument_9, void *argument_10, void *argument_11, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  int native_argument_6;
  memcpy(&native_argument_6, argument_6, sizeof(native_argument_6));
  MPI_Datatype native_argument_7;
  memcpy(&native_argument_7, argument_7, sizeof(native_argument_7));
  int native_argument_8;
  memcpy(&native_argument_8, argument_8, sizeof(native_argument_8));
  int native_argument_9;
  memcpy(&native_argument_9, argument_9, sizeof(native_argument_9));
  MPI_Comm native_argument_10;
  memcpy(&native_argument_10, argument_10, sizeof(native_argument_10));
  int native_result = MPI_Sendrecv((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, (void *)argument_5, native_argument_6, native_argument_7, native_argument_8, native_argument_9, native_argument_10, (MPI_Status *)argument_11);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Error_string(void *argument_0, void *argument_1, void *argument_2, void *result) {
  int native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = tMPI_Error_string(native_argument_0, (char *)argument_1, (int *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Error_string(void *argument_0, void *argument_1, void *argument_2, void *result) {
  int native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = MPI_Error_string(native_argument_0, (char *)argument_1, (int *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 0)))
static void plumed_mpiw_MPI_Comm_f2c(void *argument_0, void *result) {
  MPI_Fint native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  MPI_Comm native_result = MPI_Comm_f2c(native_argument_0);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Send(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  tMPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  tMPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = tMPI_Send((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Send(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  MPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = MPI_Send((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Bsend(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  MPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = MPI_Bsend((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Ssend(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  MPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = MPI_Ssend((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Rsend(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  tMPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  tMPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = tMPI_Send((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Rsend(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  MPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = MPI_Rsend((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Buffer_attach(void *argument_0, void *argument_1, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = MPI_Buffer_attach((void *)argument_0, native_argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Buffer_detach(void *argument_0, void *argument_1, void *result) {
  int native_result = MPI_Buffer_detach((void *)argument_0, (int *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Ibsend(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  MPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = MPI_Ibsend((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5, (MPI_Request *)argument_6);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Issend(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  MPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = MPI_Issend((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5, (MPI_Request *)argument_6);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Irsend(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  tMPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  tMPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = tMPI_Isend((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5, (tMPI_Request *)argument_6);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Irsend(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  MPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = MPI_Irsend((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5, (MPI_Request *)argument_6);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Irecv(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  tMPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  tMPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = tMPI_Irecv((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5, (tMPI_Request *)argument_6);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Irecv(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  MPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = MPI_Irecv((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5, (MPI_Request *)argument_6);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Test(void *argument_0, void *argument_1, void *argument_2, void *result) {
  int native_result = MPI_Test((MPI_Request *)argument_0, (int *)argument_1, (MPI_Status *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Request_free(void *argument_0, void *result) {
  int native_result = MPI_Request_free((MPI_Request *)argument_0);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Waitany(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *result) {
  int native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = tMPI_Waitany(native_argument_0, (tMPI_Request *)argument_1, (int *)argument_2, (tMPI_Status *)argument_3);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Waitany(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *result) {
  int native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = MPI_Waitany(native_argument_0, (MPI_Request *)argument_1, (int *)argument_2, (MPI_Status *)argument_3);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Testany(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *result) {
  int native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = MPI_Testany(native_argument_0, (MPI_Request *)argument_1, (int *)argument_2, (int *)argument_3, (MPI_Status *)argument_4);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Waitall(void *argument_0, void *argument_1, void *argument_2, void *result) {
  int native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = tMPI_Waitall(native_argument_0, (tMPI_Request *)argument_1, (tMPI_Status *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Waitall(void *argument_0, void *argument_1, void *argument_2, void *result) {
  int native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = MPI_Waitall(native_argument_0, (MPI_Request *)argument_1, (MPI_Status *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Testall(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *result) {
  int native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = MPI_Testall(native_argument_0, (MPI_Request *)argument_1, (int *)argument_2, (MPI_Status *)argument_3);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Waitsome(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *result) {
  int native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = MPI_Waitsome(native_argument_0, (MPI_Request *)argument_1, (int *)argument_2, (int *)argument_3, (MPI_Status *)argument_4);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Testsome(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *result) {
  int native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = MPI_Testsome(native_argument_0, (MPI_Request *)argument_1, (int *)argument_2, (int *)argument_3, (MPI_Status *)argument_4);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Iprobe(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *result) {
  int native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Comm native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_result = MPI_Iprobe(native_argument_0, native_argument_1, native_argument_2, (int *)argument_3, (MPI_Status *)argument_4);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Probe(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *result) {
  int native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Comm native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_result = MPI_Probe(native_argument_0, native_argument_1, native_argument_2, (MPI_Status *)argument_3);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Cancel(void *argument_0, void *result) {
  int native_result = MPI_Cancel((MPI_Request *)argument_0);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Test_cancelled(void *argument_0, void *argument_1, void *result) {
  int native_result = MPI_Test_cancelled((MPI_Status *)argument_0, (int *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Send_init(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  MPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = MPI_Send_init((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5, (MPI_Request *)argument_6);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Bsend_init(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  MPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = MPI_Bsend_init((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5, (MPI_Request *)argument_6);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Ssend_init(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  MPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = MPI_Ssend_init((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5, (MPI_Request *)argument_6);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Rsend_init(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  MPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = MPI_Rsend_init((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5, (MPI_Request *)argument_6);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Recv_init(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  MPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = MPI_Recv_init((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5, (MPI_Request *)argument_6);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Start(void *argument_0, void *result) {
  int native_result = MPI_Start((MPI_Request *)argument_0);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Startall(void *argument_0, void *argument_1, void *result) {
  int native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = MPI_Startall(native_argument_0, (MPI_Request *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Sendrecv_replace(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *argument_7, void *argument_8, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  int native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_argument_6;
  memcpy(&native_argument_6, argument_6, sizeof(native_argument_6));
  MPI_Comm native_argument_7;
  memcpy(&native_argument_7, argument_7, sizeof(native_argument_7));
  int native_result = MPI_Sendrecv_replace((void *)argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5, native_argument_6, native_argument_7, (MPI_Status *)argument_8);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Type_contiguous(void *argument_0, void *argument_1, void *argument_2, void *result) {
  int native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  tMPI_Datatype native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = tMPI_Type_contiguous(native_argument_0, native_argument_1, (tMPI_Datatype *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Type_contiguous(void *argument_0, void *argument_1, void *argument_2, void *result) {
  int native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  MPI_Datatype native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = MPI_Type_contiguous(native_argument_0, native_argument_1, (MPI_Datatype *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Type_vector(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *result) {
  int native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  MPI_Datatype native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_result = MPI_Type_vector(native_argument_0, native_argument_1, native_argument_2, native_argument_3, (MPI_Datatype *)argument_4);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Type_indexed(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *result) {
  int native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  MPI_Datatype native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_result = MPI_Type_indexed(native_argument_0, (int *)argument_1, (int *)argument_2, native_argument_3, (MPI_Datatype *)argument_4);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Type_commit(void *argument_0, void *result) {
  int native_result = tMPI_Type_commit((tMPI_Datatype *)argument_0);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Type_commit(void *argument_0, void *result) {
  int native_result = MPI_Type_commit((MPI_Datatype *)argument_0);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Type_free(void *argument_0, void *result) {
  int native_result = tMPI_Type_free((tMPI_Datatype *)argument_0);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Type_free(void *argument_0, void *result) {
  int native_result = MPI_Type_free((MPI_Datatype *)argument_0);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Type_size(void *argument_0, void *argument_1, void *result) {
  MPI_Datatype native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = MPI_Type_size(native_argument_0, (int *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Get_elements(void *argument_0, void *argument_1, void *argument_2, void *result) {
  MPI_Datatype native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = MPI_Get_elements((MPI_Status *)argument_0, native_argument_1, (int *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Pack(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  MPI_Comm native_argument_6;
  memcpy(&native_argument_6, argument_6, sizeof(native_argument_6));
  int native_result = MPI_Pack((void *)argument_0, native_argument_1, native_argument_2, (void *)argument_3, native_argument_4, (int *)argument_5, native_argument_6);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Unpack(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  MPI_Datatype native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  MPI_Comm native_argument_6;
  memcpy(&native_argument_6, argument_6, sizeof(native_argument_6));
  int native_result = MPI_Unpack((void *)argument_0, native_argument_1, (int *)argument_2, (void *)argument_3, native_argument_4, native_argument_5, native_argument_6);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Pack_size(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *result) {
  int native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  MPI_Datatype native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Comm native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_result = MPI_Pack_size(native_argument_0, native_argument_1, native_argument_2, (int *)argument_3);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Gatherv(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *argument_7, void *argument_8, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  tMPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  tMPI_Datatype native_argument_6;
  memcpy(&native_argument_6, argument_6, sizeof(native_argument_6));
  int native_argument_7;
  memcpy(&native_argument_7, argument_7, sizeof(native_argument_7));
  tMPI_Comm native_argument_8;
  memcpy(&native_argument_8, argument_8, sizeof(native_argument_8));
  int native_result = tMPI_Gatherv((void *)argument_0, native_argument_1, native_argument_2, (void *)argument_3, (int *)argument_4, (int *)argument_5, native_argument_6, native_argument_7, native_argument_8);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Gatherv(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *argument_7, void *argument_8, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  MPI_Datatype native_argument_6;
  memcpy(&native_argument_6, argument_6, sizeof(native_argument_6));
  int native_argument_7;
  memcpy(&native_argument_7, argument_7, sizeof(native_argument_7));
  MPI_Comm native_argument_8;
  memcpy(&native_argument_8, argument_8, sizeof(native_argument_8));
  int native_result = MPI_Gatherv((void *)argument_0, native_argument_1, native_argument_2, (void *)argument_3, (int *)argument_4, (int *)argument_5, native_argument_6, native_argument_7, native_argument_8);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Scatterv(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *argument_7, void *argument_8, void *result) {
  tMPI_Datatype native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  tMPI_Datatype native_argument_6;
  memcpy(&native_argument_6, argument_6, sizeof(native_argument_6));
  int native_argument_7;
  memcpy(&native_argument_7, argument_7, sizeof(native_argument_7));
  tMPI_Comm native_argument_8;
  memcpy(&native_argument_8, argument_8, sizeof(native_argument_8));
  int native_result = tMPI_Scatterv((void *)argument_0, (int *)argument_1, (int *)argument_2, native_argument_3, (void *)argument_4, native_argument_5, native_argument_6, native_argument_7, native_argument_8);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Scatterv(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *argument_7, void *argument_8, void *result) {
  MPI_Datatype native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  MPI_Datatype native_argument_6;
  memcpy(&native_argument_6, argument_6, sizeof(native_argument_6));
  int native_argument_7;
  memcpy(&native_argument_7, argument_7, sizeof(native_argument_7));
  MPI_Comm native_argument_8;
  memcpy(&native_argument_8, argument_8, sizeof(native_argument_8));
  int native_result = MPI_Scatterv((void *)argument_0, (int *)argument_1, (int *)argument_2, native_argument_3, (void *)argument_4, native_argument_5, native_argument_6, native_argument_7, native_argument_8);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Alltoall(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  tMPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  tMPI_Datatype native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  tMPI_Comm native_argument_6;
  memcpy(&native_argument_6, argument_6, sizeof(native_argument_6));
  int native_result = tMPI_Alltoall((void *)argument_0, native_argument_1, native_argument_2, (void *)argument_3, native_argument_4, native_argument_5, native_argument_6);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Alltoall(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *result) {
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Datatype native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  MPI_Datatype native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  MPI_Comm native_argument_6;
  memcpy(&native_argument_6, argument_6, sizeof(native_argument_6));
  int native_result = MPI_Alltoall((void *)argument_0, native_argument_1, native_argument_2, (void *)argument_3, native_argument_4, native_argument_5, native_argument_6);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Alltoallv(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *argument_7, void *argument_8, void *result) {
  tMPI_Datatype native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  tMPI_Datatype native_argument_7;
  memcpy(&native_argument_7, argument_7, sizeof(native_argument_7));
  tMPI_Comm native_argument_8;
  memcpy(&native_argument_8, argument_8, sizeof(native_argument_8));
  int native_result = tMPI_Alltoallv((void *)argument_0, (int *)argument_1, (int *)argument_2, native_argument_3, (void *)argument_4, (int *)argument_5, (int *)argument_6, native_argument_7, native_argument_8);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Alltoallv(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *argument_7, void *argument_8, void *result) {
  MPI_Datatype native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  MPI_Datatype native_argument_7;
  memcpy(&native_argument_7, argument_7, sizeof(native_argument_7));
  MPI_Comm native_argument_8;
  memcpy(&native_argument_8, argument_8, sizeof(native_argument_8));
  int native_result = MPI_Alltoallv((void *)argument_0, (int *)argument_1, (int *)argument_2, native_argument_3, (void *)argument_4, (int *)argument_5, (int *)argument_6, native_argument_7, native_argument_8);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Reduce(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *result) {
  int native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  tMPI_Datatype native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  tMPI_Op native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  int native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  tMPI_Comm native_argument_6;
  memcpy(&native_argument_6, argument_6, sizeof(native_argument_6));
  int native_result = tMPI_Reduce((void *)argument_0, (void *)argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5, native_argument_6);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Reduce(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *argument_6, void *result) {
  int native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  MPI_Datatype native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  MPI_Op native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  int native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  MPI_Comm native_argument_6;
  memcpy(&native_argument_6, argument_6, sizeof(native_argument_6));
  int native_result = MPI_Reduce((void *)argument_0, (void *)argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5, native_argument_6);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Op_free(void *argument_0, void *result) {
  int native_result = MPI_Op_free((MPI_Op *)argument_0);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Reduce_scatter(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *result) {
  MPI_Datatype native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  MPI_Op native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  MPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = MPI_Reduce_scatter((void *)argument_0, (void *)argument_1, (int *)argument_2, native_argument_3, native_argument_4, native_argument_5);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Scan(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *result) {
  int native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  tMPI_Datatype native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  tMPI_Op native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  tMPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = tMPI_Scan((void *)argument_0, (void *)argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Scan(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *result) {
  int native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  MPI_Datatype native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  MPI_Op native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  MPI_Comm native_argument_5;
  memcpy(&native_argument_5, argument_5, sizeof(native_argument_5));
  int native_result = MPI_Scan((void *)argument_0, (void *)argument_1, native_argument_2, native_argument_3, native_argument_4, native_argument_5);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Group_size(void *argument_0, void *argument_1, void *result) {
  tMPI_Group native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = tMPI_Group_size(native_argument_0, (int *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Group_size(void *argument_0, void *argument_1, void *result) {
  MPI_Group native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = MPI_Group_size(native_argument_0, (int *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Group_rank(void *argument_0, void *argument_1, void *result) {
  tMPI_Group native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = tMPI_Group_rank(native_argument_0, (int *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Group_rank(void *argument_0, void *argument_1, void *result) {
  MPI_Group native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = MPI_Group_rank(native_argument_0, (int *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Group_translate_ranks(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *result) {
  MPI_Group native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Group native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_result = MPI_Group_translate_ranks(native_argument_0, native_argument_1, (int *)argument_2, native_argument_3, (int *)argument_4);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Group_compare(void *argument_0, void *argument_1, void *argument_2, void *result) {
  MPI_Group native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  MPI_Group native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = MPI_Group_compare(native_argument_0, native_argument_1, (int *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Comm_group(void *argument_0, void *argument_1, void *result) {
  tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = tMPI_Comm_group(native_argument_0, (tMPI_Group *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Comm_group(void *argument_0, void *argument_1, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = MPI_Comm_group(native_argument_0, (MPI_Group *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Group_union(void *argument_0, void *argument_1, void *argument_2, void *result) {
  MPI_Group native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  MPI_Group native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = MPI_Group_union(native_argument_0, native_argument_1, (MPI_Group *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Group_intersection(void *argument_0, void *argument_1, void *argument_2, void *result) {
  MPI_Group native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  MPI_Group native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = MPI_Group_intersection(native_argument_0, native_argument_1, (MPI_Group *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Group_difference(void *argument_0, void *argument_1, void *argument_2, void *result) {
  MPI_Group native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  MPI_Group native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = MPI_Group_difference(native_argument_0, native_argument_1, (MPI_Group *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Group_incl(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *result) {
  tMPI_Group native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = tMPI_Group_incl(native_argument_0, native_argument_1, (int *)argument_2, (tMPI_Group *)argument_3);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Group_incl(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *result) {
  MPI_Group native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = MPI_Group_incl(native_argument_0, native_argument_1, (int *)argument_2, (MPI_Group *)argument_3);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Group_excl(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *result) {
  MPI_Group native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = MPI_Group_excl(native_argument_0, native_argument_1, (int *)argument_2, (MPI_Group *)argument_3);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Group_range_incl(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *result) {
  MPI_Group native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = MPI_Group_range_incl(native_argument_0, native_argument_1, (int (*)[3])argument_2, (MPI_Group *)argument_3);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Group_range_excl(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *result) {
  MPI_Group native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = MPI_Group_range_excl(native_argument_0, native_argument_1, (int (*)[3])argument_2, (MPI_Group *)argument_3);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Group_free(void *argument_0, void *result) {
  int native_result = tMPI_Group_free((tMPI_Group *)argument_0);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Group_free(void *argument_0, void *result) {
  int native_result = MPI_Group_free((MPI_Group *)argument_0);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Comm_compare(void *argument_0, void *argument_1, void *argument_2, void *result) {
  tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  tMPI_Comm native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = tMPI_Comm_compare(native_argument_0, native_argument_1, (int *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Comm_compare(void *argument_0, void *argument_1, void *argument_2, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  MPI_Comm native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = MPI_Comm_compare(native_argument_0, native_argument_1, (int *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Comm_create(void *argument_0, void *argument_1, void *argument_2, void *result) {
  tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  tMPI_Group native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = tMPI_Comm_create(native_argument_0, native_argument_1, (tMPI_Comm *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Comm_create(void *argument_0, void *argument_1, void *argument_2, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  MPI_Group native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = MPI_Comm_create(native_argument_0, native_argument_1, (MPI_Comm *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Comm_test_inter(void *argument_0, void *argument_1, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = MPI_Comm_test_inter(native_argument_0, (int *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Comm_remote_size(void *argument_0, void *argument_1, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = MPI_Comm_remote_size(native_argument_0, (int *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Comm_remote_group(void *argument_0, void *argument_1, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = MPI_Comm_remote_group(native_argument_0, (MPI_Group *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Intercomm_create(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  MPI_Comm native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_argument_3;
  memcpy(&native_argument_3, argument_3, sizeof(native_argument_3));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  int native_result = MPI_Intercomm_create(native_argument_0, native_argument_1, native_argument_2, native_argument_3, native_argument_4, (MPI_Comm *)argument_5);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Intercomm_merge(void *argument_0, void *argument_1, void *argument_2, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = MPI_Intercomm_merge(native_argument_0, native_argument_1, (MPI_Comm *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Cart_create(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *result) {
  tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  int native_result = tMPI_Cart_create(native_argument_0, native_argument_1, (int *)argument_2, (int *)argument_3, native_argument_4, (tMPI_Comm *)argument_5);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Cart_create(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  int native_result = MPI_Cart_create(native_argument_0, native_argument_1, (int *)argument_2, (int *)argument_3, native_argument_4, (MPI_Comm *)argument_5);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Dims_create(void *argument_0, void *argument_1, void *argument_2, void *result) {
  int native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = MPI_Dims_create(native_argument_0, native_argument_1, (int *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Graph_create(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *argument_5, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_argument_4;
  memcpy(&native_argument_4, argument_4, sizeof(native_argument_4));
  int native_result = MPI_Graph_create(native_argument_0, native_argument_1, (int *)argument_2, (int *)argument_3, native_argument_4, (MPI_Comm *)argument_5);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Topo_test(void *argument_0, void *argument_1, void *result) {
  tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = tMPI_Topo_test(native_argument_0, (int *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Topo_test(void *argument_0, void *argument_1, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = MPI_Topo_test(native_argument_0, (int *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Graphdims_get(void *argument_0, void *argument_1, void *argument_2, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = MPI_Graphdims_get(native_argument_0, (int *)argument_1, (int *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Graph_get(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_result = MPI_Graph_get(native_argument_0, native_argument_1, native_argument_2, (int *)argument_3, (int *)argument_4);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Cartdim_get(void *argument_0, void *argument_1, void *result) {
  tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = tMPI_Cartdim_get(native_argument_0, (int *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Cartdim_get(void *argument_0, void *argument_1, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = MPI_Cartdim_get(native_argument_0, (int *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Cart_get(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *result) {
  tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = tMPI_Cart_get(native_argument_0, native_argument_1, (int *)argument_2, (int *)argument_3, (int *)argument_4);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Cart_get(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = MPI_Cart_get(native_argument_0, native_argument_1, (int *)argument_2, (int *)argument_3, (int *)argument_4);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Cart_rank(void *argument_0, void *argument_1, void *argument_2, void *result) {
  tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = tMPI_Cart_rank(native_argument_0, (int *)argument_1, (int *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Cart_rank(void *argument_0, void *argument_1, void *argument_2, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = MPI_Cart_rank(native_argument_0, (int *)argument_1, (int *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Cart_coords(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *result) {
  tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_result = tMPI_Cart_coords(native_argument_0, native_argument_1, native_argument_2, (int *)argument_3);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Cart_coords(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_result = MPI_Cart_coords(native_argument_0, native_argument_1, native_argument_2, (int *)argument_3);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Graph_neighbors_count(void *argument_0, void *argument_1, void *argument_2, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = MPI_Graph_neighbors_count(native_argument_0, native_argument_1, (int *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Graph_neighbors(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_result = MPI_Graph_neighbors(native_argument_0, native_argument_1, native_argument_2, (int *)argument_3);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Cart_shift(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_argument_2;
  memcpy(&native_argument_2, argument_2, sizeof(native_argument_2));
  int native_result = MPI_Cart_shift(native_argument_0, native_argument_1, native_argument_2, (int *)argument_3, (int *)argument_4);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Cart_sub(void *argument_0, void *argument_1, void *argument_2, void *result) {
  tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = tMPI_Cart_sub(native_argument_0, (int *)argument_1, (tMPI_Comm *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Cart_sub(void *argument_0, void *argument_1, void *argument_2, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = MPI_Cart_sub(native_argument_0, (int *)argument_1, (MPI_Comm *)argument_2);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Cart_map(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *result) {
  tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = tMPI_Cart_map(native_argument_0, native_argument_1, (int *)argument_2, (int *)argument_3, (int *)argument_4);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Cart_map(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = MPI_Cart_map(native_argument_0, native_argument_1, (int *)argument_2, (int *)argument_3, (int *)argument_4);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Graph_map(void *argument_0, void *argument_1, void *argument_2, void *argument_3, void *argument_4, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = MPI_Graph_map(native_argument_0, native_argument_1, (int *)argument_2, (int *)argument_3, (int *)argument_4);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Get_processor_name(void *argument_0, void *argument_1, void *result) {
  int native_result = tMPI_Get_processor_name((char *)argument_0, (int *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Get_processor_name(void *argument_0, void *argument_1, void *result) {
  int native_result = MPI_Get_processor_name((char *)argument_0, (int *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Errhandler_free(void *argument_0, void *result) {
  int native_result = MPI_Errhandler_free((MPI_Errhandler *)argument_0);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static void plumed_mpiw_MPI_Error_class(void *argument_0, void *argument_1, void *result) {
  int native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = MPI_Error_class(native_argument_0, (int *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Wtime(void *result) {
  double native_result = tMPI_Wtime();
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Wtime(void *result) {
  double native_result = MPI_Wtime();
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Wtick(void *result) {
  double native_result = tMPI_Wtick();
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Wtick(void *result) {
  double native_result = MPI_Wtick();
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 0)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Finalized(void *argument_0, void *result) {
  int native_result = tMPI_Finalized((int *)argument_0);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Finalized(void *argument_0, void *result) {
  int native_result = MPI_Finalized((int *)argument_0);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 0)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Comm_get_errhandler(void *argument_0, void *argument_1, void *result) {
  tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = tMPI_Comm_get_errhandler(native_argument_0, (tMPI_Errhandler *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Comm_get_errhandler(void *argument_0, void *argument_1, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  int native_result = MPI_Comm_get_errhandler(native_argument_0, (MPI_Errhandler *)argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif

#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 0)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_MPI_Comm_set_errhandler(void *argument_0, void *argument_1, void *result) {
  tMPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  tMPI_Errhandler native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = tMPI_Comm_set_errhandler(native_argument_0, native_argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#else
static void plumed_mpiw_MPI_Comm_set_errhandler(void *argument_0, void *argument_1, void *result) {
  MPI_Comm native_argument_0;
  memcpy(&native_argument_0, argument_0, sizeof(native_argument_0));
  MPI_Errhandler native_argument_1;
  memcpy(&native_argument_1, argument_1, sizeof(native_argument_1));
  int native_result = MPI_Comm_set_errhandler(native_argument_0, native_argument_1);
  memcpy(result, &native_result, sizeof(native_result));
}
#endif
#define PLUMED_MPIW_ENTRY(name, id, version, argc, fn) { (uint32_t)sizeof(plumed_mpiw_function_entry), #name, id, version, argc, (plumed_mpiw_generic_fn)(fn) }
/* Rows append contiguously in declaration order, including conditional rows. The identifier is metadata only; it is never an array index or byte offset. */
static const plumed_mpiw_function_entry plumed_mpiw_functions[] = {
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Initialized, 0xd4498369u, 1u, 1u, plumed_mpiw_MPI_Initialized),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Comm_dup, 0xf0dc6608u, 1u, 2u, plumed_mpiw_MPI_Comm_dup),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Comm_free, 0x4da4e1b2u, 1u, 1u, plumed_mpiw_MPI_Comm_free),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Comm_rank, 0x4e7672aau, 1u, 2u, plumed_mpiw_MPI_Comm_rank),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Comm_size, 0xda2369ecu, 1u, 2u, plumed_mpiw_MPI_Comm_size),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Comm_split, 0x97f1adf0u, 1u, 4u, plumed_mpiw_MPI_Comm_split),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Abort, 0x84f458d6u, 1u, 2u, plumed_mpiw_MPI_Abort),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Bcast, 0xe809e1cau, 1u, 5u, plumed_mpiw_MPI_Bcast),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Scatter, 0xb2fa0f58u, 1u, 8u, plumed_mpiw_MPI_Scatter),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Gather, 0xb1f903d9u, 1u, 8u, plumed_mpiw_MPI_Gather),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Allreduce, 0x788900dbu, 1u, 6u, plumed_mpiw_MPI_Allreduce),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Isend, 0xf5e9d915u, 1u, 7u, plumed_mpiw_MPI_Isend),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Wait, 0x4b15813du, 1u, 2u, plumed_mpiw_MPI_Wait),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Allgather, 0xab8b26c9u, 1u, 7u, plumed_mpiw_MPI_Allgather),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Allgatherv, 0x0f1b22f3u, 1u, 8u, plumed_mpiw_MPI_Allgatherv),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Recv, 0xdd56a4bcu, 1u, 7u, plumed_mpiw_MPI_Recv),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Barrier, 0xe7a052f4u, 1u, 1u, plumed_mpiw_MPI_Barrier),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Get_count, 0x71967a66u, 1u, 3u, plumed_mpiw_MPI_Get_count),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Sendrecv, 0xb09a6fc1u, 1u, 12u, plumed_mpiw_MPI_Sendrecv),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Error_string, 0xd272abecu, 1u, 3u, plumed_mpiw_MPI_Error_string),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 0)))
  PLUMED_MPIW_ENTRY(MPI_Comm_f2c, 0x9a9e7fcdu, 1u, 1u, plumed_mpiw_MPI_Comm_f2c),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Send, 0xbd2b6912u, 1u, 6u, plumed_mpiw_MPI_Send),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Bsend, 0xa9dbf685u, 1u, 6u, plumed_mpiw_MPI_Bsend),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Ssend, 0x5edcf0deu, 1u, 6u, plumed_mpiw_MPI_Ssend),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Rsend, 0xbfa09b06u, 1u, 6u, plumed_mpiw_MPI_Rsend),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Buffer_attach, 0xe00bd480u, 1u, 2u, plumed_mpiw_MPI_Buffer_attach),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Buffer_detach, 0x74dec40fu, 1u, 2u, plumed_mpiw_MPI_Buffer_detach),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Ibsend, 0xa4d683afu, 1u, 7u, plumed_mpiw_MPI_Ibsend),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Issend, 0x8275e866u, 1u, 7u, plumed_mpiw_MPI_Issend),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Irsend, 0xb13d1d70u, 1u, 7u, plumed_mpiw_MPI_Irsend),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Irecv, 0xb69a0083u, 1u, 7u, plumed_mpiw_MPI_Irecv),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Test, 0x1d207e4cu, 1u, 3u, plumed_mpiw_MPI_Test),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Request_free, 0x475f38a1u, 1u, 1u, plumed_mpiw_MPI_Request_free),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Waitany, 0x5acd4436u, 1u, 4u, plumed_mpiw_MPI_Waitany),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Testany, 0xf329c686u, 1u, 5u, plumed_mpiw_MPI_Testany),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Waitall, 0x7a57ae23u, 1u, 3u, plumed_mpiw_MPI_Waitall),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Testall, 0x593a1626u, 1u, 4u, plumed_mpiw_MPI_Testall),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Waitsome, 0xb6562e1au, 1u, 5u, plumed_mpiw_MPI_Waitsome),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Testsome, 0x5a87d17bu, 1u, 5u, plumed_mpiw_MPI_Testsome),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Iprobe, 0xd1ba8db8u, 1u, 5u, plumed_mpiw_MPI_Iprobe),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Probe, 0xfed81500u, 1u, 4u, plumed_mpiw_MPI_Probe),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Cancel, 0xed47f9abu, 1u, 1u, plumed_mpiw_MPI_Cancel),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Test_cancelled, 0xec5e3092u, 1u, 2u, plumed_mpiw_MPI_Test_cancelled),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Send_init, 0xd381e137u, 1u, 7u, plumed_mpiw_MPI_Send_init),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Bsend_init, 0x16ede22eu, 1u, 7u, plumed_mpiw_MPI_Bsend_init),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Ssend_init, 0xeea1c733u, 1u, 7u, plumed_mpiw_MPI_Ssend_init),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Rsend_init, 0x018fab06u, 1u, 7u, plumed_mpiw_MPI_Rsend_init),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Recv_init, 0x9cbb4e8au, 1u, 7u, plumed_mpiw_MPI_Recv_init),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Start, 0xfb62fbe3u, 1u, 1u, plumed_mpiw_MPI_Start),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Startall, 0xc3d7cddfu, 1u, 2u, plumed_mpiw_MPI_Startall),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Sendrecv_replace, 0x7c6ffc57u, 1u, 9u, plumed_mpiw_MPI_Sendrecv_replace),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Type_contiguous, 0xe41681d2u, 1u, 3u, plumed_mpiw_MPI_Type_contiguous),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Type_vector, 0xc7b9dbb0u, 1u, 5u, plumed_mpiw_MPI_Type_vector),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Type_indexed, 0xeb606d34u, 1u, 5u, plumed_mpiw_MPI_Type_indexed),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Type_commit, 0xf848da1fu, 1u, 1u, plumed_mpiw_MPI_Type_commit),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Type_free, 0xdacde819u, 1u, 1u, plumed_mpiw_MPI_Type_free),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Type_size, 0x95660c77u, 1u, 2u, plumed_mpiw_MPI_Type_size),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Get_elements, 0xf71eaa9au, 1u, 3u, plumed_mpiw_MPI_Get_elements),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Pack, 0x5649ba7au, 1u, 7u, plumed_mpiw_MPI_Pack),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Unpack, 0x49901f8bu, 1u, 7u, plumed_mpiw_MPI_Unpack),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Pack_size, 0x24de0131u, 1u, 4u, plumed_mpiw_MPI_Pack_size),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Gatherv, 0xc90455c9u, 1u, 9u, plumed_mpiw_MPI_Gatherv),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Scatterv, 0x765eb673u, 1u, 9u, plumed_mpiw_MPI_Scatterv),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Alltoall, 0xcfe270aeu, 1u, 7u, plumed_mpiw_MPI_Alltoall),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Alltoallv, 0xa4d99ed0u, 1u, 9u, plumed_mpiw_MPI_Alltoallv),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Reduce, 0x3a815604u, 1u, 7u, plumed_mpiw_MPI_Reduce),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Op_free, 0x94fea4e6u, 1u, 1u, plumed_mpiw_MPI_Op_free),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Reduce_scatter, 0x2f70a384u, 1u, 6u, plumed_mpiw_MPI_Reduce_scatter),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Scan, 0x64e2ee54u, 1u, 6u, plumed_mpiw_MPI_Scan),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Group_size, 0x7de536d5u, 1u, 2u, plumed_mpiw_MPI_Group_size),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Group_rank, 0x668fab96u, 1u, 2u, plumed_mpiw_MPI_Group_rank),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Group_translate_ranks, 0x2ac4d2c6u, 1u, 5u, plumed_mpiw_MPI_Group_translate_ranks),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Group_compare, 0xb0d34b44u, 1u, 3u, plumed_mpiw_MPI_Group_compare),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Comm_group, 0xd11f4483u, 1u, 2u, plumed_mpiw_MPI_Comm_group),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Group_union, 0xc288bd65u, 1u, 3u, plumed_mpiw_MPI_Group_union),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Group_intersection, 0xddf0ecb4u, 1u, 3u, plumed_mpiw_MPI_Group_intersection),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Group_difference, 0x37101ec1u, 1u, 3u, plumed_mpiw_MPI_Group_difference),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Group_incl, 0xbe8bfe71u, 1u, 4u, plumed_mpiw_MPI_Group_incl),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Group_excl, 0x45fa0902u, 1u, 4u, plumed_mpiw_MPI_Group_excl),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Group_range_incl, 0x801adc54u, 1u, 4u, plumed_mpiw_MPI_Group_range_incl),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Group_range_excl, 0xc416deb6u, 1u, 4u, plumed_mpiw_MPI_Group_range_excl),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Group_free, 0xfd1c62d1u, 1u, 1u, plumed_mpiw_MPI_Group_free),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Comm_compare, 0xdacd79cdu, 1u, 3u, plumed_mpiw_MPI_Comm_compare),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Comm_create, 0x1de85fa0u, 1u, 3u, plumed_mpiw_MPI_Comm_create),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Comm_test_inter, 0x29431089u, 1u, 2u, plumed_mpiw_MPI_Comm_test_inter),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Comm_remote_size, 0xaa0c6ac3u, 1u, 2u, plumed_mpiw_MPI_Comm_remote_size),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Comm_remote_group, 0x49ee9d8fu, 1u, 2u, plumed_mpiw_MPI_Comm_remote_group),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Intercomm_create, 0x320542c0u, 1u, 6u, plumed_mpiw_MPI_Intercomm_create),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Intercomm_merge, 0x1f6a920au, 1u, 3u, plumed_mpiw_MPI_Intercomm_merge),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Cart_create, 0xfd363844u, 1u, 6u, plumed_mpiw_MPI_Cart_create),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Dims_create, 0xfb2dd1cdu, 1u, 3u, plumed_mpiw_MPI_Dims_create),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Graph_create, 0x00eb0f02u, 1u, 6u, plumed_mpiw_MPI_Graph_create),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Topo_test, 0x09cbe810u, 1u, 2u, plumed_mpiw_MPI_Topo_test),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Graphdims_get, 0xc56e163bu, 1u, 3u, plumed_mpiw_MPI_Graphdims_get),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Graph_get, 0x12829ea5u, 1u, 5u, plumed_mpiw_MPI_Graph_get),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Cartdim_get, 0xa979e88cu, 1u, 2u, plumed_mpiw_MPI_Cartdim_get),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Cart_get, 0x0a3d9750u, 1u, 5u, plumed_mpiw_MPI_Cart_get),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Cart_rank, 0xf899fe26u, 1u, 3u, plumed_mpiw_MPI_Cart_rank),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Cart_coords, 0xf0bdc30cu, 1u, 4u, plumed_mpiw_MPI_Cart_coords),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Graph_neighbors_count, 0x5a79d6b6u, 1u, 3u, plumed_mpiw_MPI_Graph_neighbors_count),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Graph_neighbors, 0x567d4ff6u, 1u, 4u, plumed_mpiw_MPI_Graph_neighbors),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Cart_shift, 0x4d227e7cu, 1u, 5u, plumed_mpiw_MPI_Cart_shift),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Cart_sub, 0x1e9d3861u, 1u, 3u, plumed_mpiw_MPI_Cart_sub),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Cart_map, 0xba5a7b4cu, 1u, 5u, plumed_mpiw_MPI_Cart_map),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Graph_map, 0xefd795ebu, 1u, 5u, plumed_mpiw_MPI_Graph_map),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Get_processor_name, 0xf8e4ccc3u, 1u, 2u, plumed_mpiw_MPI_Get_processor_name),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Errhandler_free, 0xf03f8309u, 1u, 1u, plumed_mpiw_MPI_Errhandler_free),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_ENTRY(MPI_Error_class, 0x3b9188a6u, 1u, 2u, plumed_mpiw_MPI_Error_class),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Wtime, 0x75704e2eu, 1u, 0u, plumed_mpiw_MPI_Wtime),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Wtick, 0x84f54b80u, 1u, 0u, plumed_mpiw_MPI_Wtick),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 0)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Finalized, 0x105b65c7u, 1u, 1u, plumed_mpiw_MPI_Finalized),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 0)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Comm_get_errhandler, 0x74d05159u, 1u, 2u, plumed_mpiw_MPI_Comm_get_errhandler),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 0)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
  PLUMED_MPIW_ENTRY(MPI_Comm_set_errhandler, 0xe6af3a8fu, 1u, 2u, plumed_mpiw_MPI_Comm_set_errhandler),
#endif
};
#undef PLUMED_MPIW_ENTRY

#define PLUMED_MPIW_HANDLE_CONSTANT(name, type) \
  static void plumed_mpiw_value_##name(void *out) { \
    type value = name; \
    memcpy(out, &value, sizeof(value)); \
  }
#define PLUMED_MPIW_INTEGER_CONSTANT(name) \
  static void plumed_mpiw_value_##name(void *out) { \
    int value = (int)(name); \
    memcpy(out, &value, sizeof(value)); \
  }
#define PLUMED_MPIW_POINTER_CONSTANT(name) \
  static void plumed_mpiw_value_##name(void *out) { \
    void *value = (void *)(name); \
    memcpy(out, &value, sizeof(value)); \
  }
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_FLOAT(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Datatype value = TMPI_FLOAT;
#else
  const MPI_Datatype value = MPI_FLOAT;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_DOUBLE(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Datatype value = TMPI_DOUBLE;
#else
  const MPI_Datatype value = MPI_DOUBLE;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_INT(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Datatype value = TMPI_INT;
#else
  const MPI_Datatype value = MPI_INT;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_CHAR(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Datatype value = TMPI_CHAR;
#else
  const MPI_Datatype value = MPI_CHAR;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_UNSIGNED(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Datatype value = TMPI_UNSIGNED;
#else
  const MPI_Datatype value = MPI_UNSIGNED;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_UNSIGNED_LONG(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Datatype value = TMPI_UNSIGNED_LONG;
#else
  const MPI_Datatype value = MPI_UNSIGNED_LONG;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_UNSIGNED_LONG_LONG(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Datatype value = TMPI_UNSIGNED_LONG_LONG;
#else
  const MPI_Datatype value = MPI_UNSIGNED_LONG_LONG;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_LONG_DOUBLE(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Datatype value = TMPI_LONG_DOUBLE;
#else
  const MPI_Datatype value = MPI_LONG_DOUBLE;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_COMM_WORLD(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Comm value = TMPI_COMM_WORLD;
#else
  const MPI_Comm value = MPI_COMM_WORLD;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_COMM_SELF(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Comm value = tMPI_Get_comm_self();
#else
  const MPI_Comm value = MPI_COMM_SELF;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_COMM_NULL(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Comm value = (tMPI_Comm)0;
#else
  const MPI_Comm value = MPI_COMM_NULL;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 0)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_IN_PLACE(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const void * value = TMPI_IN_PLACE;
#else
  const void * value = MPI_IN_PLACE;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_SUM(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Op value = TMPI_SUM;
#else
  const MPI_Op value = MPI_SUM;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_PROD(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Op value = TMPI_PROD;
#else
  const MPI_Op value = MPI_PROD;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_MAX(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Op value = TMPI_MAX;
#else
  const MPI_Op value = MPI_MAX;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_MIN(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Op value = TMPI_MIN;
#else
  const MPI_Op value = MPI_MIN;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_STATUS_IGNORE(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Status * value = TMPI_STATUS_IGNORE;
#else
  const MPI_Status * value = MPI_STATUS_IGNORE;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_STATUSES_IGNORE(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Status * value = TMPI_STATUSES_IGNORE;
#else
  const MPI_Status * value = MPI_STATUSES_IGNORE;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_GROUP_NULL(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Group value = (tMPI_Group)0;
#else
  const MPI_Group value = MPI_GROUP_NULL;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_GROUP_EMPTY(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Group value = TMPI_GROUP_EMPTY;
#else
  const MPI_Group value = MPI_GROUP_EMPTY;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_REQUEST_NULL(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Request value = TMPI_REQUEST_NULL;
#else
  const MPI_Request value = MPI_REQUEST_NULL;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
PLUMED_MPIW_HANDLE_CONSTANT(MPI_DATATYPE_NULL, MPI_Datatype)
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
PLUMED_MPIW_HANDLE_CONSTANT(MPI_OP_NULL, MPI_Op)
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
PLUMED_MPIW_HANDLE_CONSTANT(MPI_ERRHANDLER_NULL, MPI_Errhandler)
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_ERRORS_ARE_FATAL(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Errhandler value = TMPI_ERRORS_ARE_FATAL;
#else
  const MPI_Errhandler value = MPI_ERRORS_ARE_FATAL;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_ERRORS_RETURN(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Errhandler value = TMPI_ERRORS_RETURN;
#else
  const MPI_Errhandler value = MPI_ERRORS_RETURN;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_BYTE(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Datatype value = TMPI_BYTE;
#else
  const MPI_Datatype value = MPI_BYTE;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
PLUMED_MPIW_HANDLE_CONSTANT(MPI_PACKED, MPI_Datatype)
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_SHORT(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Datatype value = TMPI_SHORT;
#else
  const MPI_Datatype value = MPI_SHORT;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_LONG(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Datatype value = TMPI_LONG;
#else
  const MPI_Datatype value = MPI_LONG;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_UNSIGNED_CHAR(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Datatype value = TMPI_UNSIGNED_CHAR;
#else
  const MPI_Datatype value = MPI_UNSIGNED_CHAR;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_UNSIGNED_SHORT(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Datatype value = TMPI_UNSIGNED_SHORT;
#else
  const MPI_Datatype value = MPI_UNSIGNED_SHORT;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_LONG_LONG_INT(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Datatype value = TMPI_LONG_LONG_INT;
#else
  const MPI_Datatype value = MPI_LONG_LONG_INT;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_LONG_LONG(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Datatype value = TMPI_LONG_LONG;
#else
  const MPI_Datatype value = MPI_LONG_LONG;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
PLUMED_MPIW_HANDLE_CONSTANT(MPI_FLOAT_INT, MPI_Datatype)
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
PLUMED_MPIW_HANDLE_CONSTANT(MPI_DOUBLE_INT, MPI_Datatype)
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
PLUMED_MPIW_HANDLE_CONSTANT(MPI_LONG_INT, MPI_Datatype)
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
PLUMED_MPIW_HANDLE_CONSTANT(MPI_2INT, MPI_Datatype)
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
PLUMED_MPIW_HANDLE_CONSTANT(MPI_SHORT_INT, MPI_Datatype)
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
PLUMED_MPIW_HANDLE_CONSTANT(MPI_LONG_DOUBLE_INT, MPI_Datatype)
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
PLUMED_MPIW_HANDLE_CONSTANT(MPI_MAXLOC, MPI_Op)
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
PLUMED_MPIW_HANDLE_CONSTANT(MPI_MINLOC, MPI_Op)
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_LAND(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Op value = TMPI_LAND;
#else
  const MPI_Op value = MPI_LAND;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_BAND(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Op value = TMPI_BAND;
#else
  const MPI_Op value = MPI_BAND;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_LOR(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Op value = TMPI_LOR;
#else
  const MPI_Op value = MPI_LOR;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_BOR(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Op value = TMPI_BOR;
#else
  const MPI_Op value = MPI_BOR;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_LXOR(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Op value = TMPI_LXOR;
#else
  const MPI_Op value = MPI_LXOR;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_BXOR(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Op value = TMPI_BXOR;
#else
  const MPI_Op value = MPI_BXOR;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_UNDEFINED(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const int value = TMPI_UNDEFINED;
#else
  const int32_t value = MPI_UNDEFINED;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_ANY_SOURCE(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const int value = TMPI_ANY_SOURCE;
#else
  const int32_t value = MPI_ANY_SOURCE;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_ANY_TAG(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const int value = TMPI_ANY_TAG;
#else
  const int32_t value = MPI_ANY_TAG;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
PLUMED_MPIW_INTEGER_CONSTANT(MPI_PROC_NULL)
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_IDENT(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const int value = TMPI_IDENT;
#else
  const int32_t value = MPI_IDENT;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_CONGRUENT(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const int value = TMPI_CONGRUENT;
#else
  const int32_t value = MPI_CONGRUENT;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_SIMILAR(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const int value = TMPI_SIMILAR;
#else
  const int32_t value = MPI_SIMILAR;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_UNEQUAL(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const int value = TMPI_UNEQUAL;
#else
  const int32_t value = MPI_UNEQUAL;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_CART(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const int value = TMPI_CART;
#else
  const int32_t value = MPI_CART;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_GRAPH(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const int value = TMPI_GRAPH;
#else
  const int32_t value = MPI_GRAPH;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_MAX_PROCESSOR_NAME(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const int value = TMPI_MAX_PROCESSOR_NAME;
#else
  const int32_t value = MPI_MAX_PROCESSOR_NAME;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_MAX_ERROR_STRING(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const int value = TMPI_MAX_ERROR_STRING;
#else
  const int32_t value = MPI_MAX_ERROR_STRING;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
PLUMED_MPIW_INTEGER_CONSTANT(MPI_BSEND_OVERHEAD)
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_SUCCESS(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const int value = TMPI_SUCCESS;
#else
  const int32_t value = MPI_SUCCESS;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
PLUMED_MPIW_INTEGER_CONSTANT(MPI_ERR_BUFFER)
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
PLUMED_MPIW_INTEGER_CONSTANT(MPI_ERR_COUNT)
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_ERR_TYPE(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const int value = TMPI_ERR_TYPE;
#else
  const int32_t value = MPI_ERR_TYPE;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
PLUMED_MPIW_INTEGER_CONSTANT(MPI_ERR_TAG)
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_ERR_COMM(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const int value = TMPI_ERR_COMM;
#else
  const int32_t value = MPI_ERR_COMM;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
PLUMED_MPIW_INTEGER_CONSTANT(MPI_ERR_RANK)
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
PLUMED_MPIW_INTEGER_CONSTANT(MPI_ERR_REQUEST)
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
PLUMED_MPIW_INTEGER_CONSTANT(MPI_ERR_ROOT)
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_ERR_GROUP(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const int value = TMPI_ERR_GROUP;
#else
  const int32_t value = MPI_ERR_GROUP;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
PLUMED_MPIW_INTEGER_CONSTANT(MPI_ERR_OP)
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
PLUMED_MPIW_INTEGER_CONSTANT(MPI_ERR_TOPOLOGY)
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_ERR_DIMS(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const int value = TMPI_ERR_DIMS;
#else
  const int32_t value = MPI_ERR_DIMS;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
PLUMED_MPIW_INTEGER_CONSTANT(MPI_ERR_ARG)
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_ERR_UNKNOWN(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const int value = TMPI_ERR_UNKNOWN;
#else
  const int32_t value = MPI_ERR_UNKNOWN;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
PLUMED_MPIW_INTEGER_CONSTANT(MPI_ERR_TRUNCATE)
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
PLUMED_MPIW_INTEGER_CONSTANT(MPI_ERR_OTHER)
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
PLUMED_MPIW_INTEGER_CONSTANT(MPI_ERR_INTERN)
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
PLUMED_MPIW_INTEGER_CONSTANT(MPI_ERR_PENDING)
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
PLUMED_MPIW_POINTER_CONSTANT(MPI_BOTTOM)
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_ERR_IN_STATUS(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const int value = TMPI_ERR_IN_STATUS;
#else
  const int32_t value = MPI_ERR_IN_STATUS;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_SIGNED_CHAR(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Datatype value = TMPI_SIGNED_CHAR;
#else
  const MPI_Datatype value = MPI_SIGNED_CHAR;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 2)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_C_BOOL(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Datatype value = TMPI_C_BOOL;
#else
  const MPI_Datatype value = MPI_C_BOOL;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 2)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static void plumed_mpiw_value_MPI_INT64_T(void *out) {
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  const tMPI_Datatype value = TMPI_INT64_T;
#else
  const MPI_Datatype value = MPI_INT64_T;
#endif
  memcpy(out, &value, sizeof(value));
}
#endif
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
static plumed_mpiw_generic_fn plumed_mpiw_native_runtime_identity = (plumed_mpiw_generic_fn)&tMPI_Initialized;
#else
static plumed_mpiw_generic_fn plumed_mpiw_native_runtime_identity = (plumed_mpiw_generic_fn)&MPI_Initialized;
#endif
#undef PLUMED_MPIW_HANDLE_CONSTANT
#undef PLUMED_MPIW_INTEGER_CONSTANT
#undef PLUMED_MPIW_POINTER_CONSTANT

#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_FLOAT = &plumed_mpiw_value_MPI_FLOAT;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_DOUBLE = &plumed_mpiw_value_MPI_DOUBLE;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_INT = &plumed_mpiw_value_MPI_INT;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_CHAR = &plumed_mpiw_value_MPI_CHAR;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_UNSIGNED = &plumed_mpiw_value_MPI_UNSIGNED;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_UNSIGNED_LONG = &plumed_mpiw_value_MPI_UNSIGNED_LONG;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_UNSIGNED_LONG_LONG = &plumed_mpiw_value_MPI_UNSIGNED_LONG_LONG;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_LONG_DOUBLE = &plumed_mpiw_value_MPI_LONG_DOUBLE;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_COMM_WORLD = &plumed_mpiw_value_MPI_COMM_WORLD;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_COMM_SELF = &plumed_mpiw_value_MPI_COMM_SELF;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_COMM_NULL = &plumed_mpiw_value_MPI_COMM_NULL;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 0)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_IN_PLACE = &plumed_mpiw_value_MPI_IN_PLACE;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_SUM = &plumed_mpiw_value_MPI_SUM;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_PROD = &plumed_mpiw_value_MPI_PROD;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_MAX = &plumed_mpiw_value_MPI_MAX;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_MIN = &plumed_mpiw_value_MPI_MIN;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_STATUS_IGNORE = &plumed_mpiw_value_MPI_STATUS_IGNORE;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_STATUSES_IGNORE = &plumed_mpiw_value_MPI_STATUSES_IGNORE;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_GROUP_NULL = &plumed_mpiw_value_MPI_GROUP_NULL;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_GROUP_EMPTY = &plumed_mpiw_value_MPI_GROUP_EMPTY;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_REQUEST_NULL = &plumed_mpiw_value_MPI_REQUEST_NULL;
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_DATATYPE_NULL = &plumed_mpiw_value_MPI_DATATYPE_NULL;
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_OP_NULL = &plumed_mpiw_value_MPI_OP_NULL;
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_ERRHANDLER_NULL = &plumed_mpiw_value_MPI_ERRHANDLER_NULL;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_ERRORS_ARE_FATAL = &plumed_mpiw_value_MPI_ERRORS_ARE_FATAL;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_ERRORS_RETURN = &plumed_mpiw_value_MPI_ERRORS_RETURN;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_BYTE = &plumed_mpiw_value_MPI_BYTE;
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_PACKED = &plumed_mpiw_value_MPI_PACKED;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_SHORT = &plumed_mpiw_value_MPI_SHORT;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_LONG = &plumed_mpiw_value_MPI_LONG;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_UNSIGNED_CHAR = &plumed_mpiw_value_MPI_UNSIGNED_CHAR;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_UNSIGNED_SHORT = &plumed_mpiw_value_MPI_UNSIGNED_SHORT;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_LONG_LONG_INT = &plumed_mpiw_value_MPI_LONG_LONG_INT;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_LONG_LONG = &plumed_mpiw_value_MPI_LONG_LONG;
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_FLOAT_INT = &plumed_mpiw_value_MPI_FLOAT_INT;
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_DOUBLE_INT = &plumed_mpiw_value_MPI_DOUBLE_INT;
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_LONG_INT = &plumed_mpiw_value_MPI_LONG_INT;
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_2INT = &plumed_mpiw_value_MPI_2INT;
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_SHORT_INT = &plumed_mpiw_value_MPI_SHORT_INT;
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_LONG_DOUBLE_INT = &plumed_mpiw_value_MPI_LONG_DOUBLE_INT;
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_MAXLOC = &plumed_mpiw_value_MPI_MAXLOC;
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_MINLOC = &plumed_mpiw_value_MPI_MINLOC;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_LAND = &plumed_mpiw_value_MPI_LAND;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_BAND = &plumed_mpiw_value_MPI_BAND;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_LOR = &plumed_mpiw_value_MPI_LOR;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_BOR = &plumed_mpiw_value_MPI_BOR;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_LXOR = &plumed_mpiw_value_MPI_LXOR;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_BXOR = &plumed_mpiw_value_MPI_BXOR;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_UNDEFINED = &plumed_mpiw_value_MPI_UNDEFINED;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_ANY_SOURCE = &plumed_mpiw_value_MPI_ANY_SOURCE;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_ANY_TAG = &plumed_mpiw_value_MPI_ANY_TAG;
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_PROC_NULL = &plumed_mpiw_value_MPI_PROC_NULL;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_IDENT = &plumed_mpiw_value_MPI_IDENT;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_CONGRUENT = &plumed_mpiw_value_MPI_CONGRUENT;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_SIMILAR = &plumed_mpiw_value_MPI_SIMILAR;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_UNEQUAL = &plumed_mpiw_value_MPI_UNEQUAL;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_CART = &plumed_mpiw_value_MPI_CART;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_GRAPH = &plumed_mpiw_value_MPI_GRAPH;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_MAX_PROCESSOR_NAME = &plumed_mpiw_value_MPI_MAX_PROCESSOR_NAME;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_MAX_ERROR_STRING = &plumed_mpiw_value_MPI_MAX_ERROR_STRING;
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_BSEND_OVERHEAD = &plumed_mpiw_value_MPI_BSEND_OVERHEAD;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_SUCCESS = &plumed_mpiw_value_MPI_SUCCESS;
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_ERR_BUFFER = &plumed_mpiw_value_MPI_ERR_BUFFER;
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_ERR_COUNT = &plumed_mpiw_value_MPI_ERR_COUNT;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_ERR_TYPE = &plumed_mpiw_value_MPI_ERR_TYPE;
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_ERR_TAG = &plumed_mpiw_value_MPI_ERR_TAG;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_ERR_COMM = &plumed_mpiw_value_MPI_ERR_COMM;
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_ERR_RANK = &plumed_mpiw_value_MPI_ERR_RANK;
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_ERR_REQUEST = &plumed_mpiw_value_MPI_ERR_REQUEST;
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_ERR_ROOT = &plumed_mpiw_value_MPI_ERR_ROOT;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_ERR_GROUP = &plumed_mpiw_value_MPI_ERR_GROUP;
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_ERR_OP = &plumed_mpiw_value_MPI_ERR_OP;
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_ERR_TOPOLOGY = &plumed_mpiw_value_MPI_ERR_TOPOLOGY;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_ERR_DIMS = &plumed_mpiw_value_MPI_ERR_DIMS;
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_ERR_ARG = &plumed_mpiw_value_MPI_ERR_ARG;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_ERR_UNKNOWN = &plumed_mpiw_value_MPI_ERR_UNKNOWN;
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_ERR_TRUNCATE = &plumed_mpiw_value_MPI_ERR_TRUNCATE;
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_ERR_OTHER = &plumed_mpiw_value_MPI_ERR_OTHER;
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_ERR_INTERN = &plumed_mpiw_value_MPI_ERR_INTERN;
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_ERR_PENDING = &plumed_mpiw_value_MPI_ERR_PENDING;
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_BOTTOM = &plumed_mpiw_value_MPI_BOTTOM;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_ERR_IN_STATUS = &plumed_mpiw_value_MPI_ERR_IN_STATUS;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_SIGNED_CHAR = &plumed_mpiw_value_MPI_SIGNED_CHAR;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 2)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_C_BOOL = &plumed_mpiw_value_MPI_C_BOOL;
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 2)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
static const plumed_mpiw_constant_accessor_fn plumed_mpiw_accessor_MPI_INT64_T = &plumed_mpiw_value_MPI_INT64_T;
#endif
#define PLUMED_MPIW_VALUE(name, type) \
  { (uint32_t)sizeof(plumed_mpiw_metadata_entry), PLUMED_MPIW_VALUE_ACCESSOR, \
    0u, "value:" #name, sizeof(plumed_mpiw_constant_accessor_fn), sizeof(type), \
    (const void *)&plumed_mpiw_accessor_##name }
static const plumed_mpiw_metadata_entry plumed_mpiw_metadata[] = {
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_FLOAT, tMPI_Datatype),
#else
  PLUMED_MPIW_VALUE(MPI_FLOAT, MPI_Datatype),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_DOUBLE, tMPI_Datatype),
#else
  PLUMED_MPIW_VALUE(MPI_DOUBLE, MPI_Datatype),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_INT, tMPI_Datatype),
#else
  PLUMED_MPIW_VALUE(MPI_INT, MPI_Datatype),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_CHAR, tMPI_Datatype),
#else
  PLUMED_MPIW_VALUE(MPI_CHAR, MPI_Datatype),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_UNSIGNED, tMPI_Datatype),
#else
  PLUMED_MPIW_VALUE(MPI_UNSIGNED, MPI_Datatype),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_UNSIGNED_LONG, tMPI_Datatype),
#else
  PLUMED_MPIW_VALUE(MPI_UNSIGNED_LONG, MPI_Datatype),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_UNSIGNED_LONG_LONG, tMPI_Datatype),
#else
  PLUMED_MPIW_VALUE(MPI_UNSIGNED_LONG_LONG, MPI_Datatype),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_LONG_DOUBLE, tMPI_Datatype),
#else
  PLUMED_MPIW_VALUE(MPI_LONG_DOUBLE, MPI_Datatype),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_COMM_WORLD, tMPI_Comm),
#else
  PLUMED_MPIW_VALUE(MPI_COMM_WORLD, MPI_Comm),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_COMM_SELF, tMPI_Comm),
#else
  PLUMED_MPIW_VALUE(MPI_COMM_SELF, MPI_Comm),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_COMM_NULL, tMPI_Comm),
#else
  PLUMED_MPIW_VALUE(MPI_COMM_NULL, MPI_Comm),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 0)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_IN_PLACE, void *),
#else
  PLUMED_MPIW_VALUE(MPI_IN_PLACE, void *),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_SUM, tMPI_Op),
#else
  PLUMED_MPIW_VALUE(MPI_SUM, MPI_Op),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_PROD, tMPI_Op),
#else
  PLUMED_MPIW_VALUE(MPI_PROD, MPI_Op),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_MAX, tMPI_Op),
#else
  PLUMED_MPIW_VALUE(MPI_MAX, MPI_Op),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_MIN, tMPI_Op),
#else
  PLUMED_MPIW_VALUE(MPI_MIN, MPI_Op),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_STATUS_IGNORE, tMPI_Status *),
#else
  PLUMED_MPIW_VALUE(MPI_STATUS_IGNORE, MPI_Status *),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_STATUSES_IGNORE, tMPI_Status *),
#else
  PLUMED_MPIW_VALUE(MPI_STATUSES_IGNORE, MPI_Status *),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_GROUP_NULL, tMPI_Group),
#else
  PLUMED_MPIW_VALUE(MPI_GROUP_NULL, MPI_Group),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_GROUP_EMPTY, tMPI_Group),
#else
  PLUMED_MPIW_VALUE(MPI_GROUP_EMPTY, MPI_Group),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_REQUEST_NULL, tMPI_Request),
#else
  PLUMED_MPIW_VALUE(MPI_REQUEST_NULL, MPI_Request),
#endif
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_VALUE(MPI_DATATYPE_NULL, MPI_Datatype),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_VALUE(MPI_OP_NULL, MPI_Op),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_VALUE(MPI_ERRHANDLER_NULL, MPI_Errhandler),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_ERRORS_ARE_FATAL, tMPI_Errhandler),
#else
  PLUMED_MPIW_VALUE(MPI_ERRORS_ARE_FATAL, MPI_Errhandler),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_ERRORS_RETURN, tMPI_Errhandler),
#else
  PLUMED_MPIW_VALUE(MPI_ERRORS_RETURN, MPI_Errhandler),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_BYTE, tMPI_Datatype),
#else
  PLUMED_MPIW_VALUE(MPI_BYTE, MPI_Datatype),
#endif
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_VALUE(MPI_PACKED, MPI_Datatype),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_SHORT, tMPI_Datatype),
#else
  PLUMED_MPIW_VALUE(MPI_SHORT, MPI_Datatype),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_LONG, tMPI_Datatype),
#else
  PLUMED_MPIW_VALUE(MPI_LONG, MPI_Datatype),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_UNSIGNED_CHAR, tMPI_Datatype),
#else
  PLUMED_MPIW_VALUE(MPI_UNSIGNED_CHAR, MPI_Datatype),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_UNSIGNED_SHORT, tMPI_Datatype),
#else
  PLUMED_MPIW_VALUE(MPI_UNSIGNED_SHORT, MPI_Datatype),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_LONG_LONG_INT, tMPI_Datatype),
#else
  PLUMED_MPIW_VALUE(MPI_LONG_LONG_INT, MPI_Datatype),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_LONG_LONG, tMPI_Datatype),
#else
  PLUMED_MPIW_VALUE(MPI_LONG_LONG, MPI_Datatype),
#endif
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_VALUE(MPI_FLOAT_INT, MPI_Datatype),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_VALUE(MPI_DOUBLE_INT, MPI_Datatype),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_VALUE(MPI_LONG_INT, MPI_Datatype),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_VALUE(MPI_2INT, MPI_Datatype),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_VALUE(MPI_SHORT_INT, MPI_Datatype),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_VALUE(MPI_LONG_DOUBLE_INT, MPI_Datatype),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_VALUE(MPI_MAXLOC, MPI_Op),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_VALUE(MPI_MINLOC, MPI_Op),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_LAND, tMPI_Op),
#else
  PLUMED_MPIW_VALUE(MPI_LAND, MPI_Op),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_BAND, tMPI_Op),
#else
  PLUMED_MPIW_VALUE(MPI_BAND, MPI_Op),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_LOR, tMPI_Op),
#else
  PLUMED_MPIW_VALUE(MPI_LOR, MPI_Op),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_BOR, tMPI_Op),
#else
  PLUMED_MPIW_VALUE(MPI_BOR, MPI_Op),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_LXOR, tMPI_Op),
#else
  PLUMED_MPIW_VALUE(MPI_LXOR, MPI_Op),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_BXOR, tMPI_Op),
#else
  PLUMED_MPIW_VALUE(MPI_BXOR, MPI_Op),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_UNDEFINED, int32_t),
#else
  PLUMED_MPIW_VALUE(MPI_UNDEFINED, int32_t),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_ANY_SOURCE, int32_t),
#else
  PLUMED_MPIW_VALUE(MPI_ANY_SOURCE, int32_t),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_ANY_TAG, int32_t),
#else
  PLUMED_MPIW_VALUE(MPI_ANY_TAG, int32_t),
#endif
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_VALUE(MPI_PROC_NULL, int32_t),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_IDENT, int32_t),
#else
  PLUMED_MPIW_VALUE(MPI_IDENT, int32_t),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_CONGRUENT, int32_t),
#else
  PLUMED_MPIW_VALUE(MPI_CONGRUENT, int32_t),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_SIMILAR, int32_t),
#else
  PLUMED_MPIW_VALUE(MPI_SIMILAR, int32_t),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_UNEQUAL, int32_t),
#else
  PLUMED_MPIW_VALUE(MPI_UNEQUAL, int32_t),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_CART, int32_t),
#else
  PLUMED_MPIW_VALUE(MPI_CART, int32_t),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_GRAPH, int32_t),
#else
  PLUMED_MPIW_VALUE(MPI_GRAPH, int32_t),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_MAX_PROCESSOR_NAME, int32_t),
#else
  PLUMED_MPIW_VALUE(MPI_MAX_PROCESSOR_NAME, int32_t),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_MAX_ERROR_STRING, int32_t),
#else
  PLUMED_MPIW_VALUE(MPI_MAX_ERROR_STRING, int32_t),
#endif
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_VALUE(MPI_BSEND_OVERHEAD, int32_t),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_SUCCESS, int32_t),
#else
  PLUMED_MPIW_VALUE(MPI_SUCCESS, int32_t),
#endif
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_VALUE(MPI_ERR_BUFFER, int32_t),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_VALUE(MPI_ERR_COUNT, int32_t),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_ERR_TYPE, int32_t),
#else
  PLUMED_MPIW_VALUE(MPI_ERR_TYPE, int32_t),
#endif
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_VALUE(MPI_ERR_TAG, int32_t),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_ERR_COMM, int32_t),
#else
  PLUMED_MPIW_VALUE(MPI_ERR_COMM, int32_t),
#endif
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_VALUE(MPI_ERR_RANK, int32_t),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_VALUE(MPI_ERR_REQUEST, int32_t),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_VALUE(MPI_ERR_ROOT, int32_t),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_ERR_GROUP, int32_t),
#else
  PLUMED_MPIW_VALUE(MPI_ERR_GROUP, int32_t),
#endif
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_VALUE(MPI_ERR_OP, int32_t),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_VALUE(MPI_ERR_TOPOLOGY, int32_t),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_ERR_DIMS, int32_t),
#else
  PLUMED_MPIW_VALUE(MPI_ERR_DIMS, int32_t),
#endif
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_VALUE(MPI_ERR_ARG, int32_t),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_ERR_UNKNOWN, int32_t),
#else
  PLUMED_MPIW_VALUE(MPI_ERR_UNKNOWN, int32_t),
#endif
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_VALUE(MPI_ERR_TRUNCATE, int32_t),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_VALUE(MPI_ERR_OTHER, int32_t),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_VALUE(MPI_ERR_INTERN, int32_t),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_VALUE(MPI_ERR_PENDING, int32_t),
#endif
#if (!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))
  PLUMED_MPIW_VALUE(MPI_BOTTOM, void *),
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 1) || (MPI_VERSION == 1 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_ERR_IN_STATUS, int32_t),
#else
  PLUMED_MPIW_VALUE(MPI_ERR_IN_STATUS, int32_t),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 1)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_SIGNED_CHAR, tMPI_Datatype),
#else
  PLUMED_MPIW_VALUE(MPI_SIGNED_CHAR, MPI_Datatype),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 2)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_C_BOOL, tMPI_Datatype),
#else
  PLUMED_MPIW_VALUE(MPI_C_BOOL, MPI_Datatype),
#endif
#endif
#if ((!defined(GMX_THREAD_MPI) || !(GMX_THREAD_MPI)) && (defined(MPI_VERSION) && defined(MPI_SUBVERSION) && ((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 2)))) || (defined(GMX_THREAD_MPI) && (GMX_THREAD_MPI) && PLUMED_MPIW_GROMACS_VERSION >= 20260000)
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  PLUMED_MPIW_VALUE(MPI_INT64_T, tMPI_Datatype),
#else
  PLUMED_MPIW_VALUE(MPI_INT64_T, MPI_Datatype),
#endif
#endif
  { (uint32_t)sizeof(plumed_mpiw_metadata_entry), PLUMED_MPIW_VALUE_OPAQUE, 0u, "identity:runtime", sizeof(plumed_mpiw_generic_fn), 0u, &plumed_mpiw_native_runtime_identity },
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  { (uint32_t)sizeof(plumed_mpiw_metadata_entry), PLUMED_MPIW_VALUE_UTF8, 0u, "capability:intracommunicators-only", sizeof("true"), 0u, "true" },
#endif
};

typedef union {
  long double long_double_value;
  long long long_long_value;
  void *pointer_value;
} plumed_mpiw_scratch_alignment;
typedef struct {
  char c;
  plumed_mpiw_scratch_alignment value;
} plumed_mpiw_scratch_alignment_offset;
#define PLUMED_MPIW_CHECK_NATIVE_ALIGNMENT(name, type) \
  typedef struct { char c; type value; } \
  PLUMED_MPIW_JOIN(plumed_mpiw_native_align_, name); \
  PLUMED_MPIW_STATIC_ASSERT(offsetof(PLUMED_MPIW_JOIN(plumed_mpiw_native_align_, name), value) <= \
                            offsetof(plumed_mpiw_scratch_alignment_offset, value), \
                            "native MPI type over-aligns provider scratch")
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
PLUMED_MPIW_CHECK_NATIVE_ALIGNMENT(comm, tMPI_Comm);
PLUMED_MPIW_CHECK_NATIVE_ALIGNMENT(group, tMPI_Group);
PLUMED_MPIW_CHECK_NATIVE_ALIGNMENT(datatype, tMPI_Datatype);
PLUMED_MPIW_CHECK_NATIVE_ALIGNMENT(request, tMPI_Request);
PLUMED_MPIW_CHECK_NATIVE_ALIGNMENT(op, tMPI_Op);
PLUMED_MPIW_CHECK_NATIVE_ALIGNMENT(errhandler, tMPI_Errhandler);
PLUMED_MPIW_CHECK_NATIVE_ALIGNMENT(status, tMPI_Status);
#else
PLUMED_MPIW_CHECK_NATIVE_ALIGNMENT(comm, MPI_Comm);
PLUMED_MPIW_CHECK_NATIVE_ALIGNMENT(group, MPI_Group);
PLUMED_MPIW_CHECK_NATIVE_ALIGNMENT(datatype, MPI_Datatype);
PLUMED_MPIW_CHECK_NATIVE_ALIGNMENT(request, MPI_Request);
PLUMED_MPIW_CHECK_NATIVE_ALIGNMENT(op, MPI_Op);
PLUMED_MPIW_CHECK_NATIVE_ALIGNMENT(errhandler, MPI_Errhandler);
PLUMED_MPIW_CHECK_NATIVE_ALIGNMENT(status, MPI_Status);
#endif
PLUMED_MPIW_STATIC_ASSERT(sizeof(int) == sizeof(int32_t) && (int)-1 < 0,
                          "native MPI int must be signed 32-bit");
static const plumed_mpiw_provider_descriptor plumed_mpiw_descriptor = {
  (uint32_t)sizeof(plumed_mpiw_provider_descriptor),
  PLUMED_MPIW_PROVIDER_ABI_MAJOR, 0u,
  (uint32_t)(
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
    1u
#else
    0u
#endif
  ),
  "PLUMED native MPI adapter",
#if defined(GMX_THREAD_MPI) && GMX_THREAD_MPI
  (uint32_t)sizeof(tMPI_Comm), (uint32_t)sizeof(tMPI_Group),
  (uint32_t)sizeof(tMPI_Datatype), (uint32_t)sizeof(tMPI_Request),
  (uint32_t)sizeof(tMPI_Op), (uint32_t)sizeof(tMPI_Errhandler),
  (uint32_t)sizeof(tMPI_Status),
  (uint32_t)offsetof(tMPI_Status, TMPI_SOURCE),
  (uint32_t)offsetof(tMPI_Status, TMPI_TAG),
  (uint32_t)offsetof(tMPI_Status, TMPI_ERROR),
  0u, 0u,
#else
  (uint32_t)sizeof(MPI_Comm), (uint32_t)sizeof(MPI_Group),
  (uint32_t)sizeof(MPI_Datatype), (uint32_t)sizeof(MPI_Request),
  (uint32_t)sizeof(MPI_Op), (uint32_t)sizeof(MPI_Errhandler),
  (uint32_t)sizeof(MPI_Status),
  (uint32_t)offsetof(MPI_Status, MPI_SOURCE),
  (uint32_t)offsetof(MPI_Status, MPI_TAG),
  (uint32_t)offsetof(MPI_Status, MPI_ERROR),
  (uint32_t)sizeof(MPI_Fint), 1u,
#endif
  (uint64_t)(sizeof(plumed_mpiw_functions) / sizeof(plumed_mpiw_functions[0])),
  plumed_mpiw_functions,
  (uint64_t)(sizeof(plumed_mpiw_metadata) / sizeof(plumed_mpiw_metadata[0])),
  plumed_mpiw_metadata
};
#undef PLUMED_MPIW_CHECK_NATIVE_ALIGNMENT
#undef PLUMED_MPIW_VALUE
#undef PLUMED_MPIW_STATIC_ASSERT
#undef PLUMED_MPIW_JOIN
#undef PLUMED_MPIW_JOIN_INNER

static const plumed_mpiw_provider_descriptor *plumed_mpiw_get_descriptor(void) {
  return &plumed_mpiw_descriptor;
}
