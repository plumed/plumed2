// Generated stable thread-MPI constant accessors.
void get_MPI_FLOAT(void *output) {
  const thread_mpi::tMPI_Datatype value = TMPI_FLOAT;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_FLOAT = &get_MPI_FLOAT;
void get_MPI_DOUBLE(void *output) {
  const thread_mpi::tMPI_Datatype value = TMPI_DOUBLE;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_DOUBLE = &get_MPI_DOUBLE;
void get_MPI_INT(void *output) {
  const thread_mpi::tMPI_Datatype value = TMPI_INT;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_INT = &get_MPI_INT;
void get_MPI_CHAR(void *output) {
  const thread_mpi::tMPI_Datatype value = TMPI_CHAR;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_CHAR = &get_MPI_CHAR;
void get_MPI_UNSIGNED(void *output) {
  const thread_mpi::tMPI_Datatype value = TMPI_UNSIGNED;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_UNSIGNED = &get_MPI_UNSIGNED;
void get_MPI_UNSIGNED_LONG(void *output) {
  const thread_mpi::tMPI_Datatype value = TMPI_UNSIGNED_LONG;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_UNSIGNED_LONG = &get_MPI_UNSIGNED_LONG;
void get_MPI_UNSIGNED_LONG_LONG(void *output) {
  const thread_mpi::tMPI_Datatype value = TMPI_UNSIGNED_LONG_LONG;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_UNSIGNED_LONG_LONG = &get_MPI_UNSIGNED_LONG_LONG;
void get_MPI_LONG_DOUBLE(void *output) {
  const thread_mpi::tMPI_Datatype value = TMPI_LONG_DOUBLE;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_LONG_DOUBLE = &get_MPI_LONG_DOUBLE;
void get_MPI_COMM_WORLD(void *output) {
  const thread_mpi::tMPI_Comm value = TMPI_COMM_WORLD;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_COMM_WORLD = &get_MPI_COMM_WORLD;
void get_self(void *output) {
  thread_mpi::tMPI_Comm value{};
  int initialized = 0;
  if (thread_mpi::tMPI_Initialized(&initialized) == thread_mpi::TMPI_SUCCESS && initialized) {
    value = thread_mpi::tMPI_Get_comm_self();
  }
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_COMM_SELF = &get_self;
void get_MPI_COMM_NULL(void *output) {
  const thread_mpi::tMPI_Comm value = (thread_mpi::tMPI_Comm)0;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_COMM_NULL = &get_MPI_COMM_NULL;
void get_MPI_IN_PLACE(void *output) {
  const void * value = TMPI_IN_PLACE;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_IN_PLACE = &get_MPI_IN_PLACE;
void get_MPI_SUM(void *output) {
  const thread_mpi::tMPI_Op value = TMPI_SUM;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_SUM = &get_MPI_SUM;
void get_MPI_PROD(void *output) {
  const thread_mpi::tMPI_Op value = TMPI_PROD;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_PROD = &get_MPI_PROD;
void get_MPI_MAX(void *output) {
  const thread_mpi::tMPI_Op value = TMPI_MAX;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_MAX = &get_MPI_MAX;
void get_MPI_MIN(void *output) {
  const thread_mpi::tMPI_Op value = TMPI_MIN;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_MIN = &get_MPI_MIN;
void get_MPI_STATUS_IGNORE(void *output) {
  const thread_mpi::tMPI_Status * value = TMPI_STATUS_IGNORE;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_STATUS_IGNORE = &get_MPI_STATUS_IGNORE;
void get_MPI_STATUSES_IGNORE(void *output) {
  const thread_mpi::tMPI_Status * value = TMPI_STATUSES_IGNORE;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_STATUSES_IGNORE = &get_MPI_STATUSES_IGNORE;
void get_MPI_GROUP_NULL(void *output) {
  const thread_mpi::tMPI_Group value = (thread_mpi::tMPI_Group)0;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_GROUP_NULL = &get_MPI_GROUP_NULL;
void get_MPI_GROUP_EMPTY(void *output) {
  const thread_mpi::tMPI_Group value = TMPI_GROUP_EMPTY;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_GROUP_EMPTY = &get_MPI_GROUP_EMPTY;
void get_MPI_REQUEST_NULL(void *output) {
  const thread_mpi::tMPI_Request value = TMPI_REQUEST_NULL;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_REQUEST_NULL = &get_MPI_REQUEST_NULL;
void get_MPI_DATATYPE_NULL(void *output) {
  const thread_mpi::tMPI_Datatype value = (thread_mpi::tMPI_Datatype)0;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_DATATYPE_NULL = &get_MPI_DATATYPE_NULL;
void get_MPI_ERRORS_ARE_FATAL(void *output) {
  const thread_mpi::tMPI_Errhandler value = TMPI_ERRORS_ARE_FATAL;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_ERRORS_ARE_FATAL = &get_MPI_ERRORS_ARE_FATAL;
void get_MPI_ERRORS_RETURN(void *output) {
  const thread_mpi::tMPI_Errhandler value = TMPI_ERRORS_RETURN;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_ERRORS_RETURN = &get_MPI_ERRORS_RETURN;
void get_MPI_BYTE(void *output) {
  const thread_mpi::tMPI_Datatype value = TMPI_BYTE;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_BYTE = &get_MPI_BYTE;
void get_MPI_SHORT(void *output) {
  const thread_mpi::tMPI_Datatype value = TMPI_SHORT;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_SHORT = &get_MPI_SHORT;
void get_MPI_LONG(void *output) {
  const thread_mpi::tMPI_Datatype value = TMPI_LONG;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_LONG = &get_MPI_LONG;
void get_MPI_UNSIGNED_CHAR(void *output) {
  const thread_mpi::tMPI_Datatype value = TMPI_UNSIGNED_CHAR;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_UNSIGNED_CHAR = &get_MPI_UNSIGNED_CHAR;
void get_MPI_UNSIGNED_SHORT(void *output) {
  const thread_mpi::tMPI_Datatype value = TMPI_UNSIGNED_SHORT;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_UNSIGNED_SHORT = &get_MPI_UNSIGNED_SHORT;
void get_MPI_LONG_LONG_INT(void *output) {
  const thread_mpi::tMPI_Datatype value = TMPI_LONG_LONG_INT;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_LONG_LONG_INT = &get_MPI_LONG_LONG_INT;
void get_MPI_LONG_LONG(void *output) {
  const thread_mpi::tMPI_Datatype value = TMPI_LONG_LONG;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_LONG_LONG = &get_MPI_LONG_LONG;
void get_MPI_LAND(void *output) {
  const thread_mpi::tMPI_Op value = TMPI_LAND;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_LAND = &get_MPI_LAND;
void get_MPI_BAND(void *output) {
  const thread_mpi::tMPI_Op value = TMPI_BAND;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_BAND = &get_MPI_BAND;
void get_MPI_LOR(void *output) {
  const thread_mpi::tMPI_Op value = TMPI_LOR;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_LOR = &get_MPI_LOR;
void get_MPI_BOR(void *output) {
  const thread_mpi::tMPI_Op value = TMPI_BOR;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_BOR = &get_MPI_BOR;
void get_MPI_LXOR(void *output) {
  const thread_mpi::tMPI_Op value = TMPI_LXOR;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_LXOR = &get_MPI_LXOR;
void get_MPI_BXOR(void *output) {
  const thread_mpi::tMPI_Op value = TMPI_BXOR;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_BXOR = &get_MPI_BXOR;
void get_MPI_UNDEFINED(void *output) {
  const int value = TMPI_UNDEFINED;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_UNDEFINED = &get_MPI_UNDEFINED;
void get_MPI_ANY_SOURCE(void *output) {
  const int value = TMPI_ANY_SOURCE;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_ANY_SOURCE = &get_MPI_ANY_SOURCE;
void get_MPI_ANY_TAG(void *output) {
  const int value = TMPI_ANY_TAG;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_ANY_TAG = &get_MPI_ANY_TAG;
void get_MPI_IDENT(void *output) {
  const int value = TMPI_IDENT;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_IDENT = &get_MPI_IDENT;
void get_MPI_CONGRUENT(void *output) {
  const int value = TMPI_CONGRUENT;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_CONGRUENT = &get_MPI_CONGRUENT;
void get_MPI_SIMILAR(void *output) {
  const int value = TMPI_SIMILAR;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_SIMILAR = &get_MPI_SIMILAR;
void get_MPI_UNEQUAL(void *output) {
  const int value = TMPI_UNEQUAL;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_UNEQUAL = &get_MPI_UNEQUAL;
void get_MPI_CART(void *output) {
  const int value = TMPI_CART;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_CART = &get_MPI_CART;
void get_MPI_GRAPH(void *output) {
  const int value = TMPI_GRAPH;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_GRAPH = &get_MPI_GRAPH;
void get_MPI_MAX_PROCESSOR_NAME(void *output) {
  const int value = TMPI_MAX_PROCESSOR_NAME;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_MAX_PROCESSOR_NAME = &get_MPI_MAX_PROCESSOR_NAME;
void get_MPI_MAX_ERROR_STRING(void *output) {
  const int value = TMPI_MAX_ERROR_STRING;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_MAX_ERROR_STRING = &get_MPI_MAX_ERROR_STRING;
void get_MPI_SUCCESS(void *output) {
  const int value = TMPI_SUCCESS;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_SUCCESS = &get_MPI_SUCCESS;
void get_MPI_ERR_TYPE(void *output) {
  const int value = TMPI_ERR_TYPE;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_ERR_TYPE = &get_MPI_ERR_TYPE;
void get_MPI_ERR_COMM(void *output) {
  const int value = TMPI_ERR_COMM;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_ERR_COMM = &get_MPI_ERR_COMM;
void get_MPI_ERR_GROUP(void *output) {
  const int value = TMPI_ERR_GROUP;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_ERR_GROUP = &get_MPI_ERR_GROUP;
void get_MPI_ERR_DIMS(void *output) {
  const int value = TMPI_ERR_DIMS;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_ERR_DIMS = &get_MPI_ERR_DIMS;
void get_MPI_ERR_UNKNOWN(void *output) {
  const int value = TMPI_ERR_UNKNOWN;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_ERR_UNKNOWN = &get_MPI_ERR_UNKNOWN;
void get_MPI_ERR_IN_STATUS(void *output) {
  const int value = TMPI_ERR_IN_STATUS;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_ERR_IN_STATUS = &get_MPI_ERR_IN_STATUS;
void get_MPI_SIGNED_CHAR(void *output) {
  const thread_mpi::tMPI_Datatype value = TMPI_SIGNED_CHAR;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_SIGNED_CHAR = &get_MPI_SIGNED_CHAR;
void get_MPI_C_BOOL(void *output) {
  const thread_mpi::tMPI_Datatype value = TMPI_C_BOOL;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_C_BOOL = &get_MPI_C_BOOL;
void get_MPI_INT64_T(void *output) {
  const thread_mpi::tMPI_Datatype value = TMPI_INT64_T;
  std::memcpy(output, &value, sizeof(value));
}
static const plumed_mpiw_constant_accessor_fn constant_accessor_MPI_INT64_T = &get_MPI_INT64_T;
