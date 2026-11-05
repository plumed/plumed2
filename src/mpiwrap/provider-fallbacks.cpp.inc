// Kernel-side algorithms shared by embedded and injected providers.
#include "mpiwrap.h"
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace PLMD {
namespace mpiwrap_fallback_detail {

bool ready() {
  return pointers::ready.load(std::memory_order_acquire);
}

bool request_constants() {
  return pointers::constant_REQUEST_NULL && pointers::constant_UNDEFINED &&
         pointers::constant_SUCCESS && pointers::constant_ANY_SOURCE && pointers::constant_ANY_TAG;
}

bool collective_primitives() {
  return ready() && pointers::constant_SUCCESS &&
         (pointers::constant_INT || pointers::constant_CHAR || pointers::constant_FLOAT ||
          pointers::constant_DOUBLE || pointers::constant_BYTE || pointers::constant_SHORT ||
          pointers::constant_LONG || pointers::constant_LONG_LONG || pointers::constant_LONG_LONG_INT ||
          pointers::constant_SIGNED_CHAR || pointers::constant_UNSIGNED_CHAR ||
          pointers::constant_UNSIGNED_SHORT || pointers::constant_UNSIGNED || pointers::constant_UNSIGNED_LONG ||
          pointers::constant_UNSIGNED_LONG_LONG || pointers::constant_LONG_DOUBLE ||
          pointers::constant_C_BOOL || pointers::constant_INT64_T) &&
         pointers::Comm_rank && pointers::Comm_size && pointers::Bcast &&
         (pointers::Comm_test_inter ||
          mpiwrap_provider_has_capability("capability:intracommunicators-only"));
}

// Basic C datatype extents are fixed by the MPI language binding. Compare the
// provider's handles; their numeric values and representations are immaterial.
// Compound pair types and custom typemaps are deliberately not approximated.
std::size_t basic_extent(MPI_Datatype datatype) {
#define PLUMED_MPIW_BASIC_EXTENT(name, type) \
  if(pointers::constant_ ## name && datatype==mpiwrap_ ## name()) return sizeof(type)
  PLUMED_MPIW_BASIC_EXTENT(CHAR,char);
  PLUMED_MPIW_BASIC_EXTENT(SIGNED_CHAR,signed char);
  PLUMED_MPIW_BASIC_EXTENT(UNSIGNED_CHAR,unsigned char);
  PLUMED_MPIW_BASIC_EXTENT(BYTE,unsigned char);
  PLUMED_MPIW_BASIC_EXTENT(SHORT,short);
  PLUMED_MPIW_BASIC_EXTENT(UNSIGNED_SHORT,unsigned short);
  PLUMED_MPIW_BASIC_EXTENT(INT,int);
  PLUMED_MPIW_BASIC_EXTENT(UNSIGNED,unsigned);
  PLUMED_MPIW_BASIC_EXTENT(LONG,long);
  PLUMED_MPIW_BASIC_EXTENT(UNSIGNED_LONG,unsigned long);
  PLUMED_MPIW_BASIC_EXTENT(LONG_LONG,long long);
  PLUMED_MPIW_BASIC_EXTENT(LONG_LONG_INT,long long);
  PLUMED_MPIW_BASIC_EXTENT(UNSIGNED_LONG_LONG,unsigned long long);
  PLUMED_MPIW_BASIC_EXTENT(FLOAT,float);
  PLUMED_MPIW_BASIC_EXTENT(DOUBLE,double);
  PLUMED_MPIW_BASIC_EXTENT(LONG_DOUBLE,long double);
  PLUMED_MPIW_BASIC_EXTENT(C_BOOL,bool);
  PLUMED_MPIW_BASIC_EXTENT(INT64_T,std::int64_t);
#undef PLUMED_MPIW_BASIC_EXTENT
  return 0;
}

std::size_t bytes(int count,std::size_t extent) {
  if(count<0 || (extent && static_cast<std::size_t>(count)>
                 std::numeric_limits<std::size_t>::max()/extent)) {
    throw std::invalid_argument("MPI fallback count or byte extent is invalid");
  }
  return static_cast<std::size_t>(count)*extent;
}

bool in_place(const void *buffer) {
  return pointers::constant_IN_PLACE && buffer==mpiwrap_IN_PLACE();
}

void intracommunicator(MPI_Comm comm) {
  if(pointers::Comm_test_inter) {
    int inter=0;
    const int rc=mpiwrap_Comm_test_inter(comm,&inter);
    if(rc!=mpiwrap_SUCCESS() || inter) {
      throw std::runtime_error("MPI collective fallback requires an intracommunicator");
    }
  } else if(!mpiwrap_provider_has_capability("capability:intracommunicators-only")) {
    throw std::runtime_error("MPI provider cannot establish the collective fallback communicator kind");
  }
}

void rank_size(MPI_Comm comm,int &rank,int &size) {
  intracommunicator(comm);
  if(mpiwrap_Comm_rank(comm,&rank)!=mpiwrap_SUCCESS() ||
      mpiwrap_Comm_size(comm,&size)!=mpiwrap_SUCCESS() || size<1 || rank<0 || rank>=size) {
    throw std::runtime_error("MPI collective fallback could not query the communicator");
  }
}

std::size_t collective_extent(MPI_Datatype type,const char *operation) {
  const std::size_t extent=basic_extent(type);
  if(!extent) {
    throw std::runtime_error(std::string(operation)+" fallback requires a supported basic C datatype (custom/compound types are unsupported)");
  }
  return extent;
}

int some(bool wait,int count,MPI_Request requests[],int *outcount,
         int indices[],MPI_Status statuses[]) {
  if(count<0 || !outcount || (count && (!requests || !indices))) {
    throw std::invalid_argument("MPI request completion fallback requires valid arrays");
  }
  if(mpiwrap_requests_all_null(count,requests)) {
    *outcount=mpiwrap_UNDEFINED();
    return mpiwrap_SUCCESS();
  }
  std::vector<unsigned char> completed_indices(static_cast<std::size_t>(count),0);
  int completed=0,result=mpiwrap_SUCCESS();
  while(completed<count) {
    int index=mpiwrap_UNDEFINED(),flag=0;
    MPI_Status status{};
    const int rc=wait && completed==0 ?
                 mpiwrap_Waitany(count,requests,&index,&status) :
                 mpiwrap_Testany(count,requests,&index,&flag,&status);
    if(wait && completed==0) {
      flag=index!=mpiwrap_UNDEFINED();
    }
    if(index==mpiwrap_UNDEFINED() || !flag) {
      if(rc!=mpiwrap_SUCCESS()) {
        result=rc;
      }
      break;
    }
    if(index<0 || index>=count || completed_indices[static_cast<std::size_t>(index)]) {
      throw std::runtime_error("MPI provider returned an invalid completion index");
    }
    completed_indices[static_cast<std::size_t>(index)]=1;
    indices[completed]=index;
    // Waitany/Testany report a single request's error directly. Some-functions
    // collect per-request errors and return MPI_ERR_IN_STATUS instead.
    if(rc!=mpiwrap_SUCCESS()) {
      if(rc!=mpiwrap_ERR_IN_STATUS()) {
        status.MPI_ERROR=rc;
      }
      result=mpiwrap_ERR_IN_STATUS();
    } else {
      status.MPI_ERROR=mpiwrap_SUCCESS();
    }
    if(statuses) {
      statuses[completed]=status;
    }
    ++completed;
  }
  *outcount=completed;
  return result;
}

}

bool mpiwrap_can_fallback_Wait() {
  return mpiwrap_fallback_detail::ready() && pointers::Waitany &&
         mpiwrap_fallback_detail::request_constants();
}

bool mpiwrap_can_fallback_Testsome() {
  return mpiwrap_fallback_detail::ready() && pointers::Testany &&
         pointers::constant_ERR_IN_STATUS && mpiwrap_fallback_detail::request_constants();
}

bool mpiwrap_can_fallback_Waitsome() {
  return mpiwrap_can_fallback_Testsome() && pointers::Waitany;
}

bool mpiwrap_can_fallback_Allgather() {
  return mpiwrap_fallback_detail::collective_primitives();
}

bool mpiwrap_can_fallback_Allgatherv() {
  return mpiwrap_fallback_detail::collective_primitives();
}

bool mpiwrap_requests_all_null(int count,const MPI_Request requests[]) {
  if(count<0 || (count && !requests)) {
    throw std::invalid_argument("MPI request completion requires a nonnegative count and request array");
  }
  if(count==0) {
    return true;
  }
  const MPI_Request null_request=mpiwrap_REQUEST_NULL();
  for(int i=0; i<count; ++i) {
    if(requests[i]!=null_request) {
      return false;
    }
  }
  return true;
}

void mpiwrap_set_empty_status(MPI_Status *status) {
  if(status) {
    // No fabricated provider-private status layout: raw_size=0 marks a
    // synthetic empty status, interpreted by PLUMED's status-query wrappers.
    *status=MPI_Status{};
    status->MPI_SOURCE=mpiwrap_ANY_SOURCE();
    status->MPI_TAG=mpiwrap_ANY_TAG();
    status->MPI_ERROR=mpiwrap_SUCCESS();
  }
}

bool mpiwrap_is_empty_status(const MPI_Status *status) {
  return status && status->raw_size==0 && pointers::constant_ANY_SOURCE &&
         pointers::constant_ANY_TAG && pointers::constant_SUCCESS &&
         status->MPI_SOURCE==mpiwrap_ANY_SOURCE() &&
         status->MPI_TAG==mpiwrap_ANY_TAG() && status->MPI_ERROR==mpiwrap_SUCCESS();
}

int mpiwrap_fallback_Wait(MPI_Request *request,MPI_Status *status) {
  if(!request) {
    throw std::invalid_argument("MPI_Wait requires a request");
  }
  int index=mpiwrap_UNDEFINED();
  const int rc=mpiwrap_Waitany(1,request,&index,status);
  if(index==mpiwrap_UNDEFINED() && rc==mpiwrap_SUCCESS()) {
    mpiwrap_set_empty_status(status);
  }
  return rc;
}

int mpiwrap_fallback_Waitsome(int count,MPI_Request requests[],int *outcount,
                              int indices[],MPI_Status statuses[]) {
  return mpiwrap_fallback_detail::some(true,count,requests,outcount,indices,statuses);
}

int mpiwrap_fallback_Testsome(int count,MPI_Request requests[],int *outcount,
                              int indices[],MPI_Status statuses[]) {
  return mpiwrap_fallback_detail::some(false,count,requests,outcount,indices,statuses);
}

int mpiwrap_fallback_Allgather(const void *send,int sendcount,MPI_Datatype sendtype,
                               void *recv,int recvcount,MPI_Datatype recvtype,MPI_Comm comm) {
  int rank=0,size=0;
  mpiwrap_fallback_detail::rank_size(comm,rank,size);
  const std::size_t extent=mpiwrap_fallback_detail::collective_extent(recvtype,"MPI_Allgather");
  const std::size_t segment_bytes=mpiwrap_fallback_detail::bytes(recvcount,extent);
  if(segment_bytes && static_cast<std::size_t>(size)>
      std::numeric_limits<std::size_t>::max()/segment_bytes) {
    throw std::invalid_argument("MPI_Allgather fallback receive extent overflows");
  }
  if(segment_bytes && !recv) {
    throw std::invalid_argument("MPI_Allgather fallback requires a receive buffer");
  }
  if(!mpiwrap_fallback_detail::in_place(send)) {
    mpiwrap_fallback_detail::collective_extent(sendtype,"MPI_Allgather");
    if(sendcount!=recvcount || sendtype!=recvtype || (segment_bytes && !send)) {
      throw std::invalid_argument("MPI_Allgather fallback requires matching basic datatypes/counts");
    }
    if(segment_bytes) {
      std::memmove(static_cast<unsigned char *>(recv)+static_cast<std::size_t>(rank)*segment_bytes,
                   send,segment_bytes);
    }
  }
  for(int root=0; root<size; ++root) {
    void *segment=segment_bytes ? static_cast<unsigned char *>(recv)+
                  static_cast<std::size_t>(root)*segment_bytes : recv;
    const int rc=mpiwrap_Bcast(segment,recvcount,recvtype,root,comm);
    if(rc!=mpiwrap_SUCCESS()) {
      return rc;
    }
  }
  return mpiwrap_SUCCESS();
}

int mpiwrap_fallback_Allgatherv(const void *send,int sendcount,MPI_Datatype sendtype,
                                void *recv,const int counts[],const int displs[],
                                MPI_Datatype recvtype,MPI_Comm comm) {
  int rank=0,size=0;
  mpiwrap_fallback_detail::rank_size(comm,rank,size);
  const std::size_t extent=mpiwrap_fallback_detail::collective_extent(recvtype,"MPI_Allgatherv");
  if(!counts || !displs) {
    throw std::invalid_argument("MPI_Allgatherv fallback requires counts and displacements");
  }
  // Validate every segment before entering any collective. Zero-count gaps
  // must remain untouched, and pointer arithmetic must not overflow.
  for(int i=0; i<size; ++i) {
    const std::size_t length=mpiwrap_fallback_detail::bytes(counts[i],extent);
    const std::size_t offset=mpiwrap_fallback_detail::bytes(displs[i],extent);
    if(offset>std::numeric_limits<std::size_t>::max()-length || (length && !recv)) {
      throw std::invalid_argument("MPI_Allgatherv fallback receive extent is invalid");
    }
  }
  if(!mpiwrap_fallback_detail::in_place(send)) {
    mpiwrap_fallback_detail::collective_extent(sendtype,"MPI_Allgatherv");
    if(sendcount!=counts[rank] || sendtype!=recvtype || (sendcount && !send)) {
      throw std::invalid_argument("MPI_Allgatherv fallback requires matching basic datatypes/counts");
    }
    const std::size_t length=mpiwrap_fallback_detail::bytes(sendcount,extent);
    if(length) {
      std::memmove(static_cast<unsigned char *>(recv)+
                   mpiwrap_fallback_detail::bytes(displs[rank],extent),send,length);
    }
  }
  for(int root=0; root<size; ++root) {
    void *segment=counts[root] ? static_cast<unsigned char *>(recv)+
                  mpiwrap_fallback_detail::bytes(displs[root],extent) : recv;
    const int rc=mpiwrap_Bcast(segment,counts[root],recvtype,root,comm);
    if(rc!=mpiwrap_SUCCESS()) {
      return rc;
    }
  }
  return mpiwrap_SUCCESS();
}

bool mpiwrap_try_singleton_reduce(const void *send,void *recv,int count,
                                  MPI_Datatype type,int root,MPI_Comm comm,int *result) {
  if(!pointers::Comm_size || root>0 || root< -1 || count<0) {
    return false;
  }
  // An intercommunicator reduction uses MPI_ROOT/MPI_PROC_NULL semantics;
  // a local group size of one does not make it a singleton collective.
  if(pointers::Comm_test_inter) {
    int inter=0;
    if(mpiwrap_Comm_test_inter(comm,&inter)!=mpiwrap_SUCCESS() || inter) {
      return false;
    }
  } else if(!mpiwrap_provider_has_capability("capability:intracommunicators-only")) {
    return false;
  }
  const std::size_t extent=mpiwrap_fallback_detail::basic_extent(type);
  if(!extent) {
    return false;
  }
  int size=0;
  if(mpiwrap_Comm_size(comm,&size)!=mpiwrap_SUCCESS() || size!=1) {
    return false;
  }
  const std::size_t length=mpiwrap_fallback_detail::bytes(count,extent);
  if(length && !recv) {
    throw std::invalid_argument("MPI reduction requires a receive buffer");
  }
  if(!mpiwrap_fallback_detail::in_place(send)) {
    if(length && !send) {
      throw std::invalid_argument("MPI reduction requires a send buffer");
    }
    if(length) {
      std::memmove(recv,send,length);
    }
  }
  *result=mpiwrap_SUCCESS();
  return true;
}

}
