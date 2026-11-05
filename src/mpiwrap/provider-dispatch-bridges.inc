// Generated canonical-dispatch bridges for the typed frontend.
static int32_t plumed_mpiw_bridge_MPI_Initialized(int32_t * flag) {
  Registry &registry = currentRegistry();
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Initialized", 0xd4498369u, 1u, 1u);
  plumed_mpiw_transport_MPI_Initialized_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Initialized_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Initialized");
  }
  std::int32_t native_return_value = 0;
  dispatch(flag, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Comm_dup(plumed_mpiw_handle_storage comm, plumed_mpiw_handle_storage * newcomm) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Comm_dup", native_comm.data());
  NativeBuffer native_newcomm(registry.nativeHandleSize("Comm"));
  void *native_pointer_newcomm = newcomm ? native_newcomm.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Comm_dup", 0xf0dc6608u, 1u, 2u);
  plumed_mpiw_transport_MPI_Comm_dup_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Comm_dup_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Comm_dup");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm.data(), native_pointer_newcomm, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (newcomm && valid_output) {
    std::memset(newcomm->bytes, 0, sizeof(newcomm->bytes));
    std::memcpy(newcomm->bytes, native_newcomm.data(), native_newcomm.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Comm_free(plumed_mpiw_handle_storage * comm) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  if (comm) {
    std::memcpy(native_comm.data(), comm->bytes, native_comm.size());
  }
  if (comm) {
    validateInputComm(registry, "MPI_Comm_free", native_comm.data());
  }
  void *native_pointer_comm = comm ? native_comm.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Comm_free", 0x4da4e1b2u, 1u, 1u);
  plumed_mpiw_transport_MPI_Comm_free_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Comm_free_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Comm_free");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_pointer_comm, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (comm && valid_output) {
    std::memset(comm->bytes, 0, sizeof(comm->bytes));
    std::memcpy(comm->bytes, native_comm.data(), native_comm.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Comm_rank(plumed_mpiw_handle_storage comm, int32_t * rank) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Comm_rank", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Comm_rank", 0x4e7672aau, 1u, 2u);
  plumed_mpiw_transport_MPI_Comm_rank_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Comm_rank_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Comm_rank");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm.data(), rank, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Comm_size(plumed_mpiw_handle_storage comm, int32_t * size) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Comm_size", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Comm_size", 0xda2369ecu, 1u, 2u);
  plumed_mpiw_transport_MPI_Comm_size_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Comm_size_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Comm_size");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm.data(), size, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Comm_split(plumed_mpiw_handle_storage comm, int32_t color, int32_t key, plumed_mpiw_handle_storage * newcomm) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Comm_split", native_comm.data());
  int native_color = static_cast<int>(color);
  int native_key = static_cast<int>(key);
  NativeBuffer native_newcomm(registry.nativeHandleSize("Comm"));
  void *native_pointer_newcomm = newcomm ? native_newcomm.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Comm_split", 0x97f1adf0u, 1u, 4u);
  plumed_mpiw_transport_MPI_Comm_split_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Comm_split_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Comm_split");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm.data(), &native_color, &native_key, native_pointer_newcomm, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (newcomm && valid_output) {
    std::memset(newcomm->bytes, 0, sizeof(newcomm->bytes));
    std::memcpy(newcomm->bytes, native_newcomm.data(), native_newcomm.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Abort(plumed_mpiw_handle_storage comm, int32_t errorcode) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Abort", native_comm.data());
  int native_errorcode = static_cast<int>(errorcode);
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Abort", 0x84f458d6u, 1u, 2u);
  plumed_mpiw_transport_MPI_Abort_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Abort_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Abort");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm.data(), &native_errorcode, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Bcast(void * buffer, int32_t count, plumed_mpiw_handle_storage datatype, int32_t root, plumed_mpiw_handle_storage comm) {
  Registry &registry = currentRegistry();
  int native_count = static_cast<int>(count);
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_datatype.data(), datatype.bytes, native_datatype.size());
  int native_root = static_cast<int>(root);
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Bcast", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Bcast", 0xe809e1cau, 1u, 5u);
  plumed_mpiw_transport_MPI_Bcast_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Bcast_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Bcast");
  }
  std::int32_t native_return_value = 0;
  dispatch(buffer, &native_count, native_datatype.data(), &native_root, native_comm.data(), &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Scatter(const void * sendbuf, int32_t sendcount, plumed_mpiw_handle_storage sendtype, void * recvbuf, int32_t recvcount, plumed_mpiw_handle_storage recvtype, int32_t root, plumed_mpiw_handle_storage comm) {
  Registry &registry = currentRegistry();
  int native_sendcount = static_cast<int>(sendcount);
  NativeBuffer native_sendtype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_sendtype.data(), sendtype.bytes, native_sendtype.size());
  int native_recvcount = static_cast<int>(recvcount);
  NativeBuffer native_recvtype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_recvtype.data(), recvtype.bytes, native_recvtype.size());
  int native_root = static_cast<int>(root);
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Scatter", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Scatter", 0xb2fa0f58u, 1u, 8u);
  plumed_mpiw_transport_MPI_Scatter_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Scatter_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Scatter");
  }
  std::int32_t native_return_value = 0;
  dispatch(const_cast<void *>(static_cast<const void *>(sendbuf)), &native_sendcount, native_sendtype.data(), recvbuf, &native_recvcount, native_recvtype.data(), &native_root, native_comm.data(), &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Gather(const void * sendbuf, int32_t sendcount, plumed_mpiw_handle_storage sendtype, void * recvbuf, int32_t recvcount, plumed_mpiw_handle_storage recvtype, int32_t root, plumed_mpiw_handle_storage comm) {
  Registry &registry = currentRegistry();
  int native_sendcount = static_cast<int>(sendcount);
  NativeBuffer native_sendtype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_sendtype.data(), sendtype.bytes, native_sendtype.size());
  int native_recvcount = static_cast<int>(recvcount);
  NativeBuffer native_recvtype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_recvtype.data(), recvtype.bytes, native_recvtype.size());
  int native_root = static_cast<int>(root);
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Gather", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Gather", 0xb1f903d9u, 1u, 8u);
  plumed_mpiw_transport_MPI_Gather_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Gather_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Gather");
  }
  std::int32_t native_return_value = 0;
  dispatch(const_cast<void *>(static_cast<const void *>(sendbuf)), &native_sendcount, native_sendtype.data(), recvbuf, &native_recvcount, native_recvtype.data(), &native_root, native_comm.data(), &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Allreduce(const void * sendbuf, void * recvbuf, int32_t count, plumed_mpiw_handle_storage datatype, plumed_mpiw_handle_storage op, plumed_mpiw_handle_storage comm) {
  Registry &registry = currentRegistry();
  int native_count = static_cast<int>(count);
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_datatype.data(), datatype.bytes, native_datatype.size());
  NativeBuffer native_op(registry.nativeHandleSize("Op"));
  std::memcpy(native_op.data(), op.bytes, native_op.size());
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Allreduce", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Allreduce", 0x788900dbu, 1u, 6u);
  plumed_mpiw_transport_MPI_Allreduce_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Allreduce_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Allreduce");
  }
  std::int32_t native_return_value = 0;
  dispatch(const_cast<void *>(static_cast<const void *>(sendbuf)), recvbuf, &native_count, native_datatype.data(), native_op.data(), native_comm.data(), &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Isend(const void * buf, int32_t count, plumed_mpiw_handle_storage datatype, int32_t dest, int32_t tag, plumed_mpiw_handle_storage comm, plumed_mpiw_handle_storage * request) {
  Registry &registry = currentRegistry();
  int native_count = static_cast<int>(count);
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_datatype.data(), datatype.bytes, native_datatype.size());
  int native_dest = static_cast<int>(dest);
  int native_tag = static_cast<int>(tag);
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Isend", native_comm.data());
  NativeBuffer native_request(registry.nativeHandleSize("Request"));
  void *native_pointer_request = request ? native_request.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Isend", 0xf5e9d915u, 1u, 7u);
  plumed_mpiw_transport_MPI_Isend_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Isend_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Isend");
  }
  std::int32_t native_return_value = 0;
  dispatch(const_cast<void *>(static_cast<const void *>(buf)), &native_count, native_datatype.data(), &native_dest, &native_tag, native_comm.data(), native_pointer_request, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (request && valid_output) {
    std::memset(request->bytes, 0, sizeof(request->bytes));
    std::memcpy(request->bytes, native_request.data(), native_request.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Wait(plumed_mpiw_handle_storage * request, plumed_mpiw_status * status) {
  Registry &registry = currentRegistry();
  NativeBuffer native_request(registry.nativeHandleSize("Request"));
  if (request) {
    std::memcpy(native_request.data(), request->bytes, native_request.size());
  }
  void *native_pointer_request = request ? native_request.data() : nullptr;
  NativeBuffer native_status(registry.nativeStatusSize());
  if (status && false) {
    encodeNativeStatus(registry, *status, native_status.data());
  }
  void *native_pointer_status = status ? native_status.data() : nativeStatusIgnore(registry);
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Wait", 0x4b15813du, 1u, 2u);
  plumed_mpiw_transport_MPI_Wait_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Wait_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Wait");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_pointer_request, native_pointer_status, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (request && valid_output) {
    std::memset(request->bytes, 0, sizeof(request->bytes));
    std::memcpy(request->bytes, native_request.data(), native_request.size());
  }
  if (status && valid_output) {
    decodeNativeStatus(registry, native_status.data(), status);
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Allgather(const void * sendbuf, int32_t sendcount, plumed_mpiw_handle_storage sendtype, void * recvbuf, int32_t recvcount, plumed_mpiw_handle_storage recvtype, plumed_mpiw_handle_storage comm) {
  Registry &registry = currentRegistry();
  int native_sendcount = static_cast<int>(sendcount);
  NativeBuffer native_sendtype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_sendtype.data(), sendtype.bytes, native_sendtype.size());
  int native_recvcount = static_cast<int>(recvcount);
  NativeBuffer native_recvtype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_recvtype.data(), recvtype.bytes, native_recvtype.size());
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Allgather", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Allgather", 0xab8b26c9u, 1u, 7u);
  plumed_mpiw_transport_MPI_Allgather_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Allgather_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Allgather");
  }
  std::int32_t native_return_value = 0;
  dispatch(const_cast<void *>(static_cast<const void *>(sendbuf)), &native_sendcount, native_sendtype.data(), recvbuf, &native_recvcount, native_recvtype.data(), native_comm.data(), &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Allgatherv(const void * sendbuf, int32_t sendcount, plumed_mpiw_handle_storage sendtype, void * recvbuf, const int32_t recvcounts[], const int32_t displs[], plumed_mpiw_handle_storage recvtype, plumed_mpiw_handle_storage comm) {
  Registry &registry = currentRegistry();
  int native_sendcount = static_cast<int>(sendcount);
  NativeBuffer native_sendtype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_sendtype.data(), sendtype.bytes, native_sendtype.size());
  NativeBuffer native_recvtype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_recvtype.data(), recvtype.bytes, native_recvtype.size());
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Allgatherv", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Allgatherv", 0x0f1b22f3u, 1u, 8u);
  plumed_mpiw_transport_MPI_Allgatherv_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Allgatherv_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Allgatherv");
  }
  std::int32_t native_return_value = 0;
  dispatch(const_cast<void *>(static_cast<const void *>(sendbuf)), &native_sendcount, native_sendtype.data(), recvbuf, const_cast<void *>(static_cast<const void *>(recvcounts)), const_cast<void *>(static_cast<const void *>(displs)), native_recvtype.data(), native_comm.data(), &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Recv(void * buf, int32_t count, plumed_mpiw_handle_storage datatype, int32_t source, int32_t tag, plumed_mpiw_handle_storage comm, plumed_mpiw_status * status) {
  Registry &registry = currentRegistry();
  int native_count = static_cast<int>(count);
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_datatype.data(), datatype.bytes, native_datatype.size());
  int native_source = static_cast<int>(source);
  int native_tag = static_cast<int>(tag);
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Recv", native_comm.data());
  NativeBuffer native_status(registry.nativeStatusSize());
  if (status && false) {
    encodeNativeStatus(registry, *status, native_status.data());
  }
  void *native_pointer_status = status ? native_status.data() : nativeStatusIgnore(registry);
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Recv", 0xdd56a4bcu, 1u, 7u);
  plumed_mpiw_transport_MPI_Recv_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Recv_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Recv");
  }
  std::int32_t native_return_value = 0;
  dispatch(buf, &native_count, native_datatype.data(), &native_source, &native_tag, native_comm.data(), native_pointer_status, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (status && valid_output) {
    decodeNativeStatus(registry, native_status.data(), status);
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Barrier(plumed_mpiw_handle_storage comm) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Barrier", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Barrier", 0xe7a052f4u, 1u, 1u);
  plumed_mpiw_transport_MPI_Barrier_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Barrier_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Barrier");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm.data(), &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Get_count(const plumed_mpiw_status * status, plumed_mpiw_handle_storage datatype, int32_t * count) {
  Registry &registry = currentRegistry();
  NativeBuffer native_status(registry.nativeStatusSize());
  if (status && true) {
    encodeNativeStatus(registry, *status, native_status.data());
  }
  void *native_pointer_status = status ? native_status.data() : nativeStatusIgnore(registry);
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_datatype.data(), datatype.bytes, native_datatype.size());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Get_count", 0x71967a66u, 1u, 3u);
  plumed_mpiw_transport_MPI_Get_count_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Get_count_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Get_count");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_pointer_status, native_datatype.data(), count, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Sendrecv(const void * sendbuf, int32_t sendcount, plumed_mpiw_handle_storage sendtype, int32_t dest, int32_t sendtag, void * recvbuf, int32_t recvcount, plumed_mpiw_handle_storage recvtype, int32_t source, int32_t recvtag, plumed_mpiw_handle_storage comm, plumed_mpiw_status * status) {
  Registry &registry = currentRegistry();
  int native_sendcount = static_cast<int>(sendcount);
  NativeBuffer native_sendtype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_sendtype.data(), sendtype.bytes, native_sendtype.size());
  int native_dest = static_cast<int>(dest);
  int native_sendtag = static_cast<int>(sendtag);
  int native_recvcount = static_cast<int>(recvcount);
  NativeBuffer native_recvtype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_recvtype.data(), recvtype.bytes, native_recvtype.size());
  int native_source = static_cast<int>(source);
  int native_recvtag = static_cast<int>(recvtag);
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Sendrecv", native_comm.data());
  NativeBuffer native_status(registry.nativeStatusSize());
  if (status && false) {
    encodeNativeStatus(registry, *status, native_status.data());
  }
  void *native_pointer_status = status ? native_status.data() : nativeStatusIgnore(registry);
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Sendrecv", 0xb09a6fc1u, 1u, 12u);
  plumed_mpiw_transport_MPI_Sendrecv_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Sendrecv_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Sendrecv");
  }
  std::int32_t native_return_value = 0;
  dispatch(const_cast<void *>(static_cast<const void *>(sendbuf)), &native_sendcount, native_sendtype.data(), &native_dest, &native_sendtag, recvbuf, &native_recvcount, native_recvtype.data(), &native_source, &native_recvtag, native_comm.data(), native_pointer_status, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (status && valid_output) {
    decodeNativeStatus(registry, native_status.data(), status);
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Error_string(int32_t errorcode, char * string, int32_t * resultlen) {
  Registry &registry = currentRegistry();
  int native_errorcode = static_cast<int>(errorcode);
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Error_string", 0xd272abecu, 1u, 3u);
  plumed_mpiw_transport_MPI_Error_string_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Error_string_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Error_string");
  }
  std::int32_t native_return_value = 0;
  dispatch(&native_errorcode, string, resultlen, &native_return_value);
  return native_return_value;
}
static plumed_mpiw_handle_storage plumed_mpiw_bridge_MPI_Comm_f2c(int32_t comm) {
  Registry &registry = currentRegistry();
  MPI_Fint native_comm = static_cast<MPI_Fint>(comm);
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Comm_f2c", 0x9a9e7fcdu, 1u, 1u);
  plumed_mpiw_transport_MPI_Comm_f2c_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Comm_f2c_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Comm_f2c");
  }
  NativeBuffer native_result(registry.nativeHandleSize("Comm"));
  dispatch(&native_comm, native_result.data());
  plumed_mpiw_handle_storage result{};
  std::memcpy(result.bytes, native_result.data(), native_result.size());
  return result;
}
static int32_t plumed_mpiw_bridge_MPI_Send(const void * buf, int32_t count, plumed_mpiw_handle_storage datatype, int32_t dest, int32_t tag, plumed_mpiw_handle_storage comm) {
  Registry &registry = currentRegistry();
  int native_count = static_cast<int>(count);
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_datatype.data(), datatype.bytes, native_datatype.size());
  int native_dest = static_cast<int>(dest);
  int native_tag = static_cast<int>(tag);
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Send", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Send", 0xbd2b6912u, 1u, 6u);
  plumed_mpiw_transport_MPI_Send_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Send_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Send");
  }
  std::int32_t native_return_value = 0;
  dispatch(const_cast<void *>(static_cast<const void *>(buf)), &native_count, native_datatype.data(), &native_dest, &native_tag, native_comm.data(), &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Bsend(const void * buf, int32_t count, plumed_mpiw_handle_storage datatype, int32_t dest, int32_t tag, plumed_mpiw_handle_storage comm) {
  Registry &registry = currentRegistry();
  int native_count = static_cast<int>(count);
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_datatype.data(), datatype.bytes, native_datatype.size());
  int native_dest = static_cast<int>(dest);
  int native_tag = static_cast<int>(tag);
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Bsend", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Bsend", 0xa9dbf685u, 1u, 6u);
  plumed_mpiw_transport_MPI_Bsend_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Bsend_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Bsend");
  }
  std::int32_t native_return_value = 0;
  dispatch(const_cast<void *>(static_cast<const void *>(buf)), &native_count, native_datatype.data(), &native_dest, &native_tag, native_comm.data(), &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Ssend(const void * buf, int32_t count, plumed_mpiw_handle_storage datatype, int32_t dest, int32_t tag, plumed_mpiw_handle_storage comm) {
  Registry &registry = currentRegistry();
  int native_count = static_cast<int>(count);
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_datatype.data(), datatype.bytes, native_datatype.size());
  int native_dest = static_cast<int>(dest);
  int native_tag = static_cast<int>(tag);
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Ssend", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Ssend", 0x5edcf0deu, 1u, 6u);
  plumed_mpiw_transport_MPI_Ssend_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Ssend_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Ssend");
  }
  std::int32_t native_return_value = 0;
  dispatch(const_cast<void *>(static_cast<const void *>(buf)), &native_count, native_datatype.data(), &native_dest, &native_tag, native_comm.data(), &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Rsend(const void * buf, int32_t count, plumed_mpiw_handle_storage datatype, int32_t dest, int32_t tag, plumed_mpiw_handle_storage comm) {
  Registry &registry = currentRegistry();
  int native_count = static_cast<int>(count);
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_datatype.data(), datatype.bytes, native_datatype.size());
  int native_dest = static_cast<int>(dest);
  int native_tag = static_cast<int>(tag);
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Rsend", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Rsend", 0xbfa09b06u, 1u, 6u);
  plumed_mpiw_transport_MPI_Rsend_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Rsend_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Rsend");
  }
  std::int32_t native_return_value = 0;
  dispatch(const_cast<void *>(static_cast<const void *>(buf)), &native_count, native_datatype.data(), &native_dest, &native_tag, native_comm.data(), &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Buffer_attach(void * buffer, int32_t size) {
  Registry &registry = currentRegistry();
  int native_size = static_cast<int>(size);
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Buffer_attach", 0xe00bd480u, 1u, 2u);
  plumed_mpiw_transport_MPI_Buffer_attach_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Buffer_attach_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Buffer_attach");
  }
  std::int32_t native_return_value = 0;
  dispatch(buffer, &native_size, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Buffer_detach(void * buffer_addr, int32_t * size) {
  Registry &registry = currentRegistry();
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Buffer_detach", 0x74dec40fu, 1u, 2u);
  plumed_mpiw_transport_MPI_Buffer_detach_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Buffer_detach_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Buffer_detach");
  }
  std::int32_t native_return_value = 0;
  dispatch(buffer_addr, size, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Ibsend(const void * buf, int32_t count, plumed_mpiw_handle_storage datatype, int32_t dest, int32_t tag, plumed_mpiw_handle_storage comm, plumed_mpiw_handle_storage * request) {
  Registry &registry = currentRegistry();
  int native_count = static_cast<int>(count);
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_datatype.data(), datatype.bytes, native_datatype.size());
  int native_dest = static_cast<int>(dest);
  int native_tag = static_cast<int>(tag);
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Ibsend", native_comm.data());
  NativeBuffer native_request(registry.nativeHandleSize("Request"));
  void *native_pointer_request = request ? native_request.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Ibsend", 0xa4d683afu, 1u, 7u);
  plumed_mpiw_transport_MPI_Ibsend_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Ibsend_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Ibsend");
  }
  std::int32_t native_return_value = 0;
  dispatch(const_cast<void *>(static_cast<const void *>(buf)), &native_count, native_datatype.data(), &native_dest, &native_tag, native_comm.data(), native_pointer_request, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (request && valid_output) {
    std::memset(request->bytes, 0, sizeof(request->bytes));
    std::memcpy(request->bytes, native_request.data(), native_request.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Issend(const void * buf, int32_t count, plumed_mpiw_handle_storage datatype, int32_t dest, int32_t tag, plumed_mpiw_handle_storage comm, plumed_mpiw_handle_storage * request) {
  Registry &registry = currentRegistry();
  int native_count = static_cast<int>(count);
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_datatype.data(), datatype.bytes, native_datatype.size());
  int native_dest = static_cast<int>(dest);
  int native_tag = static_cast<int>(tag);
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Issend", native_comm.data());
  NativeBuffer native_request(registry.nativeHandleSize("Request"));
  void *native_pointer_request = request ? native_request.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Issend", 0x8275e866u, 1u, 7u);
  plumed_mpiw_transport_MPI_Issend_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Issend_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Issend");
  }
  std::int32_t native_return_value = 0;
  dispatch(const_cast<void *>(static_cast<const void *>(buf)), &native_count, native_datatype.data(), &native_dest, &native_tag, native_comm.data(), native_pointer_request, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (request && valid_output) {
    std::memset(request->bytes, 0, sizeof(request->bytes));
    std::memcpy(request->bytes, native_request.data(), native_request.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Irsend(const void * buf, int32_t count, plumed_mpiw_handle_storage datatype, int32_t dest, int32_t tag, plumed_mpiw_handle_storage comm, plumed_mpiw_handle_storage * request) {
  Registry &registry = currentRegistry();
  int native_count = static_cast<int>(count);
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_datatype.data(), datatype.bytes, native_datatype.size());
  int native_dest = static_cast<int>(dest);
  int native_tag = static_cast<int>(tag);
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Irsend", native_comm.data());
  NativeBuffer native_request(registry.nativeHandleSize("Request"));
  void *native_pointer_request = request ? native_request.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Irsend", 0xb13d1d70u, 1u, 7u);
  plumed_mpiw_transport_MPI_Irsend_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Irsend_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Irsend");
  }
  std::int32_t native_return_value = 0;
  dispatch(const_cast<void *>(static_cast<const void *>(buf)), &native_count, native_datatype.data(), &native_dest, &native_tag, native_comm.data(), native_pointer_request, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (request && valid_output) {
    std::memset(request->bytes, 0, sizeof(request->bytes));
    std::memcpy(request->bytes, native_request.data(), native_request.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Irecv(void * buf, int32_t count, plumed_mpiw_handle_storage datatype, int32_t source, int32_t tag, plumed_mpiw_handle_storage comm, plumed_mpiw_handle_storage * request) {
  Registry &registry = currentRegistry();
  int native_count = static_cast<int>(count);
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_datatype.data(), datatype.bytes, native_datatype.size());
  int native_source = static_cast<int>(source);
  int native_tag = static_cast<int>(tag);
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Irecv", native_comm.data());
  NativeBuffer native_request(registry.nativeHandleSize("Request"));
  void *native_pointer_request = request ? native_request.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Irecv", 0xb69a0083u, 1u, 7u);
  plumed_mpiw_transport_MPI_Irecv_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Irecv_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Irecv");
  }
  std::int32_t native_return_value = 0;
  dispatch(buf, &native_count, native_datatype.data(), &native_source, &native_tag, native_comm.data(), native_pointer_request, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (request && valid_output) {
    std::memset(request->bytes, 0, sizeof(request->bytes));
    std::memcpy(request->bytes, native_request.data(), native_request.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Test(plumed_mpiw_handle_storage * request, int32_t * flag, plumed_mpiw_status * status) {
  Registry &registry = currentRegistry();
  NativeBuffer native_request(registry.nativeHandleSize("Request"));
  if (request) {
    std::memcpy(native_request.data(), request->bytes, native_request.size());
  }
  void *native_pointer_request = request ? native_request.data() : nullptr;
  NativeBuffer native_status(registry.nativeStatusSize());
  if (status && false) {
    encodeNativeStatus(registry, *status, native_status.data());
  }
  void *native_pointer_status = status ? native_status.data() : nativeStatusIgnore(registry);
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Test", 0x1d207e4cu, 1u, 3u);
  plumed_mpiw_transport_MPI_Test_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Test_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Test");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_pointer_request, flag, native_pointer_status, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (request && valid_output) {
    std::memset(request->bytes, 0, sizeof(request->bytes));
    std::memcpy(request->bytes, native_request.data(), native_request.size());
  }
  if (status && valid_output) {
    decodeNativeStatus(registry, native_status.data(), status);
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Request_free(plumed_mpiw_handle_storage * request) {
  Registry &registry = currentRegistry();
  NativeBuffer native_request(registry.nativeHandleSize("Request"));
  if (request) {
    std::memcpy(native_request.data(), request->bytes, native_request.size());
  }
  void *native_pointer_request = request ? native_request.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Request_free", 0x475f38a1u, 1u, 1u);
  plumed_mpiw_transport_MPI_Request_free_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Request_free_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Request_free");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_pointer_request, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (request && valid_output) {
    std::memset(request->bytes, 0, sizeof(request->bytes));
    std::memcpy(request->bytes, native_request.data(), native_request.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Waitany(int32_t count, plumed_mpiw_handle_storage array_of_requests[], int32_t * indx, plumed_mpiw_status * status) {
  Registry &registry = currentRegistry();
  int native_count = static_cast<int>(count);
  if (count < 0) {
    throw std::invalid_argument("MPI_Waitany requires nonnegative count");
  }
  const std::uint32_t native_array_of_requests_size = registry.nativeHandleSize("Request");
  NativeBuffer native_array_of_requests(static_cast<std::size_t>(count) * native_array_of_requests_size);
  if (array_of_requests && count > 0)
    for (int i = 0; i < count; ++i) {
      std::memcpy(native_array_of_requests.data() + static_cast<std::size_t>(i) * native_array_of_requests_size, array_of_requests[i].bytes, native_array_of_requests_size);
    }
  void *native_pointer_array_of_requests = array_of_requests && count > 0 ? native_array_of_requests.data() : nullptr;
  NativeBuffer native_status(registry.nativeStatusSize());
  if (status && false) {
    encodeNativeStatus(registry, *status, native_status.data());
  }
  void *native_pointer_status = status ? native_status.data() : nativeStatusIgnore(registry);
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Waitany", 0x5acd4436u, 1u, 4u);
  plumed_mpiw_transport_MPI_Waitany_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Waitany_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Waitany");
  }
  std::int32_t native_return_value = 0;
  dispatch(&native_count, native_pointer_array_of_requests, indx, native_pointer_status, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (array_of_requests && valid_output)
    for (int i = 0; i < count; ++i) {
      std::memset(array_of_requests[i].bytes, 0, sizeof(array_of_requests[i].bytes));
      std::memcpy(array_of_requests[i].bytes, native_array_of_requests.data() + static_cast<std::size_t>(i) * native_array_of_requests_size, native_array_of_requests_size);
    }
  if (status && valid_output) {
    decodeNativeStatus(registry, native_status.data(), status);
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Testany(int32_t count, plumed_mpiw_handle_storage array_of_requests[], int32_t * indx, int32_t * flag, plumed_mpiw_status * status) {
  Registry &registry = currentRegistry();
  int native_count = static_cast<int>(count);
  if (count < 0) {
    throw std::invalid_argument("MPI_Testany requires nonnegative count");
  }
  const std::uint32_t native_array_of_requests_size = registry.nativeHandleSize("Request");
  NativeBuffer native_array_of_requests(static_cast<std::size_t>(count) * native_array_of_requests_size);
  if (array_of_requests && count > 0)
    for (int i = 0; i < count; ++i) {
      std::memcpy(native_array_of_requests.data() + static_cast<std::size_t>(i) * native_array_of_requests_size, array_of_requests[i].bytes, native_array_of_requests_size);
    }
  void *native_pointer_array_of_requests = array_of_requests && count > 0 ? native_array_of_requests.data() : nullptr;
  NativeBuffer native_status(registry.nativeStatusSize());
  if (status && false) {
    encodeNativeStatus(registry, *status, native_status.data());
  }
  void *native_pointer_status = status ? native_status.data() : nativeStatusIgnore(registry);
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Testany", 0xf329c686u, 1u, 5u);
  plumed_mpiw_transport_MPI_Testany_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Testany_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Testany");
  }
  std::int32_t native_return_value = 0;
  dispatch(&native_count, native_pointer_array_of_requests, indx, flag, native_pointer_status, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (array_of_requests && valid_output)
    for (int i = 0; i < count; ++i) {
      std::memset(array_of_requests[i].bytes, 0, sizeof(array_of_requests[i].bytes));
      std::memcpy(array_of_requests[i].bytes, native_array_of_requests.data() + static_cast<std::size_t>(i) * native_array_of_requests_size, native_array_of_requests_size);
    }
  if (status && valid_output) {
    decodeNativeStatus(registry, native_status.data(), status);
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Waitall(int32_t count, plumed_mpiw_handle_storage array_of_requests[], plumed_mpiw_status * array_of_statuses) {
  Registry &registry = currentRegistry();
  int native_count = static_cast<int>(count);
  if (count < 0) {
    throw std::invalid_argument("MPI_Waitall requires nonnegative count");
  }
  const std::uint32_t native_array_of_requests_size = registry.nativeHandleSize("Request");
  NativeBuffer native_array_of_requests(static_cast<std::size_t>(count) * native_array_of_requests_size);
  if (array_of_requests && count > 0)
    for (int i = 0; i < count; ++i) {
      std::memcpy(native_array_of_requests.data() + static_cast<std::size_t>(i) * native_array_of_requests_size, array_of_requests[i].bytes, native_array_of_requests_size);
    }
  void *native_pointer_array_of_requests = array_of_requests && count > 0 ? native_array_of_requests.data() : nullptr;
  if (count < 0) {
    throw std::invalid_argument("MPI_Waitall requires nonnegative count");
  }
  const std::uint32_t native_array_of_statuses_size = registry.nativeStatusSize();
  NativeBuffer native_array_of_statuses(static_cast<std::size_t>(count) * native_array_of_statuses_size);
  if (array_of_statuses && count > 0 && false)
    for (int i = 0; i < count; ++i) {
      encodeNativeStatus(registry, array_of_statuses[i], native_array_of_statuses.data() + static_cast<std::size_t>(i) * native_array_of_statuses_size);
    }
  void *native_pointer_array_of_statuses = array_of_statuses && count > 0 ? native_array_of_statuses.data() : nativeStatusesIgnore(registry);
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Waitall", 0x7a57ae23u, 1u, 3u);
  plumed_mpiw_transport_MPI_Waitall_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Waitall_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Waitall");
  }
  std::int32_t native_return_value = 0;
  dispatch(&native_count, native_pointer_array_of_requests, native_pointer_array_of_statuses, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (array_of_requests && valid_output)
    for (int i = 0; i < count; ++i) {
      std::memset(array_of_requests[i].bytes, 0, sizeof(array_of_requests[i].bytes));
      std::memcpy(array_of_requests[i].bytes, native_array_of_requests.data() + static_cast<std::size_t>(i) * native_array_of_requests_size, native_array_of_requests_size);
    }
  if (array_of_statuses && valid_output)
    for (int i = 0; i < count; ++i) {
      decodeNativeStatus(registry, native_array_of_statuses.data() + static_cast<std::size_t>(i) * native_array_of_statuses_size, &array_of_statuses[i]);
    }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Testall(int32_t count, plumed_mpiw_handle_storage array_of_requests[], int32_t * flag, plumed_mpiw_status * array_of_statuses) {
  Registry &registry = currentRegistry();
  int native_count = static_cast<int>(count);
  if (count < 0) {
    throw std::invalid_argument("MPI_Testall requires nonnegative count");
  }
  const std::uint32_t native_array_of_requests_size = registry.nativeHandleSize("Request");
  NativeBuffer native_array_of_requests(static_cast<std::size_t>(count) * native_array_of_requests_size);
  if (array_of_requests && count > 0)
    for (int i = 0; i < count; ++i) {
      std::memcpy(native_array_of_requests.data() + static_cast<std::size_t>(i) * native_array_of_requests_size, array_of_requests[i].bytes, native_array_of_requests_size);
    }
  void *native_pointer_array_of_requests = array_of_requests && count > 0 ? native_array_of_requests.data() : nullptr;
  if (count < 0) {
    throw std::invalid_argument("MPI_Testall requires nonnegative count");
  }
  const std::uint32_t native_array_of_statuses_size = registry.nativeStatusSize();
  NativeBuffer native_array_of_statuses(static_cast<std::size_t>(count) * native_array_of_statuses_size);
  if (array_of_statuses && count > 0 && false)
    for (int i = 0; i < count; ++i) {
      encodeNativeStatus(registry, array_of_statuses[i], native_array_of_statuses.data() + static_cast<std::size_t>(i) * native_array_of_statuses_size);
    }
  void *native_pointer_array_of_statuses = array_of_statuses && count > 0 ? native_array_of_statuses.data() : nativeStatusesIgnore(registry);
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Testall", 0x593a1626u, 1u, 4u);
  plumed_mpiw_transport_MPI_Testall_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Testall_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Testall");
  }
  std::int32_t native_return_value = 0;
  dispatch(&native_count, native_pointer_array_of_requests, flag, native_pointer_array_of_statuses, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (array_of_requests && valid_output)
    for (int i = 0; i < count; ++i) {
      std::memset(array_of_requests[i].bytes, 0, sizeof(array_of_requests[i].bytes));
      std::memcpy(array_of_requests[i].bytes, native_array_of_requests.data() + static_cast<std::size_t>(i) * native_array_of_requests_size, native_array_of_requests_size);
    }
  if (array_of_statuses && valid_output)
    for (int i = 0; i < count; ++i) {
      decodeNativeStatus(registry, native_array_of_statuses.data() + static_cast<std::size_t>(i) * native_array_of_statuses_size, &array_of_statuses[i]);
    }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Waitsome(int32_t incount, plumed_mpiw_handle_storage array_of_requests[], int32_t * outcount, int32_t array_of_indices[], plumed_mpiw_status * array_of_statuses) {
  Registry &registry = currentRegistry();
  int native_incount = static_cast<int>(incount);
  if (incount < 0) {
    throw std::invalid_argument("MPI_Waitsome requires nonnegative incount");
  }
  const std::uint32_t native_array_of_requests_size = registry.nativeHandleSize("Request");
  NativeBuffer native_array_of_requests(static_cast<std::size_t>(incount) * native_array_of_requests_size);
  if (array_of_requests && incount > 0)
    for (int i = 0; i < incount; ++i) {
      std::memcpy(native_array_of_requests.data() + static_cast<std::size_t>(i) * native_array_of_requests_size, array_of_requests[i].bytes, native_array_of_requests_size);
    }
  void *native_pointer_array_of_requests = array_of_requests && incount > 0 ? native_array_of_requests.data() : nullptr;
  if (incount < 0) {
    throw std::invalid_argument("MPI_Waitsome requires nonnegative incount");
  }
  const std::uint32_t native_array_of_statuses_size = registry.nativeStatusSize();
  NativeBuffer native_array_of_statuses(static_cast<std::size_t>(incount) * native_array_of_statuses_size);
  if (array_of_statuses && incount > 0 && false)
    for (int i = 0; i < incount; ++i) {
      encodeNativeStatus(registry, array_of_statuses[i], native_array_of_statuses.data() + static_cast<std::size_t>(i) * native_array_of_statuses_size);
    }
  void *native_pointer_array_of_statuses = array_of_statuses && incount > 0 ? native_array_of_statuses.data() : nativeStatusesIgnore(registry);
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Waitsome", 0xb6562e1au, 1u, 5u);
  plumed_mpiw_transport_MPI_Waitsome_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Waitsome_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Waitsome");
  }
  std::int32_t native_return_value = 0;
  dispatch(&native_incount, native_pointer_array_of_requests, outcount, array_of_indices, native_pointer_array_of_statuses, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (array_of_requests && valid_output)
    for (int i = 0; i < incount; ++i) {
      std::memset(array_of_requests[i].bytes, 0, sizeof(array_of_requests[i].bytes));
      std::memcpy(array_of_requests[i].bytes, native_array_of_requests.data() + static_cast<std::size_t>(i) * native_array_of_requests_size, native_array_of_requests_size);
    }
  if (array_of_statuses && valid_output)
    for (int i = 0; i < (outcount ? std::min(incount, std::max(0, *outcount)) : 0); ++i) {
      decodeNativeStatus(registry, native_array_of_statuses.data() + static_cast<std::size_t>(i) * native_array_of_statuses_size, &array_of_statuses[i]);
    }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Testsome(int32_t incount, plumed_mpiw_handle_storage array_of_requests[], int32_t * outcount, int32_t array_of_indices[], plumed_mpiw_status * array_of_statuses) {
  Registry &registry = currentRegistry();
  int native_incount = static_cast<int>(incount);
  if (incount < 0) {
    throw std::invalid_argument("MPI_Testsome requires nonnegative incount");
  }
  const std::uint32_t native_array_of_requests_size = registry.nativeHandleSize("Request");
  NativeBuffer native_array_of_requests(static_cast<std::size_t>(incount) * native_array_of_requests_size);
  if (array_of_requests && incount > 0)
    for (int i = 0; i < incount; ++i) {
      std::memcpy(native_array_of_requests.data() + static_cast<std::size_t>(i) * native_array_of_requests_size, array_of_requests[i].bytes, native_array_of_requests_size);
    }
  void *native_pointer_array_of_requests = array_of_requests && incount > 0 ? native_array_of_requests.data() : nullptr;
  if (incount < 0) {
    throw std::invalid_argument("MPI_Testsome requires nonnegative incount");
  }
  const std::uint32_t native_array_of_statuses_size = registry.nativeStatusSize();
  NativeBuffer native_array_of_statuses(static_cast<std::size_t>(incount) * native_array_of_statuses_size);
  if (array_of_statuses && incount > 0 && false)
    for (int i = 0; i < incount; ++i) {
      encodeNativeStatus(registry, array_of_statuses[i], native_array_of_statuses.data() + static_cast<std::size_t>(i) * native_array_of_statuses_size);
    }
  void *native_pointer_array_of_statuses = array_of_statuses && incount > 0 ? native_array_of_statuses.data() : nativeStatusesIgnore(registry);
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Testsome", 0x5a87d17bu, 1u, 5u);
  plumed_mpiw_transport_MPI_Testsome_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Testsome_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Testsome");
  }
  std::int32_t native_return_value = 0;
  dispatch(&native_incount, native_pointer_array_of_requests, outcount, array_of_indices, native_pointer_array_of_statuses, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (array_of_requests && valid_output)
    for (int i = 0; i < incount; ++i) {
      std::memset(array_of_requests[i].bytes, 0, sizeof(array_of_requests[i].bytes));
      std::memcpy(array_of_requests[i].bytes, native_array_of_requests.data() + static_cast<std::size_t>(i) * native_array_of_requests_size, native_array_of_requests_size);
    }
  if (array_of_statuses && valid_output)
    for (int i = 0; i < (outcount ? std::min(incount, std::max(0, *outcount)) : 0); ++i) {
      decodeNativeStatus(registry, native_array_of_statuses.data() + static_cast<std::size_t>(i) * native_array_of_statuses_size, &array_of_statuses[i]);
    }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Iprobe(int32_t source, int32_t tag, plumed_mpiw_handle_storage comm, int32_t * flag, plumed_mpiw_status * status) {
  Registry &registry = currentRegistry();
  int native_source = static_cast<int>(source);
  int native_tag = static_cast<int>(tag);
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Iprobe", native_comm.data());
  NativeBuffer native_status(registry.nativeStatusSize());
  if (status && false) {
    encodeNativeStatus(registry, *status, native_status.data());
  }
  void *native_pointer_status = status ? native_status.data() : nativeStatusIgnore(registry);
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Iprobe", 0xd1ba8db8u, 1u, 5u);
  plumed_mpiw_transport_MPI_Iprobe_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Iprobe_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Iprobe");
  }
  std::int32_t native_return_value = 0;
  dispatch(&native_source, &native_tag, native_comm.data(), flag, native_pointer_status, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (status && valid_output) {
    decodeNativeStatus(registry, native_status.data(), status);
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Probe(int32_t source, int32_t tag, plumed_mpiw_handle_storage comm, plumed_mpiw_status * status) {
  Registry &registry = currentRegistry();
  int native_source = static_cast<int>(source);
  int native_tag = static_cast<int>(tag);
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Probe", native_comm.data());
  NativeBuffer native_status(registry.nativeStatusSize());
  if (status && false) {
    encodeNativeStatus(registry, *status, native_status.data());
  }
  void *native_pointer_status = status ? native_status.data() : nativeStatusIgnore(registry);
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Probe", 0xfed81500u, 1u, 4u);
  plumed_mpiw_transport_MPI_Probe_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Probe_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Probe");
  }
  std::int32_t native_return_value = 0;
  dispatch(&native_source, &native_tag, native_comm.data(), native_pointer_status, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (status && valid_output) {
    decodeNativeStatus(registry, native_status.data(), status);
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Cancel(plumed_mpiw_handle_storage * request) {
  Registry &registry = currentRegistry();
  NativeBuffer native_request(registry.nativeHandleSize("Request"));
  if (request) {
    std::memcpy(native_request.data(), request->bytes, native_request.size());
  }
  void *native_pointer_request = request ? native_request.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Cancel", 0xed47f9abu, 1u, 1u);
  plumed_mpiw_transport_MPI_Cancel_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Cancel_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Cancel");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_pointer_request, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (request && valid_output) {
    std::memset(request->bytes, 0, sizeof(request->bytes));
    std::memcpy(request->bytes, native_request.data(), native_request.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Test_cancelled(const plumed_mpiw_status * status, int32_t * flag) {
  Registry &registry = currentRegistry();
  NativeBuffer native_status(registry.nativeStatusSize());
  if (status && true) {
    encodeNativeStatus(registry, *status, native_status.data());
  }
  void *native_pointer_status = status ? native_status.data() : nativeStatusIgnore(registry);
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Test_cancelled", 0xec5e3092u, 1u, 2u);
  plumed_mpiw_transport_MPI_Test_cancelled_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Test_cancelled_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Test_cancelled");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_pointer_status, flag, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Send_init(const void * buf, int32_t count, plumed_mpiw_handle_storage datatype, int32_t dest, int32_t tag, plumed_mpiw_handle_storage comm, plumed_mpiw_handle_storage * request) {
  Registry &registry = currentRegistry();
  int native_count = static_cast<int>(count);
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_datatype.data(), datatype.bytes, native_datatype.size());
  int native_dest = static_cast<int>(dest);
  int native_tag = static_cast<int>(tag);
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Send_init", native_comm.data());
  NativeBuffer native_request(registry.nativeHandleSize("Request"));
  void *native_pointer_request = request ? native_request.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Send_init", 0xd381e137u, 1u, 7u);
  plumed_mpiw_transport_MPI_Send_init_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Send_init_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Send_init");
  }
  std::int32_t native_return_value = 0;
  dispatch(const_cast<void *>(static_cast<const void *>(buf)), &native_count, native_datatype.data(), &native_dest, &native_tag, native_comm.data(), native_pointer_request, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (request && valid_output) {
    std::memset(request->bytes, 0, sizeof(request->bytes));
    std::memcpy(request->bytes, native_request.data(), native_request.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Bsend_init(const void * buf, int32_t count, plumed_mpiw_handle_storage datatype, int32_t dest, int32_t tag, plumed_mpiw_handle_storage comm, plumed_mpiw_handle_storage * request) {
  Registry &registry = currentRegistry();
  int native_count = static_cast<int>(count);
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_datatype.data(), datatype.bytes, native_datatype.size());
  int native_dest = static_cast<int>(dest);
  int native_tag = static_cast<int>(tag);
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Bsend_init", native_comm.data());
  NativeBuffer native_request(registry.nativeHandleSize("Request"));
  void *native_pointer_request = request ? native_request.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Bsend_init", 0x16ede22eu, 1u, 7u);
  plumed_mpiw_transport_MPI_Bsend_init_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Bsend_init_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Bsend_init");
  }
  std::int32_t native_return_value = 0;
  dispatch(const_cast<void *>(static_cast<const void *>(buf)), &native_count, native_datatype.data(), &native_dest, &native_tag, native_comm.data(), native_pointer_request, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (request && valid_output) {
    std::memset(request->bytes, 0, sizeof(request->bytes));
    std::memcpy(request->bytes, native_request.data(), native_request.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Ssend_init(const void * buf, int32_t count, plumed_mpiw_handle_storage datatype, int32_t dest, int32_t tag, plumed_mpiw_handle_storage comm, plumed_mpiw_handle_storage * request) {
  Registry &registry = currentRegistry();
  int native_count = static_cast<int>(count);
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_datatype.data(), datatype.bytes, native_datatype.size());
  int native_dest = static_cast<int>(dest);
  int native_tag = static_cast<int>(tag);
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Ssend_init", native_comm.data());
  NativeBuffer native_request(registry.nativeHandleSize("Request"));
  void *native_pointer_request = request ? native_request.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Ssend_init", 0xeea1c733u, 1u, 7u);
  plumed_mpiw_transport_MPI_Ssend_init_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Ssend_init_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Ssend_init");
  }
  std::int32_t native_return_value = 0;
  dispatch(const_cast<void *>(static_cast<const void *>(buf)), &native_count, native_datatype.data(), &native_dest, &native_tag, native_comm.data(), native_pointer_request, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (request && valid_output) {
    std::memset(request->bytes, 0, sizeof(request->bytes));
    std::memcpy(request->bytes, native_request.data(), native_request.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Rsend_init(const void * buf, int32_t count, plumed_mpiw_handle_storage datatype, int32_t dest, int32_t tag, plumed_mpiw_handle_storage comm, plumed_mpiw_handle_storage * request) {
  Registry &registry = currentRegistry();
  int native_count = static_cast<int>(count);
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_datatype.data(), datatype.bytes, native_datatype.size());
  int native_dest = static_cast<int>(dest);
  int native_tag = static_cast<int>(tag);
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Rsend_init", native_comm.data());
  NativeBuffer native_request(registry.nativeHandleSize("Request"));
  void *native_pointer_request = request ? native_request.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Rsend_init", 0x018fab06u, 1u, 7u);
  plumed_mpiw_transport_MPI_Rsend_init_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Rsend_init_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Rsend_init");
  }
  std::int32_t native_return_value = 0;
  dispatch(const_cast<void *>(static_cast<const void *>(buf)), &native_count, native_datatype.data(), &native_dest, &native_tag, native_comm.data(), native_pointer_request, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (request && valid_output) {
    std::memset(request->bytes, 0, sizeof(request->bytes));
    std::memcpy(request->bytes, native_request.data(), native_request.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Recv_init(void * buf, int32_t count, plumed_mpiw_handle_storage datatype, int32_t source, int32_t tag, plumed_mpiw_handle_storage comm, plumed_mpiw_handle_storage * request) {
  Registry &registry = currentRegistry();
  int native_count = static_cast<int>(count);
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_datatype.data(), datatype.bytes, native_datatype.size());
  int native_source = static_cast<int>(source);
  int native_tag = static_cast<int>(tag);
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Recv_init", native_comm.data());
  NativeBuffer native_request(registry.nativeHandleSize("Request"));
  void *native_pointer_request = request ? native_request.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Recv_init", 0x9cbb4e8au, 1u, 7u);
  plumed_mpiw_transport_MPI_Recv_init_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Recv_init_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Recv_init");
  }
  std::int32_t native_return_value = 0;
  dispatch(buf, &native_count, native_datatype.data(), &native_source, &native_tag, native_comm.data(), native_pointer_request, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (request && valid_output) {
    std::memset(request->bytes, 0, sizeof(request->bytes));
    std::memcpy(request->bytes, native_request.data(), native_request.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Start(plumed_mpiw_handle_storage * request) {
  Registry &registry = currentRegistry();
  NativeBuffer native_request(registry.nativeHandleSize("Request"));
  if (request) {
    std::memcpy(native_request.data(), request->bytes, native_request.size());
  }
  void *native_pointer_request = request ? native_request.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Start", 0xfb62fbe3u, 1u, 1u);
  plumed_mpiw_transport_MPI_Start_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Start_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Start");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_pointer_request, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (request && valid_output) {
    std::memset(request->bytes, 0, sizeof(request->bytes));
    std::memcpy(request->bytes, native_request.data(), native_request.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Startall(int32_t count, plumed_mpiw_handle_storage array_of_requests[]) {
  Registry &registry = currentRegistry();
  int native_count = static_cast<int>(count);
  if (count < 0) {
    throw std::invalid_argument("MPI_Startall requires nonnegative count");
  }
  const std::uint32_t native_array_of_requests_size = registry.nativeHandleSize("Request");
  NativeBuffer native_array_of_requests(static_cast<std::size_t>(count) * native_array_of_requests_size);
  if (array_of_requests && count > 0)
    for (int i = 0; i < count; ++i) {
      std::memcpy(native_array_of_requests.data() + static_cast<std::size_t>(i) * native_array_of_requests_size, array_of_requests[i].bytes, native_array_of_requests_size);
    }
  void *native_pointer_array_of_requests = array_of_requests && count > 0 ? native_array_of_requests.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Startall", 0xc3d7cddfu, 1u, 2u);
  plumed_mpiw_transport_MPI_Startall_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Startall_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Startall");
  }
  std::int32_t native_return_value = 0;
  dispatch(&native_count, native_pointer_array_of_requests, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (array_of_requests && valid_output)
    for (int i = 0; i < count; ++i) {
      std::memset(array_of_requests[i].bytes, 0, sizeof(array_of_requests[i].bytes));
      std::memcpy(array_of_requests[i].bytes, native_array_of_requests.data() + static_cast<std::size_t>(i) * native_array_of_requests_size, native_array_of_requests_size);
    }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Sendrecv_replace(void * buf, int32_t count, plumed_mpiw_handle_storage datatype, int32_t dest, int32_t sendtag, int32_t source, int32_t recvtag, plumed_mpiw_handle_storage comm, plumed_mpiw_status * status) {
  Registry &registry = currentRegistry();
  int native_count = static_cast<int>(count);
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_datatype.data(), datatype.bytes, native_datatype.size());
  int native_dest = static_cast<int>(dest);
  int native_sendtag = static_cast<int>(sendtag);
  int native_source = static_cast<int>(source);
  int native_recvtag = static_cast<int>(recvtag);
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Sendrecv_replace", native_comm.data());
  NativeBuffer native_status(registry.nativeStatusSize());
  if (status && false) {
    encodeNativeStatus(registry, *status, native_status.data());
  }
  void *native_pointer_status = status ? native_status.data() : nativeStatusIgnore(registry);
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Sendrecv_replace", 0x7c6ffc57u, 1u, 9u);
  plumed_mpiw_transport_MPI_Sendrecv_replace_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Sendrecv_replace_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Sendrecv_replace");
  }
  std::int32_t native_return_value = 0;
  dispatch(buf, &native_count, native_datatype.data(), &native_dest, &native_sendtag, &native_source, &native_recvtag, native_comm.data(), native_pointer_status, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (status && valid_output) {
    decodeNativeStatus(registry, native_status.data(), status);
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Type_contiguous(int32_t count, plumed_mpiw_handle_storage oldtype, plumed_mpiw_handle_storage * newtype) {
  Registry &registry = currentRegistry();
  int native_count = static_cast<int>(count);
  NativeBuffer native_oldtype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_oldtype.data(), oldtype.bytes, native_oldtype.size());
  NativeBuffer native_newtype(registry.nativeHandleSize("Datatype"));
  void *native_pointer_newtype = newtype ? native_newtype.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Type_contiguous", 0xe41681d2u, 1u, 3u);
  plumed_mpiw_transport_MPI_Type_contiguous_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Type_contiguous_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Type_contiguous");
  }
  std::int32_t native_return_value = 0;
  dispatch(&native_count, native_oldtype.data(), native_pointer_newtype, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (newtype && valid_output) {
    std::memset(newtype->bytes, 0, sizeof(newtype->bytes));
    std::memcpy(newtype->bytes, native_newtype.data(), native_newtype.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Type_vector(int32_t count, int32_t blocklength, int32_t stride, plumed_mpiw_handle_storage oldtype, plumed_mpiw_handle_storage * newtype) {
  Registry &registry = currentRegistry();
  int native_count = static_cast<int>(count);
  int native_blocklength = static_cast<int>(blocklength);
  int native_stride = static_cast<int>(stride);
  NativeBuffer native_oldtype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_oldtype.data(), oldtype.bytes, native_oldtype.size());
  NativeBuffer native_newtype(registry.nativeHandleSize("Datatype"));
  void *native_pointer_newtype = newtype ? native_newtype.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Type_vector", 0xc7b9dbb0u, 1u, 5u);
  plumed_mpiw_transport_MPI_Type_vector_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Type_vector_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Type_vector");
  }
  std::int32_t native_return_value = 0;
  dispatch(&native_count, &native_blocklength, &native_stride, native_oldtype.data(), native_pointer_newtype, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (newtype && valid_output) {
    std::memset(newtype->bytes, 0, sizeof(newtype->bytes));
    std::memcpy(newtype->bytes, native_newtype.data(), native_newtype.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Type_indexed(int32_t count, const int32_t array_of_blocklengths[], const int32_t array_of_displacements[], plumed_mpiw_handle_storage oldtype, plumed_mpiw_handle_storage * newtype) {
  Registry &registry = currentRegistry();
  int native_count = static_cast<int>(count);
  NativeBuffer native_oldtype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_oldtype.data(), oldtype.bytes, native_oldtype.size());
  NativeBuffer native_newtype(registry.nativeHandleSize("Datatype"));
  void *native_pointer_newtype = newtype ? native_newtype.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Type_indexed", 0xeb606d34u, 1u, 5u);
  plumed_mpiw_transport_MPI_Type_indexed_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Type_indexed_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Type_indexed");
  }
  std::int32_t native_return_value = 0;
  dispatch(&native_count, const_cast<void *>(static_cast<const void *>(array_of_blocklengths)), const_cast<void *>(static_cast<const void *>(array_of_displacements)), native_oldtype.data(), native_pointer_newtype, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (newtype && valid_output) {
    std::memset(newtype->bytes, 0, sizeof(newtype->bytes));
    std::memcpy(newtype->bytes, native_newtype.data(), native_newtype.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Type_commit(plumed_mpiw_handle_storage * datatype) {
  Registry &registry = currentRegistry();
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  if (datatype) {
    std::memcpy(native_datatype.data(), datatype->bytes, native_datatype.size());
  }
  void *native_pointer_datatype = datatype ? native_datatype.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Type_commit", 0xf848da1fu, 1u, 1u);
  plumed_mpiw_transport_MPI_Type_commit_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Type_commit_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Type_commit");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_pointer_datatype, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (datatype && valid_output) {
    std::memset(datatype->bytes, 0, sizeof(datatype->bytes));
    std::memcpy(datatype->bytes, native_datatype.data(), native_datatype.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Type_free(plumed_mpiw_handle_storage * datatype) {
  Registry &registry = currentRegistry();
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  if (datatype) {
    std::memcpy(native_datatype.data(), datatype->bytes, native_datatype.size());
  }
  void *native_pointer_datatype = datatype ? native_datatype.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Type_free", 0xdacde819u, 1u, 1u);
  plumed_mpiw_transport_MPI_Type_free_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Type_free_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Type_free");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_pointer_datatype, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (datatype && valid_output) {
    std::memset(datatype->bytes, 0, sizeof(datatype->bytes));
    std::memcpy(datatype->bytes, native_datatype.data(), native_datatype.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Type_size(plumed_mpiw_handle_storage datatype, int32_t * size) {
  Registry &registry = currentRegistry();
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_datatype.data(), datatype.bytes, native_datatype.size());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Type_size", 0x95660c77u, 1u, 2u);
  plumed_mpiw_transport_MPI_Type_size_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Type_size_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Type_size");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_datatype.data(), size, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Get_elements(const plumed_mpiw_status * status, plumed_mpiw_handle_storage datatype, int32_t * count) {
  Registry &registry = currentRegistry();
  NativeBuffer native_status(registry.nativeStatusSize());
  if (status && true) {
    encodeNativeStatus(registry, *status, native_status.data());
  }
  void *native_pointer_status = status ? native_status.data() : nativeStatusIgnore(registry);
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_datatype.data(), datatype.bytes, native_datatype.size());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Get_elements", 0xf71eaa9au, 1u, 3u);
  plumed_mpiw_transport_MPI_Get_elements_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Get_elements_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Get_elements");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_pointer_status, native_datatype.data(), count, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Pack(const void * inbuf, int32_t incount, plumed_mpiw_handle_storage datatype, void * outbuf, int32_t outsize, int32_t * position, plumed_mpiw_handle_storage comm) {
  Registry &registry = currentRegistry();
  int native_incount = static_cast<int>(incount);
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_datatype.data(), datatype.bytes, native_datatype.size());
  int native_outsize = static_cast<int>(outsize);
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Pack", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Pack", 0x5649ba7au, 1u, 7u);
  plumed_mpiw_transport_MPI_Pack_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Pack_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Pack");
  }
  std::int32_t native_return_value = 0;
  dispatch(const_cast<void *>(static_cast<const void *>(inbuf)), &native_incount, native_datatype.data(), outbuf, &native_outsize, position, native_comm.data(), &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Unpack(const void * inbuf, int32_t insize, int32_t * position, void * outbuf, int32_t outcount, plumed_mpiw_handle_storage datatype, plumed_mpiw_handle_storage comm) {
  Registry &registry = currentRegistry();
  int native_insize = static_cast<int>(insize);
  int native_outcount = static_cast<int>(outcount);
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_datatype.data(), datatype.bytes, native_datatype.size());
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Unpack", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Unpack", 0x49901f8bu, 1u, 7u);
  plumed_mpiw_transport_MPI_Unpack_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Unpack_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Unpack");
  }
  std::int32_t native_return_value = 0;
  dispatch(const_cast<void *>(static_cast<const void *>(inbuf)), &native_insize, position, outbuf, &native_outcount, native_datatype.data(), native_comm.data(), &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Pack_size(int32_t incount, plumed_mpiw_handle_storage datatype, plumed_mpiw_handle_storage comm, int32_t * size) {
  Registry &registry = currentRegistry();
  int native_incount = static_cast<int>(incount);
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_datatype.data(), datatype.bytes, native_datatype.size());
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Pack_size", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Pack_size", 0x24de0131u, 1u, 4u);
  plumed_mpiw_transport_MPI_Pack_size_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Pack_size_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Pack_size");
  }
  std::int32_t native_return_value = 0;
  dispatch(&native_incount, native_datatype.data(), native_comm.data(), size, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Gatherv(const void * sendbuf, int32_t sendcount, plumed_mpiw_handle_storage sendtype, void * recvbuf, const int32_t recvcounts[], const int32_t displs[], plumed_mpiw_handle_storage recvtype, int32_t root, plumed_mpiw_handle_storage comm) {
  Registry &registry = currentRegistry();
  int native_sendcount = static_cast<int>(sendcount);
  NativeBuffer native_sendtype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_sendtype.data(), sendtype.bytes, native_sendtype.size());
  NativeBuffer native_recvtype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_recvtype.data(), recvtype.bytes, native_recvtype.size());
  int native_root = static_cast<int>(root);
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Gatherv", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Gatherv", 0xc90455c9u, 1u, 9u);
  plumed_mpiw_transport_MPI_Gatherv_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Gatherv_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Gatherv");
  }
  std::int32_t native_return_value = 0;
  dispatch(const_cast<void *>(static_cast<const void *>(sendbuf)), &native_sendcount, native_sendtype.data(), recvbuf, const_cast<void *>(static_cast<const void *>(recvcounts)), const_cast<void *>(static_cast<const void *>(displs)), native_recvtype.data(), &native_root, native_comm.data(), &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Scatterv(const void * sendbuf, const int32_t sendcounts[], const int32_t displs[], plumed_mpiw_handle_storage sendtype, void * recvbuf, int32_t recvcount, plumed_mpiw_handle_storage recvtype, int32_t root, plumed_mpiw_handle_storage comm) {
  Registry &registry = currentRegistry();
  NativeBuffer native_sendtype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_sendtype.data(), sendtype.bytes, native_sendtype.size());
  int native_recvcount = static_cast<int>(recvcount);
  NativeBuffer native_recvtype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_recvtype.data(), recvtype.bytes, native_recvtype.size());
  int native_root = static_cast<int>(root);
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Scatterv", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Scatterv", 0x765eb673u, 1u, 9u);
  plumed_mpiw_transport_MPI_Scatterv_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Scatterv_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Scatterv");
  }
  std::int32_t native_return_value = 0;
  dispatch(const_cast<void *>(static_cast<const void *>(sendbuf)), const_cast<void *>(static_cast<const void *>(sendcounts)), const_cast<void *>(static_cast<const void *>(displs)), native_sendtype.data(), recvbuf, &native_recvcount, native_recvtype.data(), &native_root, native_comm.data(), &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Alltoall(const void * sendbuf, int32_t sendcount, plumed_mpiw_handle_storage sendtype, void * recvbuf, int32_t recvcount, plumed_mpiw_handle_storage recvtype, plumed_mpiw_handle_storage comm) {
  Registry &registry = currentRegistry();
  int native_sendcount = static_cast<int>(sendcount);
  NativeBuffer native_sendtype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_sendtype.data(), sendtype.bytes, native_sendtype.size());
  int native_recvcount = static_cast<int>(recvcount);
  NativeBuffer native_recvtype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_recvtype.data(), recvtype.bytes, native_recvtype.size());
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Alltoall", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Alltoall", 0xcfe270aeu, 1u, 7u);
  plumed_mpiw_transport_MPI_Alltoall_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Alltoall_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Alltoall");
  }
  std::int32_t native_return_value = 0;
  dispatch(const_cast<void *>(static_cast<const void *>(sendbuf)), &native_sendcount, native_sendtype.data(), recvbuf, &native_recvcount, native_recvtype.data(), native_comm.data(), &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Alltoallv(const void * sendbuf, const int32_t sendcounts[], const int32_t sdispls[], plumed_mpiw_handle_storage sendtype, void * recvbuf, const int32_t recvcounts[], const int32_t rdispls[], plumed_mpiw_handle_storage recvtype, plumed_mpiw_handle_storage comm) {
  Registry &registry = currentRegistry();
  NativeBuffer native_sendtype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_sendtype.data(), sendtype.bytes, native_sendtype.size());
  NativeBuffer native_recvtype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_recvtype.data(), recvtype.bytes, native_recvtype.size());
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Alltoallv", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Alltoallv", 0xa4d99ed0u, 1u, 9u);
  plumed_mpiw_transport_MPI_Alltoallv_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Alltoallv_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Alltoallv");
  }
  std::int32_t native_return_value = 0;
  dispatch(const_cast<void *>(static_cast<const void *>(sendbuf)), const_cast<void *>(static_cast<const void *>(sendcounts)), const_cast<void *>(static_cast<const void *>(sdispls)), native_sendtype.data(), recvbuf, const_cast<void *>(static_cast<const void *>(recvcounts)), const_cast<void *>(static_cast<const void *>(rdispls)), native_recvtype.data(), native_comm.data(), &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Reduce(const void * sendbuf, void * recvbuf, int32_t count, plumed_mpiw_handle_storage datatype, plumed_mpiw_handle_storage op, int32_t root, plumed_mpiw_handle_storage comm) {
  Registry &registry = currentRegistry();
  int native_count = static_cast<int>(count);
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_datatype.data(), datatype.bytes, native_datatype.size());
  NativeBuffer native_op(registry.nativeHandleSize("Op"));
  std::memcpy(native_op.data(), op.bytes, native_op.size());
  int native_root = static_cast<int>(root);
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Reduce", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Reduce", 0x3a815604u, 1u, 7u);
  plumed_mpiw_transport_MPI_Reduce_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Reduce_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Reduce");
  }
  std::int32_t native_return_value = 0;
  dispatch(const_cast<void *>(static_cast<const void *>(sendbuf)), recvbuf, &native_count, native_datatype.data(), native_op.data(), &native_root, native_comm.data(), &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Op_free(plumed_mpiw_handle_storage * op) {
  Registry &registry = currentRegistry();
  NativeBuffer native_op(registry.nativeHandleSize("Op"));
  if (op) {
    std::memcpy(native_op.data(), op->bytes, native_op.size());
  }
  void *native_pointer_op = op ? native_op.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Op_free", 0x94fea4e6u, 1u, 1u);
  plumed_mpiw_transport_MPI_Op_free_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Op_free_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Op_free");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_pointer_op, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (op && valid_output) {
    std::memset(op->bytes, 0, sizeof(op->bytes));
    std::memcpy(op->bytes, native_op.data(), native_op.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Reduce_scatter(const void * sendbuf, void * recvbuf, const int32_t recvcounts[], plumed_mpiw_handle_storage datatype, plumed_mpiw_handle_storage op, plumed_mpiw_handle_storage comm) {
  Registry &registry = currentRegistry();
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_datatype.data(), datatype.bytes, native_datatype.size());
  NativeBuffer native_op(registry.nativeHandleSize("Op"));
  std::memcpy(native_op.data(), op.bytes, native_op.size());
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Reduce_scatter", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Reduce_scatter", 0x2f70a384u, 1u, 6u);
  plumed_mpiw_transport_MPI_Reduce_scatter_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Reduce_scatter_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Reduce_scatter");
  }
  std::int32_t native_return_value = 0;
  dispatch(const_cast<void *>(static_cast<const void *>(sendbuf)), recvbuf, const_cast<void *>(static_cast<const void *>(recvcounts)), native_datatype.data(), native_op.data(), native_comm.data(), &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Scan(const void * sendbuf, void * recvbuf, int32_t count, plumed_mpiw_handle_storage datatype, plumed_mpiw_handle_storage op, plumed_mpiw_handle_storage comm) {
  Registry &registry = currentRegistry();
  int native_count = static_cast<int>(count);
  NativeBuffer native_datatype(registry.nativeHandleSize("Datatype"));
  std::memcpy(native_datatype.data(), datatype.bytes, native_datatype.size());
  NativeBuffer native_op(registry.nativeHandleSize("Op"));
  std::memcpy(native_op.data(), op.bytes, native_op.size());
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Scan", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Scan", 0x64e2ee54u, 1u, 6u);
  plumed_mpiw_transport_MPI_Scan_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Scan_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Scan");
  }
  std::int32_t native_return_value = 0;
  dispatch(const_cast<void *>(static_cast<const void *>(sendbuf)), recvbuf, &native_count, native_datatype.data(), native_op.data(), native_comm.data(), &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Group_size(plumed_mpiw_handle_storage group, int32_t * size) {
  Registry &registry = currentRegistry();
  NativeBuffer native_group(registry.nativeHandleSize("Group"));
  std::memcpy(native_group.data(), group.bytes, native_group.size());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Group_size", 0x7de536d5u, 1u, 2u);
  plumed_mpiw_transport_MPI_Group_size_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Group_size_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Group_size");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_group.data(), size, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Group_rank(plumed_mpiw_handle_storage group, int32_t * rank) {
  Registry &registry = currentRegistry();
  NativeBuffer native_group(registry.nativeHandleSize("Group"));
  std::memcpy(native_group.data(), group.bytes, native_group.size());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Group_rank", 0x668fab96u, 1u, 2u);
  plumed_mpiw_transport_MPI_Group_rank_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Group_rank_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Group_rank");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_group.data(), rank, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Group_translate_ranks(plumed_mpiw_handle_storage group1, int32_t n, const int32_t ranks1[], plumed_mpiw_handle_storage group2, int32_t ranks2[]) {
  Registry &registry = currentRegistry();
  NativeBuffer native_group1(registry.nativeHandleSize("Group"));
  std::memcpy(native_group1.data(), group1.bytes, native_group1.size());
  int native_n = static_cast<int>(n);
  NativeBuffer native_group2(registry.nativeHandleSize("Group"));
  std::memcpy(native_group2.data(), group2.bytes, native_group2.size());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Group_translate_ranks", 0x2ac4d2c6u, 1u, 5u);
  plumed_mpiw_transport_MPI_Group_translate_ranks_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Group_translate_ranks_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Group_translate_ranks");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_group1.data(), &native_n, const_cast<void *>(static_cast<const void *>(ranks1)), native_group2.data(), ranks2, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Group_compare(plumed_mpiw_handle_storage group1, plumed_mpiw_handle_storage group2, int32_t * result) {
  Registry &registry = currentRegistry();
  NativeBuffer native_group1(registry.nativeHandleSize("Group"));
  std::memcpy(native_group1.data(), group1.bytes, native_group1.size());
  NativeBuffer native_group2(registry.nativeHandleSize("Group"));
  std::memcpy(native_group2.data(), group2.bytes, native_group2.size());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Group_compare", 0xb0d34b44u, 1u, 3u);
  plumed_mpiw_transport_MPI_Group_compare_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Group_compare_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Group_compare");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_group1.data(), native_group2.data(), result, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Comm_group(plumed_mpiw_handle_storage comm, plumed_mpiw_handle_storage * group) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Comm_group", native_comm.data());
  NativeBuffer native_group(registry.nativeHandleSize("Group"));
  void *native_pointer_group = group ? native_group.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Comm_group", 0xd11f4483u, 1u, 2u);
  plumed_mpiw_transport_MPI_Comm_group_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Comm_group_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Comm_group");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm.data(), native_pointer_group, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (group && valid_output) {
    std::memset(group->bytes, 0, sizeof(group->bytes));
    std::memcpy(group->bytes, native_group.data(), native_group.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Group_union(plumed_mpiw_handle_storage group1, plumed_mpiw_handle_storage group2, plumed_mpiw_handle_storage * newgroup) {
  Registry &registry = currentRegistry();
  NativeBuffer native_group1(registry.nativeHandleSize("Group"));
  std::memcpy(native_group1.data(), group1.bytes, native_group1.size());
  NativeBuffer native_group2(registry.nativeHandleSize("Group"));
  std::memcpy(native_group2.data(), group2.bytes, native_group2.size());
  NativeBuffer native_newgroup(registry.nativeHandleSize("Group"));
  void *native_pointer_newgroup = newgroup ? native_newgroup.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Group_union", 0xc288bd65u, 1u, 3u);
  plumed_mpiw_transport_MPI_Group_union_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Group_union_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Group_union");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_group1.data(), native_group2.data(), native_pointer_newgroup, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (newgroup && valid_output) {
    std::memset(newgroup->bytes, 0, sizeof(newgroup->bytes));
    std::memcpy(newgroup->bytes, native_newgroup.data(), native_newgroup.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Group_intersection(plumed_mpiw_handle_storage group1, plumed_mpiw_handle_storage group2, plumed_mpiw_handle_storage * newgroup) {
  Registry &registry = currentRegistry();
  NativeBuffer native_group1(registry.nativeHandleSize("Group"));
  std::memcpy(native_group1.data(), group1.bytes, native_group1.size());
  NativeBuffer native_group2(registry.nativeHandleSize("Group"));
  std::memcpy(native_group2.data(), group2.bytes, native_group2.size());
  NativeBuffer native_newgroup(registry.nativeHandleSize("Group"));
  void *native_pointer_newgroup = newgroup ? native_newgroup.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Group_intersection", 0xddf0ecb4u, 1u, 3u);
  plumed_mpiw_transport_MPI_Group_intersection_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Group_intersection_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Group_intersection");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_group1.data(), native_group2.data(), native_pointer_newgroup, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (newgroup && valid_output) {
    std::memset(newgroup->bytes, 0, sizeof(newgroup->bytes));
    std::memcpy(newgroup->bytes, native_newgroup.data(), native_newgroup.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Group_difference(plumed_mpiw_handle_storage group1, plumed_mpiw_handle_storage group2, plumed_mpiw_handle_storage * newgroup) {
  Registry &registry = currentRegistry();
  NativeBuffer native_group1(registry.nativeHandleSize("Group"));
  std::memcpy(native_group1.data(), group1.bytes, native_group1.size());
  NativeBuffer native_group2(registry.nativeHandleSize("Group"));
  std::memcpy(native_group2.data(), group2.bytes, native_group2.size());
  NativeBuffer native_newgroup(registry.nativeHandleSize("Group"));
  void *native_pointer_newgroup = newgroup ? native_newgroup.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Group_difference", 0x37101ec1u, 1u, 3u);
  plumed_mpiw_transport_MPI_Group_difference_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Group_difference_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Group_difference");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_group1.data(), native_group2.data(), native_pointer_newgroup, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (newgroup && valid_output) {
    std::memset(newgroup->bytes, 0, sizeof(newgroup->bytes));
    std::memcpy(newgroup->bytes, native_newgroup.data(), native_newgroup.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Group_incl(plumed_mpiw_handle_storage group, int32_t n, const int32_t ranks[], plumed_mpiw_handle_storage * newgroup) {
  Registry &registry = currentRegistry();
  NativeBuffer native_group(registry.nativeHandleSize("Group"));
  std::memcpy(native_group.data(), group.bytes, native_group.size());
  int native_n = static_cast<int>(n);
  NativeBuffer native_newgroup(registry.nativeHandleSize("Group"));
  void *native_pointer_newgroup = newgroup ? native_newgroup.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Group_incl", 0xbe8bfe71u, 1u, 4u);
  plumed_mpiw_transport_MPI_Group_incl_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Group_incl_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Group_incl");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_group.data(), &native_n, const_cast<void *>(static_cast<const void *>(ranks)), native_pointer_newgroup, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (newgroup && valid_output) {
    std::memset(newgroup->bytes, 0, sizeof(newgroup->bytes));
    std::memcpy(newgroup->bytes, native_newgroup.data(), native_newgroup.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Group_excl(plumed_mpiw_handle_storage group, int32_t n, const int32_t ranks[], plumed_mpiw_handle_storage * newgroup) {
  Registry &registry = currentRegistry();
  NativeBuffer native_group(registry.nativeHandleSize("Group"));
  std::memcpy(native_group.data(), group.bytes, native_group.size());
  int native_n = static_cast<int>(n);
  NativeBuffer native_newgroup(registry.nativeHandleSize("Group"));
  void *native_pointer_newgroup = newgroup ? native_newgroup.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Group_excl", 0x45fa0902u, 1u, 4u);
  plumed_mpiw_transport_MPI_Group_excl_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Group_excl_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Group_excl");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_group.data(), &native_n, const_cast<void *>(static_cast<const void *>(ranks)), native_pointer_newgroup, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (newgroup && valid_output) {
    std::memset(newgroup->bytes, 0, sizeof(newgroup->bytes));
    std::memcpy(newgroup->bytes, native_newgroup.data(), native_newgroup.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Group_range_incl(plumed_mpiw_handle_storage group, int32_t n, int32_t ranges[][3], plumed_mpiw_handle_storage * newgroup) {
  Registry &registry = currentRegistry();
  NativeBuffer native_group(registry.nativeHandleSize("Group"));
  std::memcpy(native_group.data(), group.bytes, native_group.size());
  int native_n = static_cast<int>(n);
  NativeBuffer native_newgroup(registry.nativeHandleSize("Group"));
  void *native_pointer_newgroup = newgroup ? native_newgroup.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Group_range_incl", 0x801adc54u, 1u, 4u);
  plumed_mpiw_transport_MPI_Group_range_incl_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Group_range_incl_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Group_range_incl");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_group.data(), &native_n, ranges, native_pointer_newgroup, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (newgroup && valid_output) {
    std::memset(newgroup->bytes, 0, sizeof(newgroup->bytes));
    std::memcpy(newgroup->bytes, native_newgroup.data(), native_newgroup.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Group_range_excl(plumed_mpiw_handle_storage group, int32_t n, int32_t ranges[][3], plumed_mpiw_handle_storage * newgroup) {
  Registry &registry = currentRegistry();
  NativeBuffer native_group(registry.nativeHandleSize("Group"));
  std::memcpy(native_group.data(), group.bytes, native_group.size());
  int native_n = static_cast<int>(n);
  NativeBuffer native_newgroup(registry.nativeHandleSize("Group"));
  void *native_pointer_newgroup = newgroup ? native_newgroup.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Group_range_excl", 0xc416deb6u, 1u, 4u);
  plumed_mpiw_transport_MPI_Group_range_excl_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Group_range_excl_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Group_range_excl");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_group.data(), &native_n, ranges, native_pointer_newgroup, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (newgroup && valid_output) {
    std::memset(newgroup->bytes, 0, sizeof(newgroup->bytes));
    std::memcpy(newgroup->bytes, native_newgroup.data(), native_newgroup.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Group_free(plumed_mpiw_handle_storage * group) {
  Registry &registry = currentRegistry();
  NativeBuffer native_group(registry.nativeHandleSize("Group"));
  if (group) {
    std::memcpy(native_group.data(), group->bytes, native_group.size());
  }
  void *native_pointer_group = group ? native_group.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Group_free", 0xfd1c62d1u, 1u, 1u);
  plumed_mpiw_transport_MPI_Group_free_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Group_free_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Group_free");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_pointer_group, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (group && valid_output) {
    std::memset(group->bytes, 0, sizeof(group->bytes));
    std::memcpy(group->bytes, native_group.data(), native_group.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Comm_compare(plumed_mpiw_handle_storage comm1, plumed_mpiw_handle_storage comm2, int32_t * result) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm1(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm1.data(), comm1.bytes, native_comm1.size());
  validateInputComm(registry, "MPI_Comm_compare", native_comm1.data());
  NativeBuffer native_comm2(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm2.data(), comm2.bytes, native_comm2.size());
  validateInputComm(registry, "MPI_Comm_compare", native_comm2.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Comm_compare", 0xdacd79cdu, 1u, 3u);
  plumed_mpiw_transport_MPI_Comm_compare_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Comm_compare_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Comm_compare");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm1.data(), native_comm2.data(), result, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Comm_create(plumed_mpiw_handle_storage comm, plumed_mpiw_handle_storage group, plumed_mpiw_handle_storage * newcomm) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Comm_create", native_comm.data());
  NativeBuffer native_group(registry.nativeHandleSize("Group"));
  std::memcpy(native_group.data(), group.bytes, native_group.size());
  NativeBuffer native_newcomm(registry.nativeHandleSize("Comm"));
  void *native_pointer_newcomm = newcomm ? native_newcomm.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Comm_create", 0x1de85fa0u, 1u, 3u);
  plumed_mpiw_transport_MPI_Comm_create_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Comm_create_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Comm_create");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm.data(), native_group.data(), native_pointer_newcomm, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (newcomm && valid_output) {
    std::memset(newcomm->bytes, 0, sizeof(newcomm->bytes));
    std::memcpy(newcomm->bytes, native_newcomm.data(), native_newcomm.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Comm_test_inter(plumed_mpiw_handle_storage comm, int32_t * flag) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Comm_test_inter", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Comm_test_inter", 0x29431089u, 1u, 2u);
  plumed_mpiw_transport_MPI_Comm_test_inter_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Comm_test_inter_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Comm_test_inter");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm.data(), flag, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Comm_remote_size(plumed_mpiw_handle_storage comm, int32_t * size) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Comm_remote_size", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Comm_remote_size", 0xaa0c6ac3u, 1u, 2u);
  plumed_mpiw_transport_MPI_Comm_remote_size_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Comm_remote_size_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Comm_remote_size");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm.data(), size, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Comm_remote_group(plumed_mpiw_handle_storage comm, plumed_mpiw_handle_storage * group) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Comm_remote_group", native_comm.data());
  NativeBuffer native_group(registry.nativeHandleSize("Group"));
  void *native_pointer_group = group ? native_group.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Comm_remote_group", 0x49ee9d8fu, 1u, 2u);
  plumed_mpiw_transport_MPI_Comm_remote_group_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Comm_remote_group_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Comm_remote_group");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm.data(), native_pointer_group, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (group && valid_output) {
    std::memset(group->bytes, 0, sizeof(group->bytes));
    std::memcpy(group->bytes, native_group.data(), native_group.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Intercomm_create(plumed_mpiw_handle_storage local_comm, int32_t local_leader, plumed_mpiw_handle_storage peer_comm, int32_t remote_leader, int32_t tag, plumed_mpiw_handle_storage * newintercomm) {
  Registry &registry = currentRegistry();
  NativeBuffer native_local_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_local_comm.data(), local_comm.bytes, native_local_comm.size());
  validateInputComm(registry, "MPI_Intercomm_create", native_local_comm.data());
  int native_local_leader = static_cast<int>(local_leader);
  NativeBuffer native_peer_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_peer_comm.data(), peer_comm.bytes, native_peer_comm.size());
  validateInputComm(registry, "MPI_Intercomm_create", native_peer_comm.data());
  int native_remote_leader = static_cast<int>(remote_leader);
  int native_tag = static_cast<int>(tag);
  NativeBuffer native_newintercomm(registry.nativeHandleSize("Comm"));
  void *native_pointer_newintercomm = newintercomm ? native_newintercomm.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Intercomm_create", 0x320542c0u, 1u, 6u);
  plumed_mpiw_transport_MPI_Intercomm_create_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Intercomm_create_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Intercomm_create");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_local_comm.data(), &native_local_leader, native_peer_comm.data(), &native_remote_leader, &native_tag, native_pointer_newintercomm, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (newintercomm && valid_output) {
    std::memset(newintercomm->bytes, 0, sizeof(newintercomm->bytes));
    std::memcpy(newintercomm->bytes, native_newintercomm.data(), native_newintercomm.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Intercomm_merge(plumed_mpiw_handle_storage intercomm, int32_t high, plumed_mpiw_handle_storage * newintracomm) {
  Registry &registry = currentRegistry();
  NativeBuffer native_intercomm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_intercomm.data(), intercomm.bytes, native_intercomm.size());
  validateInputComm(registry, "MPI_Intercomm_merge", native_intercomm.data());
  int native_high = static_cast<int>(high);
  NativeBuffer native_newintracomm(registry.nativeHandleSize("Comm"));
  void *native_pointer_newintracomm = newintracomm ? native_newintracomm.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Intercomm_merge", 0x1f6a920au, 1u, 3u);
  plumed_mpiw_transport_MPI_Intercomm_merge_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Intercomm_merge_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Intercomm_merge");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_intercomm.data(), &native_high, native_pointer_newintracomm, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (newintracomm && valid_output) {
    std::memset(newintracomm->bytes, 0, sizeof(newintracomm->bytes));
    std::memcpy(newintracomm->bytes, native_newintracomm.data(), native_newintracomm.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Cart_create(plumed_mpiw_handle_storage comm_old, int32_t ndims, const int32_t dims[], const int32_t periods[], int32_t reorder, plumed_mpiw_handle_storage * comm_cart) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm_old(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm_old.data(), comm_old.bytes, native_comm_old.size());
  validateInputComm(registry, "MPI_Cart_create", native_comm_old.data());
  int native_ndims = static_cast<int>(ndims);
  int native_reorder = static_cast<int>(reorder);
  NativeBuffer native_comm_cart(registry.nativeHandleSize("Comm"));
  void *native_pointer_comm_cart = comm_cart ? native_comm_cart.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Cart_create", 0xfd363844u, 1u, 6u);
  plumed_mpiw_transport_MPI_Cart_create_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Cart_create_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Cart_create");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm_old.data(), &native_ndims, const_cast<void *>(static_cast<const void *>(dims)), const_cast<void *>(static_cast<const void *>(periods)), &native_reorder, native_pointer_comm_cart, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (comm_cart && valid_output) {
    std::memset(comm_cart->bytes, 0, sizeof(comm_cart->bytes));
    std::memcpy(comm_cart->bytes, native_comm_cart.data(), native_comm_cart.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Dims_create(int32_t nnodes, int32_t ndims, int32_t dims[]) {
  Registry &registry = currentRegistry();
  int native_nnodes = static_cast<int>(nnodes);
  int native_ndims = static_cast<int>(ndims);
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Dims_create", 0xfb2dd1cdu, 1u, 3u);
  plumed_mpiw_transport_MPI_Dims_create_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Dims_create_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Dims_create");
  }
  std::int32_t native_return_value = 0;
  dispatch(&native_nnodes, &native_ndims, dims, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Graph_create(plumed_mpiw_handle_storage comm_old, int32_t nnodes, const int32_t indx[], const int32_t edges[], int32_t reorder, plumed_mpiw_handle_storage * comm_graph) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm_old(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm_old.data(), comm_old.bytes, native_comm_old.size());
  validateInputComm(registry, "MPI_Graph_create", native_comm_old.data());
  int native_nnodes = static_cast<int>(nnodes);
  int native_reorder = static_cast<int>(reorder);
  NativeBuffer native_comm_graph(registry.nativeHandleSize("Comm"));
  void *native_pointer_comm_graph = comm_graph ? native_comm_graph.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Graph_create", 0x00eb0f02u, 1u, 6u);
  plumed_mpiw_transport_MPI_Graph_create_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Graph_create_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Graph_create");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm_old.data(), &native_nnodes, const_cast<void *>(static_cast<const void *>(indx)), const_cast<void *>(static_cast<const void *>(edges)), &native_reorder, native_pointer_comm_graph, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (comm_graph && valid_output) {
    std::memset(comm_graph->bytes, 0, sizeof(comm_graph->bytes));
    std::memcpy(comm_graph->bytes, native_comm_graph.data(), native_comm_graph.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Topo_test(plumed_mpiw_handle_storage comm, int32_t * status) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Topo_test", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Topo_test", 0x09cbe810u, 1u, 2u);
  plumed_mpiw_transport_MPI_Topo_test_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Topo_test_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Topo_test");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm.data(), status, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Graphdims_get(plumed_mpiw_handle_storage comm, int32_t * nnodes, int32_t * nedges) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Graphdims_get", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Graphdims_get", 0xc56e163bu, 1u, 3u);
  plumed_mpiw_transport_MPI_Graphdims_get_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Graphdims_get_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Graphdims_get");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm.data(), nnodes, nedges, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Graph_get(plumed_mpiw_handle_storage comm, int32_t maxindex, int32_t maxedges, int32_t indx[], int32_t edges[]) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Graph_get", native_comm.data());
  int native_maxindex = static_cast<int>(maxindex);
  int native_maxedges = static_cast<int>(maxedges);
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Graph_get", 0x12829ea5u, 1u, 5u);
  plumed_mpiw_transport_MPI_Graph_get_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Graph_get_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Graph_get");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm.data(), &native_maxindex, &native_maxedges, indx, edges, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Cartdim_get(plumed_mpiw_handle_storage comm, int32_t * ndims) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Cartdim_get", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Cartdim_get", 0xa979e88cu, 1u, 2u);
  plumed_mpiw_transport_MPI_Cartdim_get_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Cartdim_get_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Cartdim_get");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm.data(), ndims, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Cart_get(plumed_mpiw_handle_storage comm, int32_t maxdims, int32_t dims[], int32_t periods[], int32_t coords[]) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Cart_get", native_comm.data());
  int native_maxdims = static_cast<int>(maxdims);
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Cart_get", 0x0a3d9750u, 1u, 5u);
  plumed_mpiw_transport_MPI_Cart_get_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Cart_get_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Cart_get");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm.data(), &native_maxdims, dims, periods, coords, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Cart_rank(plumed_mpiw_handle_storage comm, const int32_t coords[], int32_t * rank) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Cart_rank", native_comm.data());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Cart_rank", 0xf899fe26u, 1u, 3u);
  plumed_mpiw_transport_MPI_Cart_rank_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Cart_rank_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Cart_rank");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm.data(), const_cast<void *>(static_cast<const void *>(coords)), rank, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Cart_coords(plumed_mpiw_handle_storage comm, int32_t rank, int32_t maxdims, int32_t coords[]) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Cart_coords", native_comm.data());
  int native_rank = static_cast<int>(rank);
  int native_maxdims = static_cast<int>(maxdims);
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Cart_coords", 0xf0bdc30cu, 1u, 4u);
  plumed_mpiw_transport_MPI_Cart_coords_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Cart_coords_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Cart_coords");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm.data(), &native_rank, &native_maxdims, coords, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Graph_neighbors_count(plumed_mpiw_handle_storage comm, int32_t rank, int32_t * nneighbors) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Graph_neighbors_count", native_comm.data());
  int native_rank = static_cast<int>(rank);
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Graph_neighbors_count", 0x5a79d6b6u, 1u, 3u);
  plumed_mpiw_transport_MPI_Graph_neighbors_count_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Graph_neighbors_count_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Graph_neighbors_count");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm.data(), &native_rank, nneighbors, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Graph_neighbors(plumed_mpiw_handle_storage comm, int32_t rank, int32_t maxneighbors, int32_t neighbors[]) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Graph_neighbors", native_comm.data());
  int native_rank = static_cast<int>(rank);
  int native_maxneighbors = static_cast<int>(maxneighbors);
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Graph_neighbors", 0x567d4ff6u, 1u, 4u);
  plumed_mpiw_transport_MPI_Graph_neighbors_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Graph_neighbors_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Graph_neighbors");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm.data(), &native_rank, &native_maxneighbors, neighbors, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Cart_shift(plumed_mpiw_handle_storage comm, int32_t direction, int32_t disp, int32_t * rank_source, int32_t * rank_dest) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Cart_shift", native_comm.data());
  int native_direction = static_cast<int>(direction);
  int native_disp = static_cast<int>(disp);
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Cart_shift", 0x4d227e7cu, 1u, 5u);
  plumed_mpiw_transport_MPI_Cart_shift_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Cart_shift_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Cart_shift");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm.data(), &native_direction, &native_disp, rank_source, rank_dest, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Cart_sub(plumed_mpiw_handle_storage comm, const int32_t remain_dims[], plumed_mpiw_handle_storage * newcomm) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Cart_sub", native_comm.data());
  NativeBuffer native_newcomm(registry.nativeHandleSize("Comm"));
  void *native_pointer_newcomm = newcomm ? native_newcomm.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Cart_sub", 0x1e9d3861u, 1u, 3u);
  plumed_mpiw_transport_MPI_Cart_sub_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Cart_sub_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Cart_sub");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm.data(), const_cast<void *>(static_cast<const void *>(remain_dims)), native_pointer_newcomm, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (newcomm && valid_output) {
    std::memset(newcomm->bytes, 0, sizeof(newcomm->bytes));
    std::memcpy(newcomm->bytes, native_newcomm.data(), native_newcomm.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Cart_map(plumed_mpiw_handle_storage comm, int32_t ndims, const int32_t dims[], const int32_t periods[], int32_t * newrank) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Cart_map", native_comm.data());
  int native_ndims = static_cast<int>(ndims);
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Cart_map", 0xba5a7b4cu, 1u, 5u);
  plumed_mpiw_transport_MPI_Cart_map_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Cart_map_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Cart_map");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm.data(), &native_ndims, const_cast<void *>(static_cast<const void *>(dims)), const_cast<void *>(static_cast<const void *>(periods)), newrank, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Graph_map(plumed_mpiw_handle_storage comm, int32_t nnodes, const int32_t indx[], const int32_t edges[], int32_t * newrank) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Graph_map", native_comm.data());
  int native_nnodes = static_cast<int>(nnodes);
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Graph_map", 0xefd795ebu, 1u, 5u);
  plumed_mpiw_transport_MPI_Graph_map_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Graph_map_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Graph_map");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm.data(), &native_nnodes, const_cast<void *>(static_cast<const void *>(indx)), const_cast<void *>(static_cast<const void *>(edges)), newrank, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Get_processor_name(char * name, int32_t * resultlen) {
  Registry &registry = currentRegistry();
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Get_processor_name", 0xf8e4ccc3u, 1u, 2u);
  plumed_mpiw_transport_MPI_Get_processor_name_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Get_processor_name_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Get_processor_name");
  }
  std::int32_t native_return_value = 0;
  dispatch(name, resultlen, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Errhandler_free(plumed_mpiw_handle_storage * errhandler) {
  Registry &registry = currentRegistry();
  NativeBuffer native_errhandler(registry.nativeHandleSize("Errhandler"));
  if (errhandler) {
    std::memcpy(native_errhandler.data(), errhandler->bytes, native_errhandler.size());
  }
  void *native_pointer_errhandler = errhandler ? native_errhandler.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Errhandler_free", 0xf03f8309u, 1u, 1u);
  plumed_mpiw_transport_MPI_Errhandler_free_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Errhandler_free_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Errhandler_free");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_pointer_errhandler, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (errhandler && valid_output) {
    std::memset(errhandler->bytes, 0, sizeof(errhandler->bytes));
    std::memcpy(errhandler->bytes, native_errhandler.data(), native_errhandler.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Error_class(int32_t errorcode, int32_t * errorclass) {
  Registry &registry = currentRegistry();
  int native_errorcode = static_cast<int>(errorcode);
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Error_class", 0x3b9188a6u, 1u, 2u);
  plumed_mpiw_transport_MPI_Error_class_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Error_class_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Error_class");
  }
  std::int32_t native_return_value = 0;
  dispatch(&native_errorcode, errorclass, &native_return_value);
  return native_return_value;
}
static double plumed_mpiw_bridge_MPI_Wtime(void) {
  Registry &registry = currentRegistry();
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Wtime", 0x75704e2eu, 1u, 0u);
  plumed_mpiw_transport_MPI_Wtime_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Wtime_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Wtime");
  }
  double native_result = 0;
  dispatch(&native_result);
  return native_result;
}
static double plumed_mpiw_bridge_MPI_Wtick(void) {
  Registry &registry = currentRegistry();
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Wtick", 0x84f54b80u, 1u, 0u);
  plumed_mpiw_transport_MPI_Wtick_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Wtick_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Wtick");
  }
  double native_result = 0;
  dispatch(&native_result);
  return native_result;
}
static int32_t plumed_mpiw_bridge_MPI_Finalized(int32_t * flag) {
  Registry &registry = currentRegistry();
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Finalized", 0x105b65c7u, 1u, 1u);
  plumed_mpiw_transport_MPI_Finalized_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Finalized_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Finalized");
  }
  std::int32_t native_return_value = 0;
  dispatch(flag, &native_return_value);
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Comm_get_errhandler(plumed_mpiw_handle_storage comm, plumed_mpiw_handle_storage * errhandler) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Comm_get_errhandler", native_comm.data());
  NativeBuffer native_errhandler(registry.nativeHandleSize("Errhandler"));
  void *native_pointer_errhandler = errhandler ? native_errhandler.data() : nullptr;
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Comm_get_errhandler", 0x74d05159u, 1u, 2u);
  plumed_mpiw_transport_MPI_Comm_get_errhandler_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Comm_get_errhandler_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Comm_get_errhandler");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm.data(), native_pointer_errhandler, &native_return_value);
  const bool valid_output = (native_return_value == 0);
  if (errhandler && valid_output) {
    std::memset(errhandler->bytes, 0, sizeof(errhandler->bytes));
    std::memcpy(errhandler->bytes, native_errhandler.data(), native_errhandler.size());
  }
  return native_return_value;
}
static int32_t plumed_mpiw_bridge_MPI_Comm_set_errhandler(plumed_mpiw_handle_storage comm, plumed_mpiw_handle_storage errhandler) {
  Registry &registry = currentRegistry();
  NativeBuffer native_comm(registry.nativeHandleSize("Comm"));
  std::memcpy(native_comm.data(), comm.bytes, native_comm.size());
  validateInputComm(registry, "MPI_Comm_set_errhandler", native_comm.data());
  NativeBuffer native_errhandler(registry.nativeHandleSize("Errhandler"));
  std::memcpy(native_errhandler.data(), errhandler.bytes, native_errhandler.size());
  plumed_mpiw_generic_fn erased = registry.dispatch("MPI_Comm_set_errhandler", 0xe6af3a8fu, 1u, 2u);
  plumed_mpiw_transport_MPI_Comm_set_errhandler_fn dispatch = reinterpret_cast<plumed_mpiw_transport_MPI_Comm_set_errhandler_fn>(erased);
  if (!dispatch) {
    mpiwrap_unavailable("MPI_Comm_set_errhandler");
  }
  std::int32_t native_return_value = 0;
  dispatch(native_comm.data(), native_errhandler.data(), &native_return_value);
  return native_return_value;
}
