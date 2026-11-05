// Generated kernel-side conversion of native provider constants.
static int32_t plumed_mpiw_bridge_constant_MPI_FLOAT(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Datatype)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_FLOAT"));
  if (!registry.copyNativeValue("MPI_FLOAT", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Datatype value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_DOUBLE(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Datatype)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_DOUBLE"));
  if (!registry.copyNativeValue("MPI_DOUBLE", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Datatype value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_INT(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Datatype)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_INT"));
  if (!registry.copyNativeValue("MPI_INT", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Datatype value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_CHAR(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Datatype)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_CHAR"));
  if (!registry.copyNativeValue("MPI_CHAR", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Datatype value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_UNSIGNED(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Datatype)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_UNSIGNED"));
  if (!registry.copyNativeValue("MPI_UNSIGNED", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Datatype value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_UNSIGNED_LONG(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Datatype)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_UNSIGNED_LONG"));
  if (!registry.copyNativeValue("MPI_UNSIGNED_LONG", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Datatype value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_UNSIGNED_LONG_LONG(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Datatype)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_UNSIGNED_LONG_LONG"));
  if (!registry.copyNativeValue("MPI_UNSIGNED_LONG_LONG", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Datatype value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_LONG_DOUBLE(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Datatype)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_LONG_DOUBLE"));
  if (!registry.copyNativeValue("MPI_LONG_DOUBLE", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Datatype value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_COMM_WORLD(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Comm)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_COMM_WORLD"));
  if (!registry.copyNativeValue("MPI_COMM_WORLD", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Comm value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_COMM_SELF(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Comm)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_COMM_SELF"));
  if (!registry.copyNativeValue("MPI_COMM_SELF", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Comm value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_COMM_NULL(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Comm)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_COMM_NULL"));
  if (!registry.copyNativeValue("MPI_COMM_NULL", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Comm value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_IN_PLACE(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(void *)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_IN_PLACE"));
  if (!registry.copyNativeValue("MPI_IN_PLACE", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(void *)) {
    return -1;
  }
  void *value = nullptr;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_SUM(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Op)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_SUM"));
  if (!registry.copyNativeValue("MPI_SUM", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Op value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_PROD(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Op)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_PROD"));
  if (!registry.copyNativeValue("MPI_PROD", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Op value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_MAX(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Op)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_MAX"));
  if (!registry.copyNativeValue("MPI_MAX", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Op value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_MIN(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Op)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_MIN"));
  if (!registry.copyNativeValue("MPI_MIN", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Op value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_STATUS_IGNORE(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Status *)) {
    return -1;
  }
  MPI_Status * value = nullptr;
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_STATUSES_IGNORE(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Status *)) {
    return -1;
  }
  MPI_Status * value = nullptr;
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_GROUP_NULL(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Group)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_GROUP_NULL"));
  if (!registry.copyNativeValue("MPI_GROUP_NULL", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Group value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_GROUP_EMPTY(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Group)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_GROUP_EMPTY"));
  if (!registry.copyNativeValue("MPI_GROUP_EMPTY", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Group value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_REQUEST_NULL(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Request)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_REQUEST_NULL"));
  if (!registry.copyNativeValue("MPI_REQUEST_NULL", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Request value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_DATATYPE_NULL(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Datatype)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_DATATYPE_NULL"));
  if (!registry.copyNativeValue("MPI_DATATYPE_NULL", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Datatype value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_OP_NULL(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Op)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_OP_NULL"));
  if (!registry.copyNativeValue("MPI_OP_NULL", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Op value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_ERRHANDLER_NULL(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Errhandler)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_ERRHANDLER_NULL"));
  if (!registry.copyNativeValue("MPI_ERRHANDLER_NULL", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Errhandler value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_ERRORS_ARE_FATAL(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Errhandler)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_ERRORS_ARE_FATAL"));
  if (!registry.copyNativeValue("MPI_ERRORS_ARE_FATAL", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Errhandler value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_ERRORS_RETURN(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Errhandler)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_ERRORS_RETURN"));
  if (!registry.copyNativeValue("MPI_ERRORS_RETURN", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Errhandler value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_BYTE(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Datatype)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_BYTE"));
  if (!registry.copyNativeValue("MPI_BYTE", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Datatype value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_PACKED(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Datatype)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_PACKED"));
  if (!registry.copyNativeValue("MPI_PACKED", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Datatype value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_SHORT(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Datatype)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_SHORT"));
  if (!registry.copyNativeValue("MPI_SHORT", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Datatype value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_LONG(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Datatype)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_LONG"));
  if (!registry.copyNativeValue("MPI_LONG", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Datatype value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_UNSIGNED_CHAR(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Datatype)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_UNSIGNED_CHAR"));
  if (!registry.copyNativeValue("MPI_UNSIGNED_CHAR", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Datatype value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_UNSIGNED_SHORT(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Datatype)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_UNSIGNED_SHORT"));
  if (!registry.copyNativeValue("MPI_UNSIGNED_SHORT", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Datatype value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_LONG_LONG_INT(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Datatype)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_LONG_LONG_INT"));
  if (!registry.copyNativeValue("MPI_LONG_LONG_INT", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Datatype value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_LONG_LONG(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Datatype)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_LONG_LONG"));
  if (!registry.copyNativeValue("MPI_LONG_LONG", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Datatype value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_FLOAT_INT(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Datatype)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_FLOAT_INT"));
  if (!registry.copyNativeValue("MPI_FLOAT_INT", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Datatype value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_DOUBLE_INT(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Datatype)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_DOUBLE_INT"));
  if (!registry.copyNativeValue("MPI_DOUBLE_INT", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Datatype value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_LONG_INT(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Datatype)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_LONG_INT"));
  if (!registry.copyNativeValue("MPI_LONG_INT", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Datatype value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_2INT(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Datatype)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_2INT"));
  if (!registry.copyNativeValue("MPI_2INT", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Datatype value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_SHORT_INT(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Datatype)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_SHORT_INT"));
  if (!registry.copyNativeValue("MPI_SHORT_INT", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Datatype value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_LONG_DOUBLE_INT(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Datatype)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_LONG_DOUBLE_INT"));
  if (!registry.copyNativeValue("MPI_LONG_DOUBLE_INT", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Datatype value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_MAXLOC(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Op)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_MAXLOC"));
  if (!registry.copyNativeValue("MPI_MAXLOC", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Op value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_MINLOC(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Op)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_MINLOC"));
  if (!registry.copyNativeValue("MPI_MINLOC", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Op value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_LAND(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Op)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_LAND"));
  if (!registry.copyNativeValue("MPI_LAND", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Op value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_BAND(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Op)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_BAND"));
  if (!registry.copyNativeValue("MPI_BAND", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Op value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_LOR(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Op)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_LOR"));
  if (!registry.copyNativeValue("MPI_LOR", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Op value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_BOR(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Op)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_BOR"));
  if (!registry.copyNativeValue("MPI_BOR", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Op value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_LXOR(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Op)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_LXOR"));
  if (!registry.copyNativeValue("MPI_LXOR", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Op value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_BXOR(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Op)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_BXOR"));
  if (!registry.copyNativeValue("MPI_BXOR", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Op value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_UNDEFINED(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_UNDEFINED"));
  if (!registry.copyNativeValue("MPI_UNDEFINED", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_ANY_SOURCE(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_ANY_SOURCE"));
  if (!registry.copyNativeValue("MPI_ANY_SOURCE", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_ANY_TAG(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_ANY_TAG"));
  if (!registry.copyNativeValue("MPI_ANY_TAG", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_PROC_NULL(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_PROC_NULL"));
  if (!registry.copyNativeValue("MPI_PROC_NULL", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_IDENT(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_IDENT"));
  if (!registry.copyNativeValue("MPI_IDENT", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_CONGRUENT(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_CONGRUENT"));
  if (!registry.copyNativeValue("MPI_CONGRUENT", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_SIMILAR(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_SIMILAR"));
  if (!registry.copyNativeValue("MPI_SIMILAR", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_UNEQUAL(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_UNEQUAL"));
  if (!registry.copyNativeValue("MPI_UNEQUAL", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_CART(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_CART"));
  if (!registry.copyNativeValue("MPI_CART", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_GRAPH(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_GRAPH"));
  if (!registry.copyNativeValue("MPI_GRAPH", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_MAX_PROCESSOR_NAME(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_MAX_PROCESSOR_NAME"));
  if (!registry.copyNativeValue("MPI_MAX_PROCESSOR_NAME", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_MAX_ERROR_STRING(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_MAX_ERROR_STRING"));
  if (!registry.copyNativeValue("MPI_MAX_ERROR_STRING", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_BSEND_OVERHEAD(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_BSEND_OVERHEAD"));
  if (!registry.copyNativeValue("MPI_BSEND_OVERHEAD", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_SUCCESS(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_SUCCESS"));
  if (!registry.copyNativeValue("MPI_SUCCESS", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_ERR_BUFFER(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_ERR_BUFFER"));
  if (!registry.copyNativeValue("MPI_ERR_BUFFER", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_ERR_COUNT(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_ERR_COUNT"));
  if (!registry.copyNativeValue("MPI_ERR_COUNT", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_ERR_TYPE(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_ERR_TYPE"));
  if (!registry.copyNativeValue("MPI_ERR_TYPE", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_ERR_TAG(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_ERR_TAG"));
  if (!registry.copyNativeValue("MPI_ERR_TAG", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_ERR_COMM(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_ERR_COMM"));
  if (!registry.copyNativeValue("MPI_ERR_COMM", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_ERR_RANK(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_ERR_RANK"));
  if (!registry.copyNativeValue("MPI_ERR_RANK", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_ERR_REQUEST(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_ERR_REQUEST"));
  if (!registry.copyNativeValue("MPI_ERR_REQUEST", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_ERR_ROOT(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_ERR_ROOT"));
  if (!registry.copyNativeValue("MPI_ERR_ROOT", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_ERR_GROUP(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_ERR_GROUP"));
  if (!registry.copyNativeValue("MPI_ERR_GROUP", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_ERR_OP(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_ERR_OP"));
  if (!registry.copyNativeValue("MPI_ERR_OP", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_ERR_TOPOLOGY(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_ERR_TOPOLOGY"));
  if (!registry.copyNativeValue("MPI_ERR_TOPOLOGY", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_ERR_DIMS(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_ERR_DIMS"));
  if (!registry.copyNativeValue("MPI_ERR_DIMS", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_ERR_ARG(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_ERR_ARG"));
  if (!registry.copyNativeValue("MPI_ERR_ARG", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_ERR_UNKNOWN(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_ERR_UNKNOWN"));
  if (!registry.copyNativeValue("MPI_ERR_UNKNOWN", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_ERR_TRUNCATE(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_ERR_TRUNCATE"));
  if (!registry.copyNativeValue("MPI_ERR_TRUNCATE", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_ERR_OTHER(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_ERR_OTHER"));
  if (!registry.copyNativeValue("MPI_ERR_OTHER", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_ERR_INTERN(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_ERR_INTERN"));
  if (!registry.copyNativeValue("MPI_ERR_INTERN", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_ERR_PENDING(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_ERR_PENDING"));
  if (!registry.copyNativeValue("MPI_ERR_PENDING", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_BOTTOM(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(void *)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_BOTTOM"));
  if (!registry.copyNativeValue("MPI_BOTTOM", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(void *)) {
    return -1;
  }
  void *value = nullptr;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_ERR_IN_STATUS(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(std::int32_t)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_ERR_IN_STATUS"));
  if (!registry.copyNativeValue("MPI_ERR_IN_STATUS", native_value.data(), native_value.size())) {
    return -1;
  }
  if (native_value.size() != sizeof(std::int32_t)) {
    return -1;
  }
  std::int32_t value = 0;
  std::memcpy(&value, native_value.data(), sizeof(value));
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_SIGNED_CHAR(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Datatype)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_SIGNED_CHAR"));
  if (!registry.copyNativeValue("MPI_SIGNED_CHAR", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Datatype value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_C_BOOL(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Datatype)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_C_BOOL"));
  if (!registry.copyNativeValue("MPI_C_BOOL", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Datatype value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
static int32_t plumed_mpiw_bridge_constant_MPI_INT64_T(void *output, uint32_t capacity) {
  if (!output || capacity < sizeof(MPI_Datatype)) {
    return -1;
  }
  Registry &registry = currentRegistry();
  NativeBuffer native_value(registry.nativeConstantSize("MPI_INT64_T"));
  if (!registry.copyNativeValue("MPI_INT64_T", native_value.data(), native_value.size())) {
    return -1;
  }
  MPI_Datatype value{};
  if (native_value.size() > sizeof(value.storage.bytes)) {
    return -1;
  }
  std::memcpy(value.storage.bytes, native_value.data(), native_value.size());
  std::memcpy(output, &value, sizeof(value));
  return 0;
}
