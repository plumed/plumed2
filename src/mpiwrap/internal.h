// PLUMED_HEADER_ONLY
// Runtime validation and typed lookup for one installed MPI provider.
#ifndef __PLUMED_mpiwrap_internal_h
#define __PLUMED_mpiwrap_internal_h

#include "mpiwrap.h"
#if !defined(__PLUMED_HAS_MPI)

#include <cstddef>
#include <cstdint>
#include <atomic>
#include <mutex>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace PLMD {
namespace mpiwrap {

/* Internal API values. The v3 provider descriptor carries only the boolean
   ranks_are_threads, not these kernel-side policy labels. */
#define PLUMED_MPIW_RUNTIME_PROCESS_MPI 1
#define PLUMED_MPIW_RUNTIME_THREAD_MPI 2

class Registry {
public:
  struct FunctionEntry {
    std::uint32_t signature_id;
    std::uint32_t signature_version;
    std::uint32_t argument_count;
    plumed_mpi_generic_pointer function;
  };
  struct ValueEntry {
    std::uint32_t value_type;
    std::uint32_t native_size;
    plumed_mpi_constant_accessor_pointer accessor;
    bool operator==(const ValueEntry &other) const {
      return value_type == other.value_type && native_size == other.native_size &&
             accessor == other.accessor;
    }
  };

  Registry();
  ~Registry();
  Registry(const Registry &) = delete;
  Registry &operator=(const Registry &) = delete;

  bool install(const plumed_mpi_provider *descriptor,
               std::string &error,void *provider_library=nullptr);
  bool has(const char *name) const;
  plumed_mpi_generic_pointer dispatch(const char *name,
                                      std::uint32_t signature_id,
                                      std::uint32_t signature_version,
                                      std::uint32_t argument_count) const;
  bool installed() const noexcept {
    return ready_.load(std::memory_order_acquire);
  }
  /* Immutable, owned copy for injection into a nested compatible kernel. */
  const plumed_mpi_provider *descriptorSnapshot() const noexcept;
  void markProviderValueUsed(const char *name);
  bool beginInternalLaunch(const plumed_mpi_provider *descriptor,
                           std::string &error);
  void finishInternalLaunch() noexcept;
  std::uint32_t runtimeKind() const;
  const ProviderBindings &bindings() const {
    if (!ready_.load(std::memory_order_acquire)) {
      throw std::runtime_error("MPI provider is not installed");
    }
    return bindings_;
  }
  std::uint32_t nativeCommSize() const;
  std::uint32_t nativeHandleSize(const char *kind) const;
  std::uint32_t nativeStatusSize() const;
  std::uint32_t nativeStatusOffset(const char *field) const;
  std::uint32_t nativeFintSize() const;
  std::uint32_t nativeConstantSize(const char *name) const;
  bool copyNativeValue(const char *name, void *output, std::size_t size) const;
  bool ranksAreThreads() const;
  bool hasCapability(const char *name) const;
  bool copyValue(const char *name, void *output, std::size_t size) const;

  std::string providerName() const;

private:
  struct DescriptorSnapshot;
  mutable std::mutex mutex_;
  std::unordered_map<std::string, FunctionEntry> functions_;
  std::unordered_map<std::string, std::string> metadata_;
  std::unordered_map<std::string, ValueEntry> values_;
  ProviderBindings bindings_;
  std::string provider_name_;
  std::vector<unsigned char> runtime_identity_;
  plumed_mpi_generic_pointer descriptor_function_comm_rank_;
  plumed_mpi_generic_pointer descriptor_function_initialized_;
  std::uint32_t runtime_kind_;
  std::uint32_t native_comm_size_;
  std::uint32_t native_handle_sizes_[6];
  std::uint32_t native_status_size_;
  std::uint32_t native_status_offsets_[3];
  std::uint32_t native_fint_size_;
  std::vector<void *> provider_libraries_;
  std::unique_ptr<DescriptorSnapshot> descriptor_snapshot_;
  bool installed_;
  bool provider_value_used_;
  bool internal_launch_started_;
  bool internal_launch_active_;
  std::atomic<bool> ready_;
};

// The loaded kernel owns one immutable provider registry shared by all its
// rank threads. Provider-dependent use fixes its selection for the kernel.
Registry &currentRegistry();
// Return true when an explicit troubleshooting override has installed its
// provider. Throws on invalid configuration or an attempted late replacement.
bool forceEnvironmentProvider();
// Native communicator import is the only implicit compiler-discovery boundary.
void ensureProviderForCommunicator();
// Honors an explicit environment provider on provider-dependent use.
// Embedded thread-MPI is installed only by its explicit launcher.
void ensureProviderSelected();
MPI_Comm importComm(Registry &registry, const void *native);
bool exportComm(Registry &registry, MPI_Comm value,
                plumed_mpi_handle_storage *storage);
bool exportNativeComm(Registry &registry, MPI_Comm value,
                      void *output,std::size_t capacity);
void validateInputComm(Registry &registry, const char *function,
                       const void *native_comm);

}  // namespace mpiwrap
}  // namespace PLMD

#endif // provider implementation
#endif
