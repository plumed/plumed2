#if defined(__PLUMED_MPIWRAP) && !defined(__PLUMED_HAS_MPI)
/* Registry and frontend support are present in every wrapper-mode kernel so
   runtime installation remains available even for an explicitly serial
   build. The embedded provider and launcher are optional. */
#include "provider-registry.cpp.inc"
#include "provider-pointer-definitions.inc"
#include "provider-frontend.cpp.inc"
#include "provider-fallbacks.cpp.inc"
#if defined(__PLUMED_HAS_THREAD_MPI)
#include "thread-mpi-adapter.cpp.inc"
#else
#include "provider-loader.cpp.inc"
namespace PLMD {
namespace mpiwrap_provider {
void ensureProviderSelected() {
  if (forceEnvironmentProvider()) {
    return;
  }
  Registry &registry = currentRegistry();
  if (registry.installed()) {
    return;
  }
  const char *selected = std::getenv("PLUMED_MPI_WRAPPER");
  if (selected && selected[0]) {
    installEnvironmentProvider(registry, selected);
  }
}

}  // namespace mpiwrap_provider
}  // namespace PLMD
#endif
#endif
