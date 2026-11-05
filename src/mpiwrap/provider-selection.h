/* Kernel-side environment selection; this is not part of the host MPI ABI. */
#ifndef PLUMED_MPIWRAP_PROVIDER_SELECTION_H
#define PLUMED_MPIWRAP_PROVIDER_SELECTION_H

namespace PLMD {
namespace mpiwrap_provider {
// Return true when an explicit troubleshooting override has installed its
// provider. Throws on invalid configuration or an attempted late replacement.
bool forceEnvironmentProvider();
}
}

#endif
