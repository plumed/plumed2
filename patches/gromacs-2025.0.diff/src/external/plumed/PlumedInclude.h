#define __PLUMED_WRAPPER_FORTRAN 0 // NOLINT(bugprone-reserved-identifier)

#define __PLUMED_WRAPPER_LINK_RUNTIME 1 // NOLINT(bugprone-reserved-identifier)
#define __PLUMED_WRAPPER_EXTERN 1       // NOLINT(bugprone-reserved-identifier)

#define __PLUMED_WRAPPER_CXX 1            // NOLINT(bugprone-reserved-identifier)
#define __PLUMED_WRAPPER_LIBCXX11 1       // NOLINT(bugprone-reserved-identifier)
#define __PLUMED_WRAPPER_LIBCXX17 1       // NOLINT(bugprone-reserved-identifier)
#ifndef __PLUMED_WRAPPER_IMPLEMENTATION  // NOLINT(bugprone-reserved-identifier)
#  define __PLUMED_WRAPPER_IMPLEMENTATION 0 // NOLINT(bugprone-reserved-identifier)
#endif
#define __PLUMED_HAS_DLOPEN               // NOLINT(bugprone-reserved-identifier)

#include "external/plumed/Plumed.h"