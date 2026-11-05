/* Build the native adapter from the public, installed Plumed.h only. */
#if defined(PLUMED_MPIW_NATIVE_ADAPTER)
#define __PLUMED_WRAPPER_IMPLEMENTATION 0
#define __PLUMED_WRAPPER_CXX 0
#include "../wrapper/Plumed.h"
#endif
