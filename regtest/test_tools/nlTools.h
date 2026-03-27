#ifndef __PLUMED_TEST_NL
#define __PLUMED_TEST_NL
#include "plumed/tools/NeighborList.h"

#include <iosfwd>

namespace PLMD {
namespace test {
void printNeighbors( std::string prefix,
                     const PLMD::NeighborList& nl,
                     unsigned const nat,
                     std::ostream& ofs);
} //namespace test
} //namespace PLMD
#endif // __PLUMED_TEST_NL
