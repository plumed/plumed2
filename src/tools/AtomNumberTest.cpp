#include "gtest/gtest.h"
#include "gmock/gmock.h" // for the matcher
#include "AtomNumber.h"

namespace PLMD {
TEST(AtomNumber,Constructor) {

  AtomNumber *atomNumber = new AtomNumber();

  EXPECT_EQ(atomNumber->index(), 0);
}
}