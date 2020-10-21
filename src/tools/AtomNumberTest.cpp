#include "gtest/gtest.h"
#include "gmock/gmock.h" // for the matcher
#include "AtomNumber.h"

namespace PLMD {
TEST(AtomNumber,Constructor) {
  AtomNumber *atomNumber = new AtomNumber();
}
TEST(AtomNumber,Index) {
  AtomNumber *atomNumber = new AtomNumber();
  EXPECT_EQ(atomNumber->index(), 0);
  atomNumber->setIndex(42);
  EXPECT_EQ(atomNumber->index(), 42);
}
TEST(AtomNumber,Serial) {
  AtomNumber *atomNumber = new AtomNumber();
  EXPECT_EQ(atomNumber->serial(), 1);
  atomNumber->setSerial(42);
  EXPECT_EQ(atomNumber->serial(), 42);
}
}
