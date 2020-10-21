#include "gtest/gtest.h"
#include "gmock/gmock.h" // for the matcher
#include "bar.h"

namespace PLMD {
TEST(Foo,Is42) {
  EXPECT_EQ(bar(), 42);
}

TEST(Foo,IsStill42) {
  // Another way to spell it is to use EXPECT_THAT with "matchers", which
  // are named predicates. You'll want to familiarize yourself with these
  // because they are composable to create very sophisticated matchers.
  // e.g. Contains(Key(Le(5))) will verify that a map/unordered_map
  // has at least one key <= 5
  EXPECT_THAT(bar(),testing::Eq(42));
}
}