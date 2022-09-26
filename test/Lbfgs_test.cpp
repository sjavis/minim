#include "Lbfgs.h"

#include "gtest/gtest.h"

TEST(LbfgsTest, TestMaxIter) {
  minim::Lbfgs lbfgs = minim::Lbfgs().setMaxIter(10);
  EXPECT_EQ(lbfgs.maxIter, 10);
}
