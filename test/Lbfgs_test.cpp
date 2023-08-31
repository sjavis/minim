#include "test_main.cpp"
#include "Lbfgs.h"

using namespace minim;


TEST(LbfgsTest, TestMaxIter) {
  minim::Lbfgs lbfgs = minim::Lbfgs().setMaxIter(10);
  EXPECT_EQ(lbfgs.maxIter, 10);
}
