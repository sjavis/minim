#include "Lbfgs.h"

#include "gtest/gtest.h"

#include "utils/mpi.h"

using namespace minim;


TEST(LbfgsTest, TestMaxIter) {
  minim::Lbfgs lbfgs = minim::Lbfgs().setMaxIter(10);
  EXPECT_EQ(lbfgs.maxIter, 10);
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  mpiInit(&argc, &argv);

  // Ensure only one processor prints
  ::testing::TestEventListeners& listeners =
      ::testing::UnitTest::GetInstance()->listeners();
  if (mpi.rank != 0) {
      delete listeners.Release(listeners.default_result_printer());
  }

  return RUN_ALL_TESTS();
}
