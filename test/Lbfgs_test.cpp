#include "Lbfgs.h"

#include "gtest/gtest.h"
#include "gtest-mpi-listener.hpp"

#include "utils/mpi.h"

using namespace minim;


TEST(LbfgsTest, TestMaxIter) {
  minim::Lbfgs lbfgs = minim::Lbfgs().setMaxIter(10);
  EXPECT_EQ(lbfgs.maxIter, 10);
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  mpiInit(&argc, &argv);

  // Add an MPI listener (https://github.com/LLNL/gtest-mpi-listener)
  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
  ::testing::TestEventListener *l = listeners.Release(listeners.default_result_printer());
  listeners.Append(new GTestMPIListener::MPIWrapperPrinter(l, MPI_COMM_WORLD));

  return RUN_ALL_TESTS();
}
