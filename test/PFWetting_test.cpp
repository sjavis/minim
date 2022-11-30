#include "PFWetting.h"

#include "gtest/gtest.h"
#include "gtest-mpi-listener.hpp"
#include "ArraysMatch.h"

#include "State.h"
#include "utils/vec.h"
#include "utils/mpi.h"

using namespace minim;


TEST(PFWettingTest, TestBulkEnergy) {
  PFWetting pot;

  // Constant bulk fluid
  pot.setGridSize({1,1,1});
  State s1 = pot.newState({1});
  EXPECT_FLOAT_EQ(s1.energy(), 0);
  EXPECT_TRUE(ArraysNear(s1.gradient(), {0}, 1e-6));
  s1.coords({0.1});
  EXPECT_FLOAT_EQ(s1.energy(), 0.245025);
  EXPECT_TRUE(ArraysNear(s1.gradient(), {-0.099}, 1e-6));

  // Bulk fluid gradient
  State s2x = pot.setGridSize({5,1,1}).setSolid({1,0,0,0,1}).newState({0,-1,0,1,0});
  EXPECT_FLOAT_EQ(s2x.energy(), 1.25);
  EXPECT_TRUE(ArraysNear(s2x.gradient(), {0,-1,0,1,0}, 1e-6));
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  MPI_Init(&argc, &argv);
  mpi.getSizeRank(MPI_COMM_WORLD);

  // Add an MPI listener (https://github.com/LLNL/gtest-mpi-listener)
  ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
  ::testing::TestEventListener *l = listeners.Release(listeners.default_result_printer());
  listeners.Append(new GTestMPIListener::MPIWrapperPrinter(l, MPI_COMM_WORLD));

  return RUN_ALL_TESTS();
}
