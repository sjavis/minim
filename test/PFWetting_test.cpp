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
  pot.setInterfaceSize(1/sqrt(2.0));
  pot.setSurfaceTension(sqrt(8.0/9));

  // Constant bulk fluid
  State s1 = pot.setGridSize({1,1,1}).newState({1});
  EXPECT_FLOAT_EQ(s1.energy(), 0);
  EXPECT_TRUE(ArraysNear(s1.gradient(), {0}, 1e-6));
  s1.coords({0.1});
  EXPECT_FLOAT_EQ(s1.energy(), 0.245025);
  EXPECT_TRUE(ArraysNear(s1.gradient(), {-0.099}, 1e-6));

  // Bulk fluid gradient
  State s2x = pot.setGridSize({5,1,1}).setSolid({1,0,0,0,1}).newState({0,-1,0,1,0});
  State s2y = pot.setGridSize({1,5,1}).setSolid({1,0,0,0,1}).newState({0,-1,0,1,0});
  State s2z = pot.setGridSize({1,1,5}).setSolid({1,0,0,0,1}).newState({0,-1,0,1,0});
  EXPECT_FLOAT_EQ(s2x.energy(), 1.25);
  EXPECT_FLOAT_EQ(s2y.energy(), 1.25);
  EXPECT_FLOAT_EQ(s2z.energy(), 1.25);
  EXPECT_TRUE(ArraysNear(s2x.gradient(), {0,-1,0,1,0}, 1e-6));
  EXPECT_TRUE(ArraysNear(s2y.gradient(), {0,-1,0,1,0}, 1e-6));
  EXPECT_TRUE(ArraysNear(s2z.gradient(), {0,-1,0,1,0}, 1e-6));
  // Test periodic boundary
  State s3x = pot.setGridSize({3,1,1}).setSolid({0,1,0}).newState({-1,0,1});
  State s3y = pot.setGridSize({1,3,1}).setSolid({0,1,0}).newState({-1,0,1});
  State s3z = pot.setGridSize({1,1,3}).setSolid({0,1,0}).newState({-1,0,1});
  EXPECT_FLOAT_EQ(s3x.energy(), 2);
  EXPECT_FLOAT_EQ(s3y.energy(), 2);
  EXPECT_FLOAT_EQ(s3z.energy(), 2);
  EXPECT_TRUE(ArraysNear(s3x.gradient(), {-2,0,2}, 1e-6));
  EXPECT_TRUE(ArraysNear(s3y.gradient(), {-2,0,2}, 1e-6));
  EXPECT_TRUE(ArraysNear(s3z.gradient(), {-2,0,2}, 1e-6));
}


TEST(PFWettingTest, TestExternalForce) {
  PFWetting pot;
  pot.setInterfaceSize(1/sqrt(2.0));
  pot.setSurfaceTension(sqrt(8.0/9));

  // Test force in x, y, z directions
  auto stateForceX = pot.setGridSize({2,2,2}).setForce({-2,0,0}).newState({1,1,1,1,1,1,1,1});
  auto stateForceY = pot.setGridSize({2,2,2}).setForce({0,-2,0}).newState({1,1,1,1,1,1,1,1});
  auto stateForceZ = pot.setGridSize({2,2,2}).setForce({0,0,-2}).newState({1,1,1,1,1,1,1,1});
  EXPECT_FLOAT_EQ(stateForceX.energy(), 0);
  EXPECT_FLOAT_EQ(stateForceY.energy(), 0);
  EXPECT_FLOAT_EQ(stateForceZ.energy(), 0);
  EXPECT_TRUE(ArraysNear(stateForceX.gradient(), {-.5,-.5,-.5,-.5,  .5, .5, .5, .5}, 1e-6));
  EXPECT_TRUE(ArraysNear(stateForceY.gradient(), {-.5,-.5, .5, .5, -.5,-.5, .5, .5}, 1e-6));
  EXPECT_TRUE(ArraysNear(stateForceZ.gradient(), {-.5, .5,-.5, .5, -.5, .5,-.5, .5}, 1e-6));

  // Test force with non-constant phi
  auto stateForce1 = pot.setGridSize({2,2,2}).setForce({-4,0,0}).newState({-1,-1,-1,-1, 1,1,1,1});
  for (auto el=stateForce1.pot->elements.begin(); el!=stateForce1.pot->elements.end(); el++) {
    if (el->type==0) stateForce1.pot->elements.erase(el); // Remove the bulk fluid energy elements
  }
  EXPECT_FLOAT_EQ(stateForce1.energy(), 8);
  EXPECT_TRUE(ArraysNear(stateForce1.gradient(), {-1,-1,-1,-1,  1, 1, 1, 1}, 1e-6));
}


TEST(PFWettingTest, TestSurfaceEnergy) {
  PFWetting pot;
  pot.setInterfaceSize(1/sqrt(2.0));
  pot.setSurfaceTension(sqrt(8.0/9));
  pot.setGridSize({2,1,1}).setContactAngle({90,60}).setSolid({1,0});
  auto state = pot.newState({0, 0.5});
  for (auto el=state.pot->elements.begin(); el!=state.pot->elements.end(); el++) {
    if (el->type==0) state.pot->elements.erase(el); // Remove the bulk fluid energy elements
  }
  EXPECT_FLOAT_EQ(state.energy(), 0.5*sqrt(2.0)*(-1.0/12));
  EXPECT_TRUE(ArraysNear(state.gradient(), {0, 0.5*sqrt(2.0)*(-0.25)}, 1e-6));
}


TEST(PFWettingTest, TestPressureConstraint) {
  PFWetting pot;
  pot.setInterfaceSize(1/sqrt(2.0));
  pot.setSurfaceTension(sqrt(8.0/9));
  pot.setGridSize({6,1,1}).setSolid({1,0,0,0,0,1}).setPressure(10);
  auto state = pot.newState({1,1,1,1,1,1});
  EXPECT_FLOAT_EQ(state.energy(), -30);
  EXPECT_TRUE(ArraysNear(state.gradient(), {0, -2.5, -5, -5, -2.5, 0}, 1e-6));
}


TEST(PFWettingTest, TestVolumeConstraint) {
  PFWetting pot;
  pot.setInterfaceSize(1/sqrt(2.0));
  pot.setSurfaceTension(sqrt(8.0/9));
  pot.setGridSize({6,1,1}).setSolid({1,0,0,0,0,1}).setVolume(1, 100);
  auto state = pot.newState({1,1,1,1,1,1});
  EXPECT_FLOAT_EQ(state.energy(), 400);
  EXPECT_TRUE(ArraysNear(state.gradient(), {0, 100, 200, 200, 100, 0}, 1e-6));
}


TEST(PFWettingTest, TestResolution) {
  PFWetting pot;
  pot.setInterfaceSize(1/sqrt(2.0));
  pot.setSurfaceTension(sqrt(8.0/9));
  EXPECT_FLOAT_EQ(pot.resolution, 1);

  pot.setResolution(2);
  ASSERT_FLOAT_EQ(pot.resolution, 2);

  // Bulk fluid
  State s1 = pot.setGridSize({5,1,1}).setSolid({1,0,0,0,1}).newState({0,-1,0,1,0});
  EXPECT_FLOAT_EQ(s1.energy(), 4); // Bulk: 2, Gradient: 2
  EXPECT_TRUE(ArraysNear(s1.gradient(), {0,-2,0,2,0}, 1e-6));

  // External Force
  pot.setGridSize({2,2,2}).setSolid({0,0,0,0,0,0,0,0}).setForce({-4,0,0});
  auto s2 = pot.newState({-1,-1,-1,-1, 1,1,1,1});
  for (auto el=s2.pot->elements.begin(); el!=s2.pot->elements.end(); el++) {
    if (el->type==0) s2.pot->elements.erase(el); // Remove the bulk fluid energy elements
  }
  EXPECT_FLOAT_EQ(s2.energy(), 128);
  EXPECT_TRUE(ArraysNear(s2.gradient(), {-16,-16,-16,-16,  16, 16, 16, 16}, 1e-6));
  pot.setForce({0,0,0});

  // Surface energy
  pot.setGridSize({2,1,1}).setSolid({1,0}).setContactAngle({90,60});
  auto s3 = pot.newState({0, 0.5});
  for (auto el=s3.pot->elements.begin(); el!=s3.pot->elements.end(); el++) {
    if (el->type==0) s3.pot->elements.erase(el); // Remove the bulk fluid energy elements
  }
  EXPECT_FLOAT_EQ(s3.energy(), 0.5*sqrt(2.0)*(-1.0/12)*4);
  EXPECT_TRUE(ArraysNear(s3.gradient(), {0, 0.5*sqrt(2.0)*(-0.25)*4}, 1e-6));

  // Pressure
  pot.setGridSize({6,1,1}).setSolid({1,0,0,0,0,1}).setContactAngle({90,90,90,90,90,90}).setPressure(10);
  auto s4 = pot.newState({1,1,1,1,1,1});
  EXPECT_FLOAT_EQ(s4.energy(), -240);
  EXPECT_TRUE(ArraysNear(s4.gradient(), {0, -20, -40, -40, -20, 0}, 1e-6));
  pot.setPressure(0);

  // Volume
  pot.setGridSize({6,1,1}).setSolid({1,0,0,0,0,1}).setContactAngle({90,90,90,90,90,90}).setVolume(8, 100);
  auto s5 = pot.newState({1,1,1,1,1,1});
  EXPECT_FLOAT_EQ(s5.energy(), 64*400);
  EXPECT_TRUE(ArraysNear(s5.gradient(), {0, 64*100, 64*200, 64*200, 64*100, 0}, 1e-6));
}


TEST(PFWettingTest, TestNFluid) {
  PFWetting pot;
  EXPECT_FLOAT_EQ(pot.nFluid, 2);
  pot.setNFluid(3);
  EXPECT_FLOAT_EQ(pot.nFluid, 3);
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
