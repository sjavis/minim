#include "test_main.cpp"
#include "potentials/PhaseField.h"

#include "State.h"
#include "communicators/CommGrid.h"
#include "utils/vec.h"

using namespace minim;

void tmp(const Communicator& comm) {
    auto commGrid = static_cast<const CommGrid&>(comm);
}

TEST(PhaseFieldTest, gridSizeMpi) {
  PhaseField pot;
  EXPECT_NO_THROW({
    pot.setGridSize({2,1,2});
    State s(pot, {0,0,0,0});
  });
  EXPECT_ANY_THROW({
    pot.setGridSize({3,1,1});
    State s(pot, {0,0,0});
  });
  EXPECT_NO_THROW({
    pot.setGridSize({3,1,1});
    State s(pot, {0,0,0}, {0});
  });
}


TEST(PhaseFieldTest, TestBulkEnergy) {
  PhaseField pot;

  // Constant bulk fluid
  State s1(pot.setGridSize({1,1,1}), {1}, {0});
  ASSERT_FLOAT_EQ(static_cast<PhaseField&>(*s1.pot).kappa[0], 3);
  ASSERT_FLOAT_EQ(static_cast<PhaseField&>(*s1.pot).kappaP[0], 3);
  EXPECT_FLOAT_EQ(s1.allEnergy(), 0);
  EXPECT_TRUE(ArraysNear(s1.allGradient(), {0}));
  s1.coords({0.1});
  EXPECT_FLOAT_EQ(s1.allEnergy(), 3.0/16*1.21*0.81);
  EXPECT_TRUE(ArraysNear(s1.allGradient(), {-3.0/4*0.099}));

  // Bulk fluid gradient
  // Gradient is 1 over the 1 and 2 half fluid nodes
  State s2x = pot.setGridSize({5,1,1}).setSolid({1,0,0,0,1}).newState({0,-1,0,1,0}, {0});
  State s2y = pot.setGridSize({1,5,1}).setSolid({1,0,0,0,1}).newState({0,-1,0,1,0}, {0});
  State s2z = pot.setGridSize({1,1,5}).setSolid({1,0,0,0,1}).newState({0,-1,0,1,0}, {0});
  double ebulk = 3.0 / 16;
  double egrad = 2 * (3.0/4);
  EXPECT_FLOAT_EQ(s2x.allEnergy(), ebulk+egrad);
  EXPECT_FLOAT_EQ(s2y.allEnergy(), ebulk+egrad);
  EXPECT_FLOAT_EQ(s2z.allEnergy(), ebulk+egrad);
  EXPECT_TRUE(ArraysNear(s2x.allGradient(), {0,-1.5,0,1.5,0}));
  EXPECT_TRUE(ArraysNear(s2y.allGradient(), {0,-1.5,0,1.5,0}));
  EXPECT_TRUE(ArraysNear(s2z.allGradient(), {0,-1.5,0,1.5,0}));
  // Test periodic boundary
  // Gradient is 2 over the 2 half fluid nodes
  State s3x = pot.setGridSize({3,1,1}).setSolid({0,1,0}).newState({-1,0,1}, {0});
  State s3y = pot.setGridSize({1,3,1}).setSolid({0,1,0}).newState({-1,0,1}, {0});
  State s3z = pot.setGridSize({1,1,3}).setSolid({0,1,0}).newState({-1,0,1}, {0});
  egrad = 3.0/4 * 4;
  EXPECT_FLOAT_EQ(s3x.allEnergy(), egrad);
  EXPECT_FLOAT_EQ(s3y.allEnergy(), egrad);
  EXPECT_FLOAT_EQ(s3z.allEnergy(), egrad);
  EXPECT_TRUE(ArraysNear(s3x.allGradient(), {-3,0,3}));
  EXPECT_TRUE(ArraysNear(s3y.allGradient(), {-3,0,3}));
  EXPECT_TRUE(ArraysNear(s3z.allGradient(), {-3,0,3}));
}


TEST(PhaseFieldTest, TestExternalForce) {
  PhaseField pot;

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
  auto stateForce1 = pot.setGridSize({2,2,2}).setSurfaceTension(0).setForce({-4,0,0}).newState({-1,-1,-1,-1, 1,1,1,1});
  EXPECT_FLOAT_EQ(stateForce1.energy(), 8);
  EXPECT_TRUE(ArraysNear(stateForce1.gradient(), {0,0,0,0, 1,1,1,1}, 1e-6));
}


TEST(PhaseFieldTest, TestSurfaceEnergy) {
  PhaseField pot;
  pot.setGridSize({2,1,1}).setSurfaceTension(0).setContactAngle({90,60}).setSolid({1,0});
  auto state = pot.newState({0.0, 0.5});
  EXPECT_FLOAT_EQ(state.energy(), 0.5/sqrt(2.0)*(-27.0/24));
  EXPECT_TRUE(ArraysNear(state.gradient(), {0, 0.5/sqrt(2.0)*(-0.75)}, 1e-6));
}


TEST(PhaseFieldTest, TestPressureConstraint) {
  PhaseField pot;
  pot.setGridSize({6,1,1}).setSolid({1,0,0,0,0,1}).setPressure({10});
  auto state = pot.newState({1,1,1,1,1,1});
  EXPECT_FLOAT_EQ(state.energy(), -30);
  EXPECT_TRUE(ArraysNear(state.gradient(), {0, -2.5, -5, -5, -2.5, 0}, 1e-6));
}


TEST(PhaseFieldTest, TestVolumeConstraint) {
  PhaseField pot;
  pot.setGridSize({6,1,1});
  pot.setSolid({1,0,0,0,0,1});

  pot.setVolumeFixed(true, 1);
  State state1(pot, {1,1,1,1,1,1});
  EXPECT_TRUE(static_cast<PhaseField&>(*state1.pot).volumeFixed);
  EXPECT_TRUE(ArraysNear(static_cast<PhaseField&>(*state1.pot).volume, {3}, 1e-6));
  EXPECT_FLOAT_EQ(state1.energy(), 0);

  pot.setVolume({1}, 100);
  State state2(pot, {1,1,1,1,1,1});
  EXPECT_FLOAT_EQ(state2.energy(), 100 * pow(2, 2));
  // EXPECT_TRUE(ArraysNear(state2.gradient(), {0, 100, 200, 200, 100, 0}));
  EXPECT_TRUE(ArraysNear(state2.gradient(), {0, 0, 0, 0, 0, 0})); // Zero gradient when not on interfaces
}


TEST(PhaseFieldTest, TestResolution) {
  PhaseField pot;
  EXPECT_FLOAT_EQ(pot.resolution, 1);

  pot.setResolution(2);
  ASSERT_FLOAT_EQ(pot.resolution, 2);

  // Bulk fluid
  State s1 = pot.setGridSize({5,1,1}).setSolid({1,0,0,0,1}).newState({0,-1,0,1,0}, {0});
  ASSERT_FLOAT_EQ(static_cast<PhaseField&>(*s1.pot).kappa[0], 1.5);
  ASSERT_FLOAT_EQ(static_cast<PhaseField&>(*s1.pot).kappaP[0], 6);
  double ebulk = 1.5/16 * 8;
  double egrad = 2 * (6.0/4) * 2;
  EXPECT_FLOAT_EQ(s1.allEnergy(), ebulk+egrad); // Bulk: 2, Gradient: 2
  EXPECT_TRUE(ArraysNear(s1.allGradient(), {0,-6,0,6,0}));

  // External Force
  pot.setGridSize({2,2,2}).setSurfaceTension(0).setSolid({0,0,0,0,0,0,0,0}).setForce({-4,0,0});
  auto s2 = pot.newState({-1,-1,-1,-1, 1,1,1,1});
  EXPECT_FLOAT_EQ(s2.energy(), 128);
  EXPECT_TRUE(ArraysNear(s2.gradient(), {0,0,0,0, 16,16,16,16}, 1e-6));
  pot.setForce({0,0,0});

  // Surface energy
  pot.setGridSize({2,1,1}).setSurfaceTension(0).setSolid({1,0}).setContactAngle({90,60});
  auto s3 = pot.newState({0, 0.5});
  EXPECT_FLOAT_EQ(s3.energy(), 0.5/sqrt(2.0)*(-27.0/24)*4);
  EXPECT_TRUE(ArraysNear(s3.gradient(), {0, 0.5/sqrt(2.0)*(-0.75)*4}, 1e-6));

  // Pressure
  pot.setGridSize({6,1,1}).setSolid({1,0,0,0,0,1}).setContactAngle({90,90,90,90,90,90}).setPressure({10});
  auto s4 = pot.newState({1,1,1,1,1,1});
  EXPECT_FLOAT_EQ(s4.energy(), -240);
  EXPECT_TRUE(ArraysNear(s4.gradient(), {0, -20, -40, -40, -20, 0}, 1e-6));
  pot.setPressure({0});

  // Volume
  pot.setGridSize({6,1,1}).setSurfaceTension(1).setSolid({1,0,0,0,0,1}).setContactAngle({90,90,90,90,90,90}).setVolume({8}, 100);
  auto s5 = pot.newState({1,1,1,1,1,1});
  EXPECT_FLOAT_EQ(s5.energy(), 100.0/16 * pow(2*8, 2));
  // EXPECT_TRUE(ArraysNear(s5.gradient(), {0, 64*100, 64*200, 64*200, 64*100, 0}));
  EXPECT_TRUE(ArraysNear(s5.gradient(), {0,0,0,0,0,0})); // Zero gradient unless at interface
}


TEST(PhaseFieldTest, TestNFluid) {
  PhaseField pot;
  EXPECT_FLOAT_EQ(pot.nFluid, 1);
  pot.setNFluid(3);
  EXPECT_FLOAT_EQ(pot.nFluid, 3);
}

TEST(PhaseFieldTest, TestFixFluid) {
  int nGrid = 4;
  PhaseField pot;
  pot.setNFluid(3).setGridSize({2,2,1});
  pot.setDensityConstraint("none");

  pot.init(vector<double>(3*nGrid));
  EXPECT_TRUE(ArraysMatch(pot.fixFluid, {false,false,false}));
  pot.setFixFluid(1);
  EXPECT_TRUE(ArraysMatch(pot.fixFluid, {false,true,false}));
  pot.setFixFluid(1, false);
  EXPECT_TRUE(ArraysMatch(pot.fixFluid, {false,false,false}));

  pot.setFixFluid(0);
  State s(pot, vector<double>(3*nGrid, 1));
  auto g = s.gradient();
  for (int iGrid=0; iGrid<nGrid; iGrid++) {
    EXPECT_FLOAT_EQ(g[3*iGrid], 0);
  }
}
