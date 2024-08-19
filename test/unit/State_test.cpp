#include "test_main.cpp"
#include "State.h"

#include "potentials/LjNd.h"
#include "communicators/CommUnstructured.h"
#include "utils/mpi.h"

using namespace minim;
typedef std::vector<double> Vector;


TEST(StateTest, TestInit) {
  Lj3d pot;
  State s = pot.newState({0,0,0, 1,0,0});
  EXPECT_EQ(s.ndof, 6);
  EXPECT_TRUE(ArraysMatch(s.coords(), {0,0,0, 1,0,0}));
  EXPECT_NO_THROW(dynamic_cast<LjNd&>(*s.pot));
  EXPECT_NO_THROW(dynamic_cast<CommUnstructured&>(*s.comm));
}


TEST(StateTest, TestClone) {
  Lj3d pot;
  State s = pot.newState({0,0,0, 1,0,0});
  State sClone = s;
  EXPECT_EQ(sClone.ndof, 6);
  EXPECT_TRUE(ArraysMatch(sClone.coords(), {0,0,0, 1,0,0}));
  EXPECT_NO_THROW(dynamic_cast<LjNd&>(*sClone.pot));
  EXPECT_NO_THROW(dynamic_cast<CommUnstructured&>(*sClone.comm));
}


TEST(StateTest, TestEnergyGradient) {
  Lj3d pot;
  State s = pot.newState({0,0,0, 2,0,0});
  double e1 = s.energy();
  Vector g1 = s.gradient();
  double e2;
  Vector g2;
  s.energyGradient(&e2, &g2);
  EXPECT_FLOAT_EQ(e1, -63./1024);
  EXPECT_FLOAT_EQ(e2, -63./1024);
  EXPECT_TRUE(ArraysNear(g1, {-93./512,0,0, 93./512,0,0}, 1e-6));
  EXPECT_TRUE(ArraysNear(g2, {-93./512,0,0, 93./512,0,0}, 1e-6));

  Vector coords = {0,0,0, 2,0,0};
  double e3 = s.energy(coords);
  Vector g3 = s.gradient(coords);
  EXPECT_FLOAT_EQ(e3, -63./1024);
  EXPECT_TRUE(ArraysNear(g3, {-93./512,0,0, 93./512,0,0}, 1e-6));
}


TEST(StateTest, TestBlockEnergyGradient) {
  Lj3d pot;
  State s = pot.newState({0,0,0, 2,0,0});
  int nElements = (mpi.rank==0) ? 1 : 0;
  ASSERT_EQ(s.pot->elements.size(), nElements);

  double e1 = s.blockEnergy();
  Vector g1 = s.blockGradient();
  double e2;
  Vector g2;
  s.blockEnergyGradient(&e2, &g2);
  EXPECT_FLOAT_EQ(e1, (mpi.rank==0) ? -63./1024 : 0);
  EXPECT_FLOAT_EQ(e2, (mpi.rank==0) ? -63./1024 : 0);
  Vector gBlock = (mpi.rank==0) ? Vector({-93./512,0,0, 93./512,0,0}) : Vector();
  EXPECT_TRUE(ArraysNear(g1, gBlock, 1e-6));
  EXPECT_TRUE(ArraysNear(g2, gBlock, 1e-6));
}


TEST(StateTest, TestProcEnergyGradient) {
  // Note: This example does not have a halo region. It would be better to add an example with one.
  Lj3d pot;
  State s = pot.newState({0,0,0, 2,0,0});
  int nElements = (mpi.rank==0) ? 1 : 0;
  ASSERT_EQ(s.pot->elements.size(), nElements);

  double e1 = s.procEnergy();
  Vector g1 = s.procGradient();
  double e2;
  Vector g2;
  s.procEnergyGradient(&e2, &g2);
  EXPECT_FLOAT_EQ(e1, (mpi.rank==0) ? -63./1024 : 0);
  EXPECT_FLOAT_EQ(e2, (mpi.rank==0) ? -63./1024 : 0);
  Vector gProc = (mpi.rank==0) ? Vector({-93./512,0,0, 93./512,0,0}) : Vector();
  EXPECT_TRUE(ArraysNear(g1, gProc, 1e-6));
  EXPECT_TRUE(ArraysNear(g2, gProc, 1e-6));
}


TEST(StateTest, TestAssignProc) {
  Lj3d pot;
  State s1 = pot.newState({0,0,0, 2,0,0}, {0});
  State s2 = pot.newState({0,0,0, 2,0,0}, {1});
  EXPECT_EQ(s1.usesThisProc, (mpi.rank==0)?true:false);
  EXPECT_EQ(s2.usesThisProc, (mpi.rank==1)?true:false);
  EXPECT_TRUE(ArraysMatch(s1.comm->ranks, {0}));
  EXPECT_TRUE(ArraysMatch(s2.comm->ranks, {1}));
  double e1 = s1.energy();
  double e2 = s2.energy();
  Vector g1 = s1.gradient();
  Vector g2 = s2.gradient();
  EXPECT_FLOAT_EQ(e1, (mpi.rank==0)?-63./1024:0);
  EXPECT_FLOAT_EQ(e2, (mpi.rank==1)?-63./1024:0);
  EXPECT_TRUE(ArraysNear(g1, (mpi.rank==0)?Vector{-93./512,0,0, 93./512,0,0}:Vector(), 1e-6));
  EXPECT_TRUE(ArraysNear(g2, (mpi.rank==1)?Vector{-93./512,0,0, 93./512,0,0}:Vector(), 1e-6));
}


TEST(StateTest, TestAllEnergy) {
  Lj3d pot;
  State s1 = pot.newState({0,0,0, 2,0,0}, {0});
  State s2 = pot.newState({0,0,0, 2,0,0}, {1});
  double e1 = s1.allEnergy();
  double e2 = s2.allEnergy();
  Vector g1 = s1.allGradient();
  Vector g2 = s2.allGradient();
  EXPECT_FLOAT_EQ(e1, -63./1024);
  EXPECT_FLOAT_EQ(e2, -63./1024);
  EXPECT_TRUE(ArraysNear(g1, {-93./512,0,0, 93./512,0,0}, 1e-6));
  EXPECT_TRUE(ArraysNear(g2, {-93./512,0,0, 93./512,0,0}, 1e-6));
  double e12, e22;
  Vector g12, g22;
  s1.allEnergyGradient(&e12, &g12);
  s2.allEnergyGradient(&e22, &g22);
  EXPECT_FLOAT_EQ(e1, -63./1024);
  EXPECT_FLOAT_EQ(e2, -63./1024);
  EXPECT_TRUE(ArraysNear(g1, {-93./512,0,0, 93./512,0,0}, 1e-6));
  EXPECT_TRUE(ArraysNear(g2, {-93./512,0,0, 93./512,0,0}, 1e-6));
}
