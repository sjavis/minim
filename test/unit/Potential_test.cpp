#include "test_main.cpp"
#include "Potential.h"

#include "State.h"
#include "utils/vec.h"

using namespace minim;


TEST(PotentialTest, TestConstraints) {
  auto efunc = [](const vector<double>& x){ return vec::dotProduct(x, x); };
  auto gfunc = [](const vector<double>& x){ return 2*x; };
  Potential pot(efunc, gfunc);
  vector<double> coords{1, 2};
  std::unique_ptr<Communicator> comm = pot.newComm();

  // Single degrees of freedom
  pot.setConstraints({1});
  EXPECT_EQ(pot.constraints.size(), 1);
  EXPECT_TRUE(ArraysMatch(pot.constraints[0].idof, {1}));
  auto g = pot.gradient(coords);
  EXPECT_TRUE(ArraysMatch(g, {2,4}));
  pot.applyConstraints(coords, *comm, g);
  EXPECT_TRUE(ArraysMatch(g, {2,0}));

  // Vector constraint
  pot.constraints = {};
  pot.setConstraints({{0,1}}, {1,-1});
  EXPECT_TRUE(ArraysMatch(pot.constraints[0].idof, {0,1}));
  EXPECT_TRUE(ArraysMatch(pot.constraints[0].normal(coords), {1,-1}));
  g = pot.gradient(coords);
  EXPECT_TRUE(ArraysMatch(g, {2,4}));
  pot.applyConstraints(coords, *comm, g);
  EXPECT_TRUE(ArraysMatch(g, {3,3}));

  // Function constraint
  pot.constraints = {};
  pot.setConstraints({{0,1}}, [](const vector<int>& idof, const vector<double>& x){ return x-1; });
  EXPECT_TRUE(ArraysMatch(pot.constraints[0].idof, {0,1}));
  EXPECT_TRUE(ArraysMatch(pot.constraints[0].normal(coords), {0,1}));
  g = pot.gradient(coords);
  EXPECT_TRUE(ArraysMatch(g, {2,4}));
  pot.applyConstraints(coords, *comm, g);
  EXPECT_TRUE(ArraysMatch(g, {2,0}));

  // State apply constraints
  State s(pot, coords);
  g = s.blockGradient();
  s.applyConstraints(g);
  EXPECT_TRUE(ArraysMatch(g, {2,0}));

  // Check if fixed
  pot.constraints = {};
  pot.setConstraints({1});
  EXPECT_TRUE(pot.isFixed(1));
  EXPECT_FALSE(pot.isFixed(0));
  EXPECT_TRUE(ArraysMatch(pot.isFixed({0,1}), {false,true}));
}
