#include "test_main.cpp"
#include "Potential.h"

#include "State.h"
#include "utils/vec.h"

using namespace minim;


TEST(PotentialTest, TestConstraints) {
  auto efunc = [](const vector<double>& x){ return vec::dotProduct(x, x); };
  auto gfunc = [](const vector<double>& x){ return 2*x; };
  Potential pot(efunc, gfunc);

  // Single degrees of freedom
  pot.setConstraints({1});
  EXPECT_EQ(pot.constraints.size(), 1);
  EXPECT_TRUE(ArraysMatch(pot.constraints[0].idof, {1}));

  // Vector constraint
  pot.constraints = {};
  pot.setConstraints({{0,1}}, {1,-1});
  EXPECT_TRUE(ArraysMatch(pot.constraints[0].idof, {0,1}));
  EXPECT_TRUE(ArraysMatch(pot.constraints[0].normal({1,2}), {1,-1}));

  // Function constraint
  pot.constraints = {};
  pot.setConstraints({{0,1}}, [](const vector<double>& x){ return x; });
  EXPECT_TRUE(ArraysMatch(pot.constraints[0].idof, {0,1}));
  EXPECT_TRUE(ArraysMatch(pot.constraints[0].normal({1,2}), {1,2}));
}
