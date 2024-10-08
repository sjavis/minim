#include "test_main.cpp"
#include "minimisers/Fire.h"

#include <math.h>
#include "State.h"
#include "Potential.h"
#include "utils/vec.h"
#include "utils/mpi.h"

using namespace minim;
typedef std::vector<double> Vector;

class Toy2d : public NewPotential<Toy2d> {
  public:
    double energy(const std::vector<double>& coords) const override {
      return coords[0]*coords[0] + coords[1]*coords[1];
    }

    std::vector<double> gradient(const std::vector<double>& coords) const override {
      return 2 * coords;
    }
};


TEST(FireTest, TestStep) {
  Vector coords = {1, 4};
  Vector g = {2, 8};

  Toy2d pot;
  State state = pot.newState(coords);
  ASSERT_FLOAT_EQ(state.energy(), 17);
  ASSERT_TRUE(ArraysNear(state.gradient(), g, 1e-6));

  Fire min = Fire().setMaxIter(0);
  min.minimise(state);

  double dt = 0.1 / pow(68, 0.25);
  Vector step = -0.5 * dt * dt * g;
  EXPECT_TRUE(ArraysMatch(state.coords(), coords+step));
}


TEST(FireTest, TestConvergence) {
  Toy2d pot;
  State state = pot.newState({1, 4});
  state.convergence = 2e-6/sqrt(2); // Designed to make a radius of convergence of 1e-6

  Fire min = Fire();
  min.minimise(state);

  EXPECT_LT(vec::norm(state.coords()), 1e-6);
}


TEST(FireTest, TestSingleProc) {
  Toy2d pot;
  State state = pot.newState({1, 4}, {0});
  state.convergence = 2e-6/sqrt(2); // Designed to make a radius of convergence of 1e-6
  Fire min = Fire();
  auto result = min.minimise(state);
  if (mpi.rank == 0) {
    EXPECT_LT(vec::norm(result), 1e-6);
    EXPECT_LT(vec::norm(state.coords()), 1e-6);
  } else {
    EXPECT_TRUE(result.empty());
  }
}
