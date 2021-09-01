#include <vector>
#include "minim.h"


int main(int argc, char **argv) {
  mpiInit(&argc, &argv);
  minim::print();

  Lj3d potential = Lj3d();
  State state = potential.newState({0,0,0, 2,0,0, 1,1,0, 5,1,1});
  state.convergence = 1e-4;

  //GradDescent min = GradDescent(state);
  Lbfgs min = Lbfgs(state);
  auto result = min.minimise();

  minim::print("Complete after", min.iter, "iterations.");
  minim::print(result[0], result[1], result[2]);
  minim::print(result[3], result[4], result[5]);
  minim::print(result[6], result[7], result[8]);
  minim::print(result[9], result[10], result[11]);
}
