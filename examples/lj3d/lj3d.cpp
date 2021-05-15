#include <iostream>
#include <vector>
#include "minim.h"


int main(int argc, char **argv) {
  mpiInit(&argc, &argv);
  std::cout << "Number of processors: " << minim::mpi.size << "; Rank: " << minim::mpi.rank << std::endl;

  Lj3d potential = Lj3d();
  State state = potential.newState({0,0,0,2,0,0});
  state.convergence = 1e-4;

  GradDescent min = GradDescent(state);
  min.minimise();

  std::cout << min.iter << std::endl;
  std::cout << state[0] << " " << state[1] << " " << state[2] << std::endl;
  std::cout << state[3] << " " << state[4] << " " << state[5] << std::endl;
}
