#include <iostream>
#include <vector>
#include "minim.h"


int main(int argc, char **argv) {
  mpiInit(&argc, &argv);
  std::cout << "Number of processors: " << minim::mpi.size << "; Rank: " << minim::mpi.rank << std::endl;

  Lj3d potential = Lj3d();
  State state = potential.newState({0,0,0, 2,0,0, 1,1,0});
  state.convergence = 1e-4;

  GradDescent min = GradDescent(state);
  auto result = min.minimise();

  if (minim::mpi.rank == 0) {
    std::cout << min.iter << std::endl;
    std::cout << result[0] << " " << result[1] << " " << result[2] << std::endl;
    std::cout << result[3] << " " << result[4] << " " << result[5] << std::endl;
    std::cout << result[6] << " " << result[7] << " " << result[8] << std::endl;
  }
}
