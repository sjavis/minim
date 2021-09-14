#include <iostream>
#include <vector>
#include "minim.h"


class Toy2d : public Potential {
  public:
    double energy(const std::vector<double> &coords, const Args &args) override {
      return coords[0]*coords[0] + coords[1]*coords[1] + coords[0]*coords[0]*coords[0]*coords[0] + coords[1]*coords[1]*coords[1]*coords[1];
    }

    std::vector<double> gradient(const std::vector<double> &coords, const Args &args) override {
      std::vector<double> g(args.ndof);
      g[0] = 2*coords[0] + 4*coords[0]*coords[0]*coords[0];
      g[1] = 2*coords[1] + 4*coords[1]*coords[1]*coords[1];
      return g;
    }
};


int main(int argc, char **argv) {
  mpiInit(&argc, &argv);
  std::cout << "Number of processors: " << minim::mpi.size << "; Rank: " << minim::mpi.rank << std::endl;

  Toy2d potential = Toy2d();
  State state = potential.newState({1,4});
  state.convergence = 1e-4;

  //GradDescent min = GradDescent(state);
  //Lbfgs min = Lbfgs(state).setMaxIter(100);
  //Fire min = Fire(state);
  Anneal min = Anneal(state, 100, 1);

  auto result = min.minimise();

  std::cout << min.iter << " " << result[0] << " " << result[1] << std::endl;
}
