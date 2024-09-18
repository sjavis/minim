#include <iostream>
#include <vector>
#include "minim.h"

using namespace minim;

class Toy2d : public NewPotential<Toy2d> {
  public:
    double energy(const std::vector<double>& coords) const override {
      return coords[0]*coords[0] + coords[1]*coords[1] + coords[0]*coords[0]*coords[0]*coords[0] + coords[1]*coords[1]*coords[1]*coords[1];
    }

    std::vector<double> gradient(const std::vector<double>& coords) const override {
      std::vector<double> g(coords.size());
      g[0] = 2*coords[0] + 4*coords[0]*coords[0]*coords[0];
      g[1] = 2*coords[1] + 4*coords[1]*coords[1]*coords[1];
      return g;
    }
};


int main(int argc, char** argv) {
  mpi.init(&argc, &argv);
  print("Number of processors:", mpi.size, "; Rank:", mpi.rank);

  Toy2d potential = Toy2d();
  State state = potential.newState({1,4});
  state.convergence = 1e-4;

  //GradDescent min = GradDescent();
  //Lbfgs min = Lbfgs().setMaxIter(100);
  //Fire min = Fire();
  Anneal min = Anneal(100, 1);

  auto result = min.minimise(state);

  print(min.iter, result);
}
