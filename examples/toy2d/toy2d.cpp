#include <iostream>
#include <vector>
#include "minim.h"

class Toy2d : public Potential {
  public:
    double energy(std::vector<double> coords, Args &args) override {
      return coords[0]*coords[0] + coords[1]*coords[1] + coords[0]*coords[0]*coords[0]*coords[0] + coords[1]*coords[1]*coords[1]*coords[1];
    }

    std::vector<double> gradient(std::vector<double> coords, Args &args) override {
      std::vector<double> g(args.ndof);
      g[0] = 2*coords[0] + 4*coords[0]*coords[0]*coords[0];
      g[1] = 2*coords[1] + 4*coords[1]*coords[1]*coords[1];
      return g;
    }
};

int main() {
  Toy2d potential = Toy2d();
  State state = potential.newState({1,4});
  state.convergence = 1e-4;

  //GradDescent min = GradDescent(state);
  Lbfgs min = Lbfgs(state);
  min.minimise();

  std::cout << min.iter << " " << state[0] << " " << state[1] << std::endl;
}
