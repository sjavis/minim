#include <iostream>
#include <vector>
#include "minim.h"

class Toy2d : public System {
  public:
    Toy2d(std::vector<double> state) : System(2, state) {};
    ~Toy2d() {};

    double energy(std::vector<double> state){
      return state[0]*state[0] + state[1]*state[1] + state[0]*state[0]*state[0]*state[0] + state[1]*state[1]*state[1]*state[1];
    }

    void gradient(std::vector<double> state, std::vector<double> &g){
      g[0] = 2*state[0] + 4*state[0]*state[0]*state[0];
      g[1] = 2*state[1] + 4*state[1]*state[1]*state[1];
    }
};

int main() {
  Toy2d sys = Toy2d({1,4});
  sys.convergence_rms = 1e-4;
  //GradDescent min = GradDescent(sys, 1e-3);
  Lbfgs min = Lbfgs(sys);

  min.minimise();
  std::cout << min.iter << " " << sys.state[0] << " " << sys.state[1] << std::endl;
}
