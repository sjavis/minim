#include <iostream>
#include "Minimiser.h"
#include "GradDescent.h"
#include "System.h"

class Toy2d : public System {
  public:
    Toy2d(double* state) : System(2, state) {};
    ~Toy2d() {};

    double energy(double* state){
      return state[0]*state[0] + state[1]*state[1];
    }

    void gradient(double* state, double* g){
      g[0] = 2*state[0];
      g[1] = 2*state[1];
    }
};

int main() {
  double state[2] = {1,4};
  Toy2d sys = Toy2d(state);
  GradDescent min = GradDescent(sys, 1e-3);

  min.minimise();
  std::cout << sys.state[0] << " " << sys.state[1] << std::endl;
}
