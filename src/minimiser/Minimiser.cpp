#include <cmath>
#include "Minimiser.h"
#include "State.h"
#include "vec.h"


Minimiser::Minimiser(State &state, int maxIter)
 : state(state), maxIter(maxIter)
{}


void Minimiser::minimise() {
  for (iter=0; iter<maxIter; iter++) {
    iteration();
    if (checkConvergence()) break;
  }
}


bool Minimiser::checkConvergence() {
  std::vector<double> g = state.gradient();
  double sum = vec::dotProduct(g, g);
  double rms = sqrt(sum/state.ndof);
  return (rms < state.convergence);
}
