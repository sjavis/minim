#include <cmath>
#include "Minimiser.h"
#include "State.h"
#include "vec.h"


Minimiser::Minimiser(State &state)
 : state(state), _adjustModel(NULL)
{}

Minimiser::Minimiser(State &state, AdjustFunc adjustModel)
 : state(state), _adjustModel(adjustModel)
{}


void Minimiser::setMaxIter(int maxIter) {
  this->maxIter = maxIter;
}


void Minimiser::minimise() {
  for (iter=0; iter<maxIter; iter++) {
    adjustModel();
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


void Minimiser::adjustModel() {
  if (_adjustModel) {
    _adjustModel(iter, state.coords, state.args);
  }
}
