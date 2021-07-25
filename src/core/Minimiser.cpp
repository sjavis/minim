#include <math.h>
#include "Minimiser.h"
#include "State.h"
#include "vec.h"


Minimiser::Minimiser(State &state) : state(state), _adjustModel(NULL) {}
Minimiser::Minimiser(State &state, AdjustFunc adjustModel) : state(state), _adjustModel(adjustModel) {}


void Minimiser::setMaxIter(int maxIter) {
  this->maxIter = maxIter;
}


std::vector<double> Minimiser::minimise() {
  for (iter=0; iter<maxIter; iter++) {
    adjustModel();
    iteration();
    if (checkConvergence()) break;
  }
  return state.getCoords();
}


void Minimiser::adjustModel() {
  if (_adjustModel) {
    _adjustModel(iter, state);
  }
}
