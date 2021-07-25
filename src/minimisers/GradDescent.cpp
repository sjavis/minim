#include <math.h>
#include <vector>
#include "State.h"
#include "Minimiser.h"
#include "GradDescent.h"
#include "vec.h"


GradDescent::GradDescent(State &state)
: Minimiser(state), _g(state.ndof)
{}


GradDescent& GradDescent::setAlpha(double alpha) {
  _alpha = alpha;
  return *this;
}


GradDescent& GradDescent::setMaxIter(int maxIter) {
  Minimiser::setMaxIter(maxIter);
  return *this;
}


void GradDescent::iteration() {
  _g = state.gradient()
  auto step = vec::multiply(-_alpha, _g);
  state.blockCoords(vec::sum(state.blockCoords(), step));
}


bool GradDescent::checkConvergence() {
  double sum = vec::dotProduct(_g, _g);
  double rms = sqrt(sum/state.ndof);
  return (rms < state.convergence);
}
