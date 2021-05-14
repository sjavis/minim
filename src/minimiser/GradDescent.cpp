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
  auto step = vec::multiply(-_alpha, state.gradient());
  state.blockCoords(vec::sum(state.blockCoords(), step));
}
